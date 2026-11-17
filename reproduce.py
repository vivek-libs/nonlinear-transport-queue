#!/usr/bin/env python3
"""
Reproduction script for:
  "A Reproducible Finite-Volume Coupling of a
   Nonlinear Transport Equation to a Periodic Queue"

Produces:
  - baseline peak queues and conservation residuals (linear vs nonlinear)
  - matched-flux control: linear reference rerun with u_in scaled so that
    f_lin(u_in_lin) == f_nl(0.06); isolates the plume-compression mechanism
    from the boundary-flux confound
  - mesh-refinement table rows for N in {100, 250, 500, 1000}
  - observed order of accuracy and Richardson-extrapolated asymptotic gap
  - parameter-sweep table rows over plume amplitude A and duty ratio W/P
  - splitting-error check against a 1000x finer queue reference
  - buffer-sizing exceedance fractions for a range of B
  - results.pdf (two-panel figure)

All numbers are deterministic. No random seed is used.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------- baseline parameters ---------------------------
L, v, K = 1.0, 1.0, 1.0
P, W, C = 2.0, 0.7, 0.22
T = 24.0
N_BASE = 250
CFL = 0.35
U_IN = 0.06
PLUME_AMP = 0.52
PLUME_CENTER = 0.24
PLUME_WIDTH = 0.085

def f_nonlinear(u): return v * u * (1.0 - u / K)
def f_linear(u):    return v * u

def rusanov_flux(a, b, f, speed_bound):
    return 0.5 * (f(a) + f(b)) - 0.5 * speed_bound * (b - a)

def service_availability(t0, t1, P, W, C):
    total = 0.0
    k0 = int(np.floor(t0 / P))
    k1 = int(np.floor(t1 / P))
    for k in range(k0, k1 + 1):
        a = max(t0, k * P)
        b = min(t1, k * P + W)
        if b > a:
            total += C * (b - a)
    return total

def run_coupled(f, speed_bound, u0, N, T, P, W, C,
                cfl=CFL, u_in=U_IN, q_substeps=1):
    dx = L / N
    dt = cfl * dx / speed_bound
    nsteps = int(np.ceil(T / dt))
    dt = T / nsteps
    u = u0.copy(); q = 0.0; t = 0.0
    q_hist = np.empty(nsteps + 1); q_hist[0] = q
    flux_in_hist = np.empty(nsteps); flux_out_hist = np.empty(nsteps)
    mass_hist = np.empty(nsteps + 1); mass_hist[0] = dx * u.sum()

    for n in range(nsteps):
        u_ext = np.concatenate(([u_in], u, [u[-1]]))
        F = rusanov_flux(u_ext[:-1], u_ext[1:], f, speed_bound)
        flux_in_hist[n] = F[0]
        flux_out_hist[n] = f(u[-1])
        u = u - (dt / dx) * (F[1:] - F[:-1])

        a = flux_out_hist[n] * dt
        if q_substeps == 1:
            A = service_availability(t, t + dt, P, W, C)
            d = min(q + a, A); q = q + a - d
        else:
            sub_dt = dt / q_substeps
            remaining = a
            for s in range(q_substeps):
                ts = t + s * sub_dt; te = ts + sub_dt
                A_s = service_availability(ts, te, P, W, C)
                d_s = min(q + remaining, A_s)
                q = q + remaining - d_s; remaining = 0.0
        t += dt
        q_hist[n + 1] = q
        mass_hist[n + 1] = dx * u.sum()

    return dict(u_final=u, q_hist=q_hist, flux_in_hist=flux_in_hist,
                flux_out_hist=flux_out_hist, mass_hist=mass_hist, dt=dt)

def conservation_residual(run):
    M = run["mass_hist"] + run["q_hist"]
    a = run["flux_out_hist"] * run["dt"]
    d = run["q_hist"][:-1] + a - run["q_hist"][1:]
    increments = run["dt"] * run["flux_in_hist"] - d
    return M[1:] - M[0] - np.cumsum(increments)

def make_u0(N, amp=PLUME_AMP, center=PLUME_CENTER, width=PLUME_WIDTH):
    x = (np.arange(N) + 0.5) * (L / N)
    return U_IN + amp * np.exp(-((x - center) / width) ** 2)

def baseline():
    u0 = make_u0(N_BASE)
    rn = run_coupled(f_nonlinear, v, u0, N_BASE, T, P, W, C, u_in=U_IN)
    rl = run_coupled(f_linear,    v, u0, N_BASE, T, P, W, C, u_in=U_IN)
    return dict(rn=rn, rl=rl,
                peak_nl=float(np.max(rn["q_hist"])),
                peak_lin=float(np.max(rl["q_hist"])),
                u0=u0)

def matched_flux_control():
    """Linear reference rerun with u_in scaled so that
       f_lin(u_in_lin) == f_nl(U_IN).  Isolates plume compression
       from the boundary-flux difference."""
    u_in_lin = f_nonlinear(U_IN) / v          # = 0.0564
    u0 = make_u0(N_BASE)
    rn = run_coupled(f_nonlinear, v, u0, N_BASE, T, P, W, C, u_in=U_IN)
    rl = run_coupled(f_linear,    v, u0, N_BASE, T, P, W, C, u_in=u_in_lin)
    pn = float(np.max(rn["q_hist"]))
    pl = float(np.max(rl["q_hist"]))
    return dict(u_in_lin=u_in_lin,
                peak_nl=pn, peak_lin_mf=pl,
                gap_pct=100*(pn-pl)/pl)

def mesh_refinement():
    rows = []
    for N in (100, 250, 500, 1000):
        u0 = make_u0(N)
        rn = run_coupled(f_nonlinear, v, u0, N, T, P, W, C, u_in=U_IN)
        rl = run_coupled(f_linear,    v, u0, N, T, P, W, C, u_in=U_IN)
        rows.append(dict(
            N=N,
            peak_nl=float(np.max(rn["q_hist"])),
            peak_lin=float(np.max(rl["q_hist"])),
            res_nl=float(np.max(np.abs(conservation_residual(rn)))),
            res_lin=float(np.max(np.abs(conservation_residual(rl)))),
        ))
    return rows

def observed_orders(rows):
    Ns = [r["N"] for r in rows]
    nl  = {r["N"]: r["peak_nl"] for r in rows}
    lin = {r["N"]: r["peak_lin"] for r in rows}
    out = {}
    for key, d in (("nl", nl), ("lin", lin)):
        num = d[500] - d[250]; den = d[1000] - d[500]
        out[key] = np.log2(abs(num/den)) if abs(den) > 0 else float("nan")
    g = {N: nl[N] - lin[N] for N in Ns}
    num = g[500] - g[250]; den = g[1000] - g[500]
    out["gap"] = np.log2(abs(num/den)) if abs(den) > 0 else float("nan")
    return out

def richardson(rows):
    nl  = {r["N"]: r["peak_nl"] for r in rows}
    lin = {r["N"]: r["peak_lin"] for r in rows}
    g   = {N: nl[N] - lin[N] for N in nl}
    p_nl  = np.log2(abs((nl[500]  - nl[250])  / (nl[1000]  - nl[500])))
    p_lin = np.log2(abs((lin[500] - lin[250]) / (lin[1000] - lin[500])))
    p_g   = np.log2(abs((g[500]   - g[250])   / (g[1000]   - g[500])))
    nl_inf  = nl[1000]  + (nl[1000]  - nl[500])  / (2**p_nl  - 1)
    lin_inf = lin[1000] + (lin[1000] - lin[500]) / (2**p_lin - 1)
    g_inf   = g[1000]   + (g[1000]   - g[500])   / (2**p_g   - 1)
    # robustness range over p in [1, 2]
    def g_at_p(p):
        return (g[1000] + (g[1000] - g[500])/(2**p - 1))
    lo = g_at_p(2.0); hi = g_at_p(1.0)
    lin_at_p1 = lin[1000] + (lin[1000] - lin[500])  # p = 1
    return dict(p_nl=p_nl, p_lin=p_lin, p_g=p_g,
                nl_inf=nl_inf, lin_inf=lin_inf,
                g_inf=g_inf, gap_pct=100*g_inf/lin_inf,
                gap_lo_pct=100*lo/lin_at_p1,
                gap_hi_pct=100*hi/lin_at_p1)

def parameter_sweep():
    rows = []
    for amp in (0.26, 0.52, 0.78):
        for duty in (0.35, 0.50, 0.70):
            W_sweep = duty * P
            u0 = make_u0(N_BASE, amp=amp)
            rn = run_coupled(f_nonlinear, v, u0, N_BASE, T, P, W_sweep, C, u_in=U_IN)
            rl = run_coupled(f_linear,    v, u0, N_BASE, T, P, W_sweep, C, u_in=U_IN)
            pnl = float(np.max(rn["q_hist"]))
            pln = float(np.max(rl["q_hist"]))
            rows.append(dict(amp=amp, duty=duty, peak_nl=pnl, peak_lin=pln,
                             ordering="nl>lin" if pnl > pln else "lin>nl"))
    return rows

def splitting_error():
    u0 = make_u0(N_BASE)
    run = run_coupled(f_nonlinear, v, u0, N_BASE, T, P, W, C, u_in=U_IN, q_substeps=1)
    dt = run["dt"]; flux_out = run["flux_out_hist"]
    q_step = run["q_hist"]
    n_fine = 1000; sub_dt = dt / n_fine
    q = 0.0; q_ref_peak = 0.0
    for n, a_flux in enumerate(flux_out):
        a = a_flux * dt
        for s in range(n_fine):
            ts = n * dt + s * sub_dt; te = ts + sub_dt
            A_s = service_availability(ts, te, P, W, C)
            d_s = min(q + a, A_s); q = q + a - d_s; a = 0.0
            if q > q_ref_peak: q_ref_peak = q
    return (float(np.max(q_step)), float(q_ref_peak),
            float(abs(np.max(q_step) - q_ref_peak)))

def buffer_exceedance(run, Bs):
    q = run["q_hist"]
    return {B: float(np.mean(q > B)) for B in Bs}

def make_figure(base, res_nl, res_lin, path="results.pdf"):
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6.4, 6.4), sharex=True,
                                   constrained_layout=True)
    dt = base["rn"]["dt"]
    t = np.arange(len(base["rn"]["q_hist"])) * dt
    ax1.plot(t, base["rl"]["q_hist"], label="linear flux")
    ax1.plot(t, base["rn"]["q_hist"], label="nonlinear flux")
    ax1.set_ylabel("Downstream queue (normalized units)")
    ax1.legend(loc="upper right"); ax1.grid(True, alpha=0.3)
    ax2.plot(t[1:], res_lin, label="linear flux")
    ax2.plot(t[1:], res_nl,  label="nonlinear flux")
    ax2.set_xlabel("Time (normalized units)")
    # Fix: Removed the hardcoded \times 10^{-13} so matplotlib's
    # automatic offset (1e-14) handles the scaling cleanly.
    ax2.set_ylabel("Mass-balance residual")
    ax2.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))
    ax2.grid(True, alpha=0.3)
    fig.savefig(path); plt.close(fig)

if __name__ == "__main__":
    print("=" * 72); print("BASELINE"); print("=" * 72)
    base = baseline()
    r_nl  = conservation_residual(base["rn"])
    r_lin = conservation_residual(base["rl"])
    print(f"peak nonlinear     = {base['peak_nl']:.6f}")
    print(f"peak linear        = {base['peak_lin']:.6f}")
    print(f"gap %              = {100*(base['peak_nl']-base['peak_lin'])/base['peak_lin']:.2f}")
    print(f"residual nonlinear = {np.max(np.abs(r_nl)):.3e}")
    print(f"residual linear    = {np.max(np.abs(r_lin)):.3e}")
    print(f"dt = {base['rn']['dt']:.3e}, steps = {len(base['rn']['q_hist']) - 1}")

    print("\n" + "=" * 72); print("MATCHED-FLUX CONTROL"); print("=" * 72)
    mc = matched_flux_control()
    print(f"u_in_lin (matched) = {mc['u_in_lin']:.4f}")
    print(f"peak nonlinear     = {mc['peak_nl']:.6f}")
    print(f"peak linear (match)= {mc['peak_lin_mf']:.6f}")
    print(f"matched-flux gap % = {mc['gap_pct']:.2f}  (> 5.3 expected)")

    print("\n" + "=" * 72); print("MESH REFINEMENT (Table 1)"); print("=" * 72)
    rows = mesh_refinement()
    for r in rows:
        print(f"{r['N']:4d} & {r['peak_nl']:.6f} & {r['peak_lin']:.6f} & "
              f"${r['res_nl']:.2e}$ & ${r['res_lin']:.2e}$ \\\\")

    print("\n" + "=" * 72); print("OBSERVED ORDER AND RICHARDSON"); print("=" * 72)
    oo = observed_orders(rows)
    print(f"observed p (nonlinear peak) = {oo['nl']:.3f}")
    print(f"observed p (linear peak)    = {oo['lin']:.3f}")
    print(f"observed p (gap)            = {oo['gap']:.3f}")
    ri = richardson(rows)
    print(f"Richardson nl_inf  = {ri['nl_inf']:.6f}")
    print(f"Richardson lin_inf = {ri['lin_inf']:.6f}")
    print(f"Richardson g_inf   = {ri['g_inf']:.6f}  ({ri['gap_pct']:.2f}%)")
    print(f"gap robustness over p in [1,2]: "
          f"[{min(ri['gap_lo_pct'], ri['gap_hi_pct']):.2f}%, "
          f"{max(ri['gap_lo_pct'], ri['gap_hi_pct']):.2f}%]")

    print("\n" + "=" * 72); print("PARAMETER SWEEP (Table 2)"); print("=" * 72)
    for r in parameter_sweep():
        print(f"A={r['amp']:.2f} W/P={r['duty']:.2f} : "
              f"nl={r['peak_nl']:.6f} lin={r['peak_lin']:.6f} [{r['ordering']}]")

    print("\n" + "=" * 72); print("SPLITTING ERROR"); print("=" * 72)
    p1, pref, diff = splitting_error()
    print(f"step-aggregated peak = {p1:.6f}")
    print(f"fine-grid peak       = {pref:.6f}")
    print(f"absolute difference  = {diff:.3e}")

    print("\n" + "=" * 72); print("BUFFER EXCEEDANCE"); print("=" * 72)
    Bs = (0.02, 0.05, 0.10, 0.15, 0.152)
    for label, run in (("nonlinear", base["rn"]), ("linear", base["rl"])):
        ex = buffer_exceedance(run, Bs)
        for B in Bs:
            print(f"  {label:9s} B={B:.3f}  exceedance={100*ex[B]:5.2f}%")

    print("\n" + "=" * 72); print("WRITING results.pdf"); print("=" * 72)
    make_figure(base, r_nl, r_lin)
    print("done.")