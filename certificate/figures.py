"""Figures for the manuscript. Run after exps.py (needs results/trace.json)."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exps

R = exps.OUT
FIG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "figures"); os.makedirs(FIG, exist_ok=True)
plt.rcParams.update({"font.family": "serif", "font.size": 10, "axes.linewidth": 0.8,
                     "mathtext.fontset": "cm"})
COL = {"tangent": "#1f4e79", "outward": "#2e8b57", "normal": "#7f7f7f",
       "mixed": "#b8860b", "inward": "#b22222"}

# ---------------- Figure 1: certificate envelope (raw trials, same generator as exp_cert)
n, p = 100, 10
pts = {k: [] for k in COL}
for kind in COL:
    for eps in [1e-1, 3e-2, 1e-2, 1e-3, 1e-4, 1e-6]:
        for s in range(30):
            rng = np.random.default_rng(1000 * s + 7)
            Q = exps.random_frame(rng, n, p)
            E, _, _ = exps.perturbation(rng, Q, kind, eps)
            A = Q @ Q.T + E
            rho = exps.residual(A)
            if rho >= 0.25:
                continue
            d, _ = exps.dist_to_projectors(A)
            pts[kind].append((rho, d / rho))
np.savez(os.path.join(R, "fig1_points.npz"), **{k: np.array(v) for k, v in pts.items()})

fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
r = np.logspace(-14, np.log10(0.2499), 400)
ax[0].plot(r, exps.phi(r) / r, "k-", lw=1.2, label=r"$\varphi(\rho)/\rho$ (upper)")
ax[0].plot(r, exps.psi(r) / r, "k--", lw=1.2, label=r"$\psi(\rho)/\rho$ (lower)")
for k, v in pts.items():
    v = np.array(v)
    ax[0].plot(v[:, 0], v[:, 1], "o", ms=2.5, color=COL[k], alpha=0.7, label=k, mec="none")
ax[0].set_xscale("log"); ax[0].set_xlabel(r"residual $\rho(A)=\|A^2-A\|_F$")
ax[0].set_ylabel(r"$\delta(A)/\rho(A)$")
ax[0].set_ylim(0.75, 1.35); ax[0].set_xlim(1e-14, 0.3)
ax[0].legend(fontsize=7, frameon=False, loc="upper left", ncol=2)
ax[0].set_title("(a) distance / residual, 5 perturbation types", fontsize=9)

r2 = np.linspace(1e-4, 0.25, 400)
ax[1].plot(r2, exps.phi(r2) / r2, "k-", lw=1.4, label=r"$\varphi(\rho)/\rho$ (Thm 3.2, sharp)")
ax[1].plot(r2, 1 + 4 * r2, "-", color="#b8860b", lw=1.2, label=r"$1+4\rho$ (Tang)")
ax[1].plot(r2, 2 + 0 * r2, "-", color="#b22222", lw=1.2, label=r"$2$ (Weyl-type)")
ax[1].plot(r2, exps.psi(r2) / r2, "k--", lw=1.4, label=r"$\psi(\rho)/\rho$ (Thm 3.2, sharp)")
ax[1].plot(r2, 1 / (1 + r2), ":", color="#1f4e79", lw=1.4, label=r"$1/(1+\rho)$")
ax[1].set_xlabel(r"residual $\rho$"); ax[1].set_ylabel("bound / $\\rho$")
ax[1].set_ylim(0.6, 2.15)
ax[1].legend(fontsize=7, frameon=False, loc="upper left")
ax[1].set_title("(b) upper and lower certificates compared", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig1_certificate.png"), dpi=300)

# ---------------- Figure 2: certified lazy re-projection traces
tr = json.load(open(os.path.join(R, "trace.json")))
fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.0))
sty = {"apriori_noshift": ("A priori test, no shift", "#b22222"),
       "certified_noshift": ("Certified test, no shift", "#b8860b"),
       "certified_shift": ("Certified test, with shift", "#1f4e79")}
for k, (lab, c) in sty.items():
    a = np.array(tr[k])
    ax[0].semilogy(a[:, 0], np.maximum(a[:, 1], 1e-16), "-", color=c, lw=1.0, label=lab)
    ax[1].plot(a[:, 0], a[:, 3], "-", color=c, lw=1.2, label=lab)
a = np.array(tr["certified_shift"])
ax[0].semilogy(a[:, 0], a[:, 2], "-", color="0.6", lw=0.8, label=r"single-step bound $\varepsilon_k^2/\sqrt{2}$")
ax[0].axhline(1e-10, color="k", ls="--", lw=0.8, label=r"tolerance $\tau$")
ax[0].set_xlabel("iteration $k$"); ax[0].set_ylabel(r"$\delta(P_k)$ (distance to $\Pi_p$)")
ax[0].set_ylim(1e-16, 1e0)
ax[0].legend(fontsize=6.5, frameon=True, framealpha=0.92, edgecolor="none", loc="upper right")
ax[0].set_title("(a) distance of the iterate to the manifold", fontsize=9)
ax[1].plot(a[:, 0], a[:, 0] + 1, ":", color="0.4", lw=1.0, label="always re-project")
ax[1].set_xlabel("iteration $k$"); ax[1].set_ylabel("cumulative eigendecompositions")
ax[1].legend(fontsize=7, frameon=False, loc="upper left")
ax[1].set_title("(b) cost", fontsize=9)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "fig2_lazy.png"), dpi=300)
print("figures written")
