"""
Numerical experiments for quadratic perturbation bounds of orthogonal
projectors on the Grassmann manifold.

Reproduces:
  1. Residual scaling: tangent perturbations give an O(eps^2) idempotency
     defect, generic symmetric perturbations give O(eps).
  2. Rank dependence of the defect at fixed eps.
  3. (optional, --cifar) The same test around a PCA projector built from
     CIFAR-10 images.

Run from the repository root:
    python python/experiments.py            # synthetic experiments
    python python/experiments.py --cifar    # also the CIFAR-10 check

All outputs (CSV tables and PNG figures) are written to ./results/.
"""

import argparse
import os

import matplotlib

matplotlib.use("Agg")  # write files only; works without a display
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from projector_utils import (
    backward_error,
    nontangent_perturbation,
    random_projector,
    residual,
    tangent_perturbation,
)

SEED = 0
RESULTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")


def experiment(n=100, p=10, epsilons=np.logspace(-6, -1, 6), trials=30):
    """Average residuals and backward errors over random projectors."""
    rows = []
    for eps in epsilons:
        rt, et, rn, en = [], [], [], []
        for _ in range(trials):
            P = random_projector(n, p)
            Pt = P + tangent_perturbation(P, p, n, eps)
            Pn = P + nontangent_perturbation(n, eps)
            rt.append(residual(Pt))
            et.append(backward_error(Pt, p))
            rn.append(residual(Pn))
            en.append(backward_error(Pn, p))
        rows.append({
            "epsilon": eps,
            "tangent_residual": np.mean(rt),
            "tangent_backward_error": np.mean(et),
            "nontangent_residual": np.mean(rn),
            "nontangent_backward_error": np.mean(en),
            "tangent_residual_over_eps2": np.mean(rt) / eps**2,
            "nontangent_residual_over_eps": np.mean(rn) / eps,
        })
    return pd.DataFrame(rows)


def rank_dependence(n=100, p_values=(1, 5, 10, 20, 40), eps=1e-6, trials=30):
    """Residual / eps^2 as a function of the rank p (tangent perturbations)."""
    rows = []
    for p in p_values:
        r, e = [], []
        for _ in range(trials):
            P = random_projector(n, p)
            Pt = P + tangent_perturbation(P, p, n, eps)
            r.append(residual(Pt))
            e.append(backward_error(Pt, p))
        rows.append({
            "p": p,
            "residual": np.mean(r),
            "backward_error": np.mean(e),
            "residual_over_eps2": np.mean(r) / eps**2,
        })
    return pd.DataFrame(rows)


def cifar10_projector(p=10, samples=8000):
    """Rank-p PCA projector of grayscale CIFAR-10 images (n = 1024).

    Raises an error if the data cannot be downloaded; no silent fallback.
    """
    from sklearn.datasets import fetch_openml
    from sklearn.decomposition import PCA

    X, _ = fetch_openml("CIFAR_10", version=1, return_X_y=True,
                        as_frame=False, parser="auto")
    X = np.asarray(X[:samples], dtype=float)
    # OpenML stores each image as 3072 values: 1024 red, 1024 green, 1024 blue.
    X = X.reshape(-1, 3, 1024).mean(axis=1)  # grayscale, 1024 pixels
    X -= X.mean(axis=0)
    Q = PCA(n_components=p).fit(X).components_.T
    return Q @ Q.T


def cifar10_experiment(epsilons=np.logspace(-6, -2, 3), p=10, trials=50):
    P = cifar10_projector(p)
    n = P.shape[0]
    rows = []
    for eps in epsilons:
        rt, et = [], []
        for _ in range(trials):
            Pt = P + tangent_perturbation(P, p, n, eps)
            rt.append(residual(Pt))
            et.append(backward_error(Pt, p))
        rows.append({
            "epsilon": eps,
            "tangent_residual": np.mean(rt),
            "tangent_backward_error": np.mean(et),
            "tangent_residual_over_eps2": np.mean(rt) / eps**2,
        })
    return pd.DataFrame(rows)


def plot_residual_scaling(df, path):
    plt.figure(figsize=(7, 5))
    plt.loglog(df.epsilon, df.tangent_residual, "o-", label="Tangent perturbation")
    plt.loglog(df.epsilon, df.nontangent_residual, "s-", label="Generic symmetric perturbation")
    plt.loglog(df.epsilon, df.epsilon**2 / np.sqrt(2), "k--", label=r"$\varepsilon^2/\sqrt{2}$ (bound)")
    plt.loglog(df.epsilon, df.epsilon, "r:", label=r"$\varepsilon$ reference")
    plt.xlabel(r"Perturbation norm $\varepsilon = \|E\|_F$")
    plt.ylabel(r"Idempotency defect $\|\tilde P^2 - \tilde P\|_F$")
    plt.title("Quadratic vs linear scaling (n = 100, p = 10)")
    plt.legend()
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def plot_rank_dependence(df, path):
    p = np.array(df.p, dtype=float)
    plt.figure(figsize=(7, 5))
    plt.plot(p, df.residual_over_eps2, "o-", label="Observed (mean over trials)")
    plt.plot(p, np.full_like(p, 1 / np.sqrt(2)), "k--", label=r"Upper bound $1/\sqrt{2}$")
    plt.plot(p, 1 / np.sqrt(2 * p), "g:", label=r"Lower bound $1/\sqrt{2r}$, $r = p$")
    plt.xlabel("Rank p")
    plt.ylabel(r"Defect / $\varepsilon^2$")
    plt.title(r"Rank dependence ($\varepsilon = 10^{-6}$, n = 100)")
    plt.legend()
    plt.grid(True, alpha=0.5)
    plt.tight_layout()
    plt.savefig(path, dpi=200)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cifar", action="store_true",
                        help="also run the CIFAR-10 check (downloads data from OpenML)")
    args = parser.parse_args()

    np.random.seed(SEED)
    os.makedirs(RESULTS_DIR, exist_ok=True)
    fmt = "%.3e"

    df = experiment()
    df.to_csv(os.path.join(RESULTS_DIR, "residual_scaling.csv"), index=False)
    plot_residual_scaling(df, os.path.join(RESULTS_DIR, "residual_scaling.png"))
    print("Residual scaling (n=100, p=10, 30 trials)")
    print(df.to_string(index=False, float_format=fmt), "\n")

    rk = rank_dependence()
    rk.to_csv(os.path.join(RESULTS_DIR, "rank_dependence.csv"), index=False)
    plot_rank_dependence(rk, os.path.join(RESULTS_DIR, "rank_dependence.png"))
    print("Rank dependence (eps=1e-6)")
    print(rk.to_string(index=False, float_format=fmt), "\n")

    if args.cifar:
        cf = cifar10_experiment()
        cf.to_csv(os.path.join(RESULTS_DIR, "cifar10_check.csv"), index=False)
        print("CIFAR-10 PCA projector (n=1024, p=10)")
        print(cf.to_string(index=False, float_format=fmt), "\n")

    print(f"Results written to {os.path.abspath(RESULTS_DIR)}")


if __name__ == "__main__":
    main()
