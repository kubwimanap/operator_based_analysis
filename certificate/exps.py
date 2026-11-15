"""
Reproducible experiments for
"Idempotency residuals as sharp certificates of distance to the Grassmann manifold"
(Kubwimana, Reddy, Mkenyeleye).

Usage:  python exps.py <experiment> [--seeds 30]
Experiments: cert, sharp, sens, floor, computed, lazy_counts, lazy_time, all
Each experiment writes results/<name>.json. All randomness is seeded.
"""
import argparse, json, time, os, sys
import numpy as np

U64 = 2.0 ** -53            # unit roundoff, IEEE double
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "results")
os.makedirs(OUT, exist_ok=True)


# ---------------------------------------------------------------- certificate
def phi(r):   # upper certificate (inward eigenvalues), valid for r < 1/4
    return 2.0 * r / (1.0 + np.sqrt(1.0 - 4.0 * r))


def psi(r):   # lower certificate (outward eigenvalues)
    return 2.0 * r / (1.0 + np.sqrt(1.0 + 4.0 * r))


def gamma(k, u=U64):
    return k * u / (1.0 - k * u)


def residual(A):
    return np.linalg.norm(A @ A - A)


def fp_floor(A, u=U64):
    """Worst-case error of the computed residual (Proposition 3.8)."""
    n = A.shape[0]
    return (1 + u) * gamma(n, u) * np.linalg.norm(A) ** 2


def dist_to_projectors(A):
    """Exact Frobenius distance from symmetric A to the set of orthogonal projectors."""
    w = np.linalg.eigvalsh(A)
    return float(np.linalg.norm(np.minimum(np.abs(w), np.abs(w - 1.0)))), w


# ---------------------------------------------------------------- generators
def random_frame(rng, n, p):
    return np.linalg.qr(rng.standard_normal((n, p)))[0]


def adapted_basis(rng, Q):
    n, p = Q.shape
    return np.linalg.qr(np.hstack([Q, rng.standard_normal((n, n - p))]))[0]


def perturbation(rng, Q, kind, eps):
    """Symmetric perturbation of P=QQ^T of Frobenius norm eps.
    kind: tangent | normal | inward | outward | mixed"""
    n, p = Q.shape
    U = adapted_basis(rng, Q)
    B = np.zeros((n, n))
    if kind in ("tangent", "mixed"):
        D = rng.standard_normal((p, n - p))
        B[:p, p:] = D
        B[p:, :p] = D.T
    if kind in ("normal", "mixed"):
        M = rng.standard_normal((n, n)); M = (M + M.T) / 2
        B[:p, :p] += M[:p, :p]; B[p:, p:] += M[p:, p:]
    if kind in ("inward", "outward"):   # definite diagonal blocks
        Ma = rng.standard_normal((p, p)); Ma = Ma @ Ma.T
        Mc = rng.standard_normal((n - p, n - p)); Mc = Mc @ Mc.T
        s = -1.0 if kind == "inward" else 1.0   # inward: eigenvalues move into (0,1)
        B[:p, :p] = s * Ma; B[p:, p:] = -s * Mc
    B *= eps / np.linalg.norm(B)
    return U @ B @ U.T, U, B


def ms(x):
    x = np.asarray(x, float)
    return {"mean": float(x.mean()), "sd": float(x.std(ddof=1)) if x.size > 1 else 0.0,
            "min": float(x.min()), "max": float(x.max())}


# ---------------------------------------------------------------- experiments
def exp_cert(seeds):
    """Certificate validity and tightness for five perturbation types (n=100, p=10)."""
    n, p = 100, 10
    eps_grid = [1e-1, 3e-2, 1e-2, 1e-3, 1e-4, 1e-6]
    kinds = ["tangent", "outward", "normal", "mixed", "inward"]
    res = []
    for kind in kinds:
        for eps in eps_grid:
            rows = []
            for s in range(seeds):
                rng = np.random.default_rng(1000 * s + 7)
                Q = random_frame(rng, n, p)
                E, _, _ = perturbation(rng, Q, kind, eps)
                A = Q @ Q.T + E
                rho = residual(A)
                if rho >= 0.25:
                    continue
                d, _ = dist_to_projectors(A)
                fl = fp_floor(A)
                rows.append(dict(rho=rho, d=d, lo=psi(rho), hi=phi(rho),
                                 ok=bool(psi(max(rho - fl, 0)) - 1e-15 <= d <= phi(rho + fl) + 1e-15),
                                 d_over_hi=d / phi(rho), d_over_lo=d / psi(rho),
                                 d_over_2rho=d / (2 * rho), d_over_tang=d / (rho + 4 * rho ** 2)))
            if not rows:
                continue
            res.append(dict(kind=kind, eps=eps, n_trials=len(rows),
                            all_ok=all(r["ok"] for r in rows),
                            rho=ms([r["rho"] for r in rows]), d=ms([r["d"] for r in rows]),
                            d_over_phi=ms([r["d_over_hi"] for r in rows]),
                            d_over_psi=ms([r["d_over_lo"] for r in rows]),
                            d_over_2rho=ms([r["d_over_2rho"] for r in rows]),
                            d_over_tang=ms([r["d_over_tang"] for r in rows])))
            print(kind, eps, res[-1]["all_ok"], "%.4f" % res[-1]["d_over_phi"]["mean"], flush=True)
    json.dump(res, open(os.path.join(OUT, "cert.json"), "w"), indent=1)


def exp_sharp(seeds):
    """Controlled equality cases: rank-one Delta (upper constant 1/sqrt2) and equal singular values
    (lower constant 1/sqrt(2r)); single inward / outward eigenvalue (phi / psi attained)."""
    n, p = 100, 10
    out = {}
    for case in ["rank_one", "equal_sv"]:
        for eps in [1e-2, 1e-4]:
            vals = []
            for s in range(seeds):
                rng = np.random.default_rng(31 * s + 3)
                Q = random_frame(rng, n, p)
                U = adapted_basis(rng, Q)
                a = random_frame(rng, p, 1)[:, 0]; b = random_frame(rng, n - p, 1)[:, 0]
                if case == "rank_one":
                    D = np.outer(a, b)
                    r = 1
                else:
                    W = random_frame(rng, p, p); V = random_frame(rng, n - p, p)
                    D = W @ V.T
                    r = p
                B = np.zeros((n, n)); B[:p, p:] = D; B[p:, :p] = D.T
                B *= eps / np.linalg.norm(B)
                A = Q @ Q.T + U @ B @ U.T
                vals.append(residual(A) / eps ** 2)
            pred = 1 / np.sqrt(2) if case == "rank_one" else 1 / np.sqrt(2 * r)
            out[f"{case}_eps{eps:g}"] = dict(ratio=ms(vals), predicted=pred,
                                            max_rel_dev=float(np.max(np.abs(np.array(vals) / pred - 1))))
    # scalar equality in the certificate
    lam_in, lam_out = [], []
    for s in range(seeds):
        rng = np.random.default_rng(77 + s)
        delta = 10 ** rng.uniform(-6, np.log10(0.45))
        Q = random_frame(rng, n, p); U = adapted_basis(rng, Q)
        for sign, store in [(+1, lam_in), (-1, lam_out)]:
            lam = np.r_[np.ones(p), np.zeros(n - p)]
            lam[p] = sign * delta          # first zero moved inward (+) or outward (-)
            A = U @ np.diag(lam) @ U.T
            rho = residual(A); d, _ = dist_to_projectors(A)
            store.append(d / (phi(rho) if sign > 0 else psi(rho)))
    out["single_inward_d_over_phi"] = ms(lam_in)
    out["single_outward_d_over_psi"] = ms(lam_out)
    json.dump(out, open(os.path.join(OUT, "sharp.json"), "w"), indent=1)
    print(json.dumps(out, indent=1))


def exp_sens(seeds):
    """Sensitivity of residual constant and certificate to n and p (tangent Gaussian Delta)."""
    eps = 1e-4
    res = []
    for n in [50, 100, 500, 1024]:
        for p in [2, 5, 10, 20]:
            ratio, tight = [], []
            for s in range(seeds):
                rng = np.random.default_rng(10007 * s + n + p)
                Q = random_frame(rng, n, p)
                E, _, _ = perturbation(rng, Q, "tangent", eps)
                A = Q @ Q.T + E
                rho = residual(A); d, _ = dist_to_projectors(A)
                ratio.append(rho / eps ** 2); tight.append(d / rho)
            res.append(dict(n=n, p=p, ratio=ms(ratio), pred=float(np.sqrt(n / (2 * p * (n - p)))),
                            d_over_rho=ms(tight)))
            print(n, p, "%.4f %.4f" % (res[-1]["ratio"]["mean"], res[-1]["pred"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "sens.json"), "w"), indent=1)


def exp_floor(seeds):
    """Rounding floor of the computed residual for exact projectors vs. the worst-case bound."""
    res = []
    for n in [50, 100, 500, 1024]:
        p = 10
        obs = []
        for s in range(seeds):
            rng = np.random.default_rng(5 * s + n)
            Q = random_frame(rng, n, p)
            A = Q @ Q.T
            A = (A + A.T) / 2
            obs.append(residual(A))
        A = Q @ Q.T
        res.append(dict(n=n, p=p, observed=ms(obs), worst_case=float(fp_floor(A)),
                        sqrt_n_u_p=float(np.sqrt(n) * U64 * p)))
        print(res[-1], flush=True)
    json.dump(res, open(os.path.join(OUT, "floor.json"), "w"), indent=1)


def _digits_cov(rng, m=1500):
    from sklearn.datasets import load_digits
    X = load_digits().data.astype(float)
    idx = rng.choice(X.shape[0], m, replace=False)
    X = X[idx]; X = X - X.mean(0)
    C = X.T @ X / (m - 1)
    return C / np.linalg.norm(C, 2)


def _synthetic_cov(rng, n, p, gap):
    lam = np.r_[np.linspace(2, 1 + gap, p), rng.uniform(0, 1, n - p)]
    Uo = random_frame(rng, n, n)
    return (Uo * lam) @ Uo.T


def _cgs(X):
    """Classical Gram-Schmidt (one pass) -- loses orthogonality like u*kappa(X)^2."""
    n, p = X.shape
    Q = np.zeros_like(X)
    for j in range(p):
        v = X[:, j] - Q[:, :j] @ (Q[:, :j].T @ X[:, j])
        Q[:, j] = v / np.linalg.norm(v)
    return Q


def exp_computed(seeds):
    """Projectors computed by practical routines; the certificate is checked against the exact
    distance, and the frame drift D(X)=||X^T X-I||_F is compared with the distance (Lemma 3.9)."""
    problems = [("digits", 64, 5), ("digits", 64, 10), ("synthetic", 500, 10)]
    methods = ["eigh64", "eigh32", "store16", "cgs_power"]
    res = []
    for prob, n, p in problems:
        for meth in methods:
            rows = []
            for s in range(seeds):
                rng = np.random.default_rng(900 + s)
                C = _digits_cov(rng) if prob == "digits" else _synthetic_cov(rng, n, p, 0.1)
                w, V = np.linalg.eigh(C); Pstar = V[:, -p:] @ V[:, -p:].T
                if meth == "eigh64":
                    X = V[:, -p:]
                elif meth == "eigh32":
                    w32, V32 = np.linalg.eigh(C.astype(np.float32)); X = V32[:, -p:].astype(np.float64)
                elif meth == "store16":
                    X = V[:, -p:].astype(np.float16).astype(np.float64)
                else:   # subspace iteration without re-orthonormalisation, then one CGS pass
                    Y = rng.standard_normal((n, p))
                    for _ in range(6):
                        Y = C @ Y; Y /= np.linalg.norm(Y, axis=0)
                    X = _cgs(Y)
                    Qy = np.linalg.qr(Y)[0]; Pstar = Qy @ Qy.T   # intended subspace = span(Y)
                A = X @ X.T
                rho = residual(A); fl = fp_floor(A)
                d, ev = dist_to_projectors(A)
                drift = np.linalg.norm(X.T @ X - np.eye(p))
                Ps = X @ np.linalg.solve(X.T @ X, X.T)
                rows.append(dict(rho=rho, d=d, lo=psi(max(rho - fl, 0.0)), hi=phi(rho + fl),
                                 drift=drift, sub_err=np.linalg.norm(Ps - Pstar) / np.sqrt(2),
                                 kappa=float(np.linalg.cond(X)),
                                 inward=int(np.sum((ev > 1e-12) & (ev < 1 - 1e-12)))))
            ok = all(r["lo"] - 1e-15 <= r["d"] <= r["hi"] + 1e-15 for r in rows if r["rho"] < 0.25)
            res.append(dict(problem=prob, n=n, p=p, method=meth, all_ok=ok,
                            rho=ms([r["rho"] for r in rows]), d=ms([r["d"] for r in rows]),
                            hi=ms([r["hi"] for r in rows]), drift=ms([r["drift"] for r in rows]),
                            d_minus_drift=ms([abs(r["d"] - r["drift"]) for r in rows]),
                            sub_err=ms([r["sub_err"] for r in rows])))
            print(prob, n, p, meth, ok, "%.2e %.2e %.2e" % (res[-1]["d"]["mean"], res[-1]["hi"]["mean"],
                                                           res[-1]["sub_err"]["mean"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "computed.json"), "w"), indent=1)


# ---------------------------------------------------------------- lazy re-projection
def _mu_shift(A, Q, rng):
    p = Q.shape[1]
    ritz_min = np.linalg.eigvalsh(Q.T @ A @ Q)[0]
    v = rng.standard_normal(A.shape[0])
    for _ in range(30):
        v = v - Q @ (Q.T @ v); v = A @ v; v = v - Q @ (Q.T @ v); v /= np.linalg.norm(v)
    lam = v @ (A @ v)
    return 0.5 * (ritz_min + lam)


def lazy_rgd(A, p, P0, t, tau, policy, shift, tol, maxit, Pstar, rng, monitor=True):
    """Riemannian gradient ascent for max tr(AP) on Gr(n,p) in the projector representation.
    policy: 'always' | 'apriori' (accumulated eps^2/sqrt2) | 'certified' (phi of computed residual)"""
    n = A.shape[0]
    P = P0.copy()
    w, V = np.linalg.eigh(P); Q = V[:, -p:]
    mu = _mu_shift(A, Q, rng) if shift else 0.0
    beta = 0.0; neig = 0; ngemm_cert = 0; maxrho = 0.0; maxd = 0.0
    k = 0
    for k in range(maxit):
        AP = A @ P
        PA = AP.T
        G = PA - PA @ P - mu * (P - P @ P)
        E = t * (G + G.T)
        eps = np.linalg.norm(E)
        if eps / t < tol:
            break
        P = P + E
        P = 0.5 * (P + P.T)
        if policy == "always":
            do = True
        elif policy == "apriori":
            beta += eps ** 2 / np.sqrt(2)
            do = beta > tau
        else:
            R = P @ P - P; ngemm_cert += 1
            rho = np.linalg.norm(R) + fp_floor(P)
            do = (rho >= 0.25) or (phi(rho) > tau)
        if do:
            w, V = np.linalg.eigh(P); Q = V[:, -p:]; P = Q @ Q.T; neig += 1; beta = 0.0
            if shift:
                mu = _mu_shift(A, Q, rng)
        if monitor:   # defect of the iterate actually carried to the next step
            maxrho = max(maxrho, residual(P))
            maxd = max(maxd, dist_to_projectors(P)[0])
    w, V = np.linalg.eigh(P); Q = V[:, -p:]; Pf = Q @ Q.T
    return dict(iters=k, eig=neig, gemm_cert=ngemm_cert, max_rho=maxrho, max_dist=maxd,
                final_err=float(np.linalg.norm(Pf - Pstar) / np.sqrt(2)))


POLICIES = [("always", True), ("apriori", True), ("certified", True), ("apriori", False), ("certified", False)]


def _problem(name, rng):
    if name.startswith("syn"):
        _, n, g = name.split("_"); n = int(n); g = float(g); p = 10
        A = _synthetic_cov(rng, n, p, g)
    else:
        p = int(name.split("_")[1]); A = _digits_cov(rng, 1797); n = A.shape[0]
    w, V = np.linalg.eigh(A)
    return A, p, V[:, -p:] @ V[:, -p:].T, 1.0 / (w[-1] - w[0])


def exp_lazy_counts(seeds, problems=("syn_300_0.5", "syn_300_0.1", "digits_5", "digits_10"), tau=1e-10):
    res = []
    for prob in problems:
        for pol, sh in POLICIES:
            rows = []
            for s in range(seeds):
                rng = np.random.default_rng(4242 + s)
                A, p, Pstar, t = _problem(prob, rng)
                Q0 = random_frame(rng, A.shape[0], p)
                rows.append(lazy_rgd(A, p, Q0 @ Q0.T, t, tau, pol, sh, 1e-10, 5000, Pstar,
                                     np.random.default_rng(s)))
            res.append(dict(problem=prob, policy=pol, shift=sh,
                            **{k: ms([r[k] for r in rows]) for k in rows[0]}))
            r = res[-1]
            print(prob, pol, sh, "it %.1f eig %.1f maxd %.2e err %.2e" % (
                r["iters"]["mean"], r["eig"]["mean"], r["max_dist"]["max"], r["final_err"]["max"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "lazy_counts.json"), "w"), indent=1)


def exp_lazy_time(seeds, problems=("syn_1000_0.5", "syn_1000_0.1"), tau=1e-10):
    res = []
    for prob in problems:
        for pol, sh in POLICIES[:3]:
            rows = []
            for s in range(seeds):
                rng = np.random.default_rng(4242 + s)
                A, p, Pstar, t = _problem(prob, rng)
                Q0 = random_frame(rng, A.shape[0], p)
                t0 = time.perf_counter()
                r = lazy_rgd(A, p, Q0 @ Q0.T, t, tau, pol, sh, 1e-10, 5000, Pstar,
                             np.random.default_rng(s), monitor=False)
                r["time"] = time.perf_counter() - t0
                rows.append(r)
            res.append(dict(problem=prob, policy=pol, shift=sh,
                            **{k: ms([r[k] for r in rows]) for k in rows[0]}))
            r = res[-1]
            print(prob, pol, "it %.1f eig %.1f time %.1f err %.2e" % (
                r["iters"]["mean"], r["eig"]["mean"], r["time"]["mean"], r["final_err"]["max"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "lazy_time.json"), "w"), indent=1)


def exp_trace(seeds=1, prob="syn_300_0.1", tau=1e-10):
    """Per-iteration defect traces for Figure 3 (one start)."""
    out = {}
    for pol, sh in [("apriori", False), ("certified", False), ("certified", True)]:
        rng = np.random.default_rng(4242)
        A, p, Pstar, t = _problem(prob, rng)
        Q0 = random_frame(rng, A.shape[0], p)
        P = Q0 @ Q0.T
        tr = []
        # re-implementation with per-step logging (same logic as lazy_rgd)
        w, V = np.linalg.eigh(P); Q = V[:, -p:]
        r2 = np.random.default_rng(0)
        mu = _mu_shift(A, Q, r2) if sh else 0.0
        beta = 0.0; neig = 0
        for k in range(5000):
            PA = (A @ P).T
            G = PA - PA @ P - mu * (P - P @ P)
            E = t * (G + G.T); eps = np.linalg.norm(E)
            if eps / t < 1e-10:
                break
            P = 0.5 * ((P + E) + (P + E).T)
            if pol == "apriori":
                beta += eps ** 2 / np.sqrt(2); do = beta > tau
            else:
                rho = np.linalg.norm(P @ P - P) + fp_floor(P); do = rho >= 0.25 or phi(rho) > tau
            if do:
                w, V = np.linalg.eigh(P); Q = V[:, -p:]; P = Q @ Q.T; neig += 1; beta = 0.0
                if sh:
                    mu = _mu_shift(A, Q, r2)
            tr.append([k, dist_to_projectors(P)[0], eps ** 2 / np.sqrt(2), neig])
        out[f"{pol}_{'shift' if sh else 'noshift'}"] = tr
        print(pol, sh, len(tr), neig, max(x[1] for x in tr), flush=True)
    json.dump(out, open(os.path.join(OUT, "trace.json"), "w"))


# ---------------------------------------------------------------- template benchmarks (B1, B2) and sensitivity
U32 = 2.0 ** -24


def lazy_rgd_p(A, p, P0, eta, tau, policy, shift, tol, maxit, Pstar, rng, dtype=np.float64):
    """lazy_rgd in a given working precision (float64 or float32); returns iterations, eigendecompositions,
    final subspace distance, wall time and whether the tolerance was reached."""
    u = U64 if dtype == np.float64 else U32
    A = A.astype(dtype); P = P0.astype(dtype)
    n = A.shape[0]
    w, V = np.linalg.eigh(P); Q = V[:, -p:]
    mu = dtype(_mu_shift(A.astype(np.float64), Q.astype(np.float64), rng)) if shift else dtype(0)
    beta = 0.0; neig = 0; conv = False
    t0 = time.perf_counter()
    for k in range(maxit):
        PA = (A @ P).T
        G = PA - PA @ P - mu * (P - P @ P)
        E = dtype(eta) * (G + G.T)
        eps = float(np.linalg.norm(E))
        if not np.isfinite(eps):
            break
        if eps / eta < tol:
            conv = True
            break
        P = P + E; P = (P + P.T) / dtype(2)
        if policy == "always":
            do = True
        elif policy == "apriori":
            beta += eps ** 2 / np.sqrt(2); do = beta > tau
        else:
            R = P @ P - P
            nrm = float(np.linalg.norm(P))
            rho = (float(np.linalg.norm(R)) + (1 + u) * gamma(n, u) * nrm ** 2) / (1 - u)
            do = (rho >= 0.25) or (phi(rho) > tau)
        if do:
            w, V = np.linalg.eigh(P); Q = V[:, -p:]; P = Q @ Q.T; neig += 1; beta = 0.0
            if shift:
                mu = dtype(_mu_shift(A.astype(np.float64), Q.astype(np.float64), rng))
    el = time.perf_counter() - t0
    P64 = P.astype(np.float64); P64 = (P64 + P64.T) / 2
    w, V = np.linalg.eigh(P64); Q = V[:, -p:]
    return dict(iters=k, eig=neig, final_dist=float(np.linalg.norm(Q @ Q.T - Pstar) / np.sqrt(2)),
                time=el, conv=conv)


def _bench_problem(name, rng):
    if name == "B1":
        A = _synthetic_cov(rng, 100, 10, 0.1); p = 10
    else:   # B2 substitute: UCI digits covariance
        A = _digits_cov(rng, 1797); p = 10
    w, V = np.linalg.eigh(A)
    return A, p, V[:, -p:] @ V[:, -p:].T, 1.0 / (w[-1] - w[0])


def exp_bench(seeds):
    """Template Table 2: benchmark suite, float64 and float32."""
    res = []
    settings = {np.float64: (1e-10, 1e-10), np.float32: (2e-4, 1e-5)}
    pols = [("always", True, 1.0), ("apriori", True, 1.0), ("certified", True, 1.0), ("certified", False, 1.0)]
    for prob in ["B1", "B2"]:
        for dt in [np.float64, np.float32]:
            tau, tol = settings[dt]
            for pol, sh, _ in pols + [("certified0", True, 1.0)]:
                rows = []
                for s in range(seeds):
                    rng = np.random.default_rng(7000 + s)
                    A, p, Pstar, eta = _bench_problem(prob, rng)
                    Q0 = random_frame(rng, A.shape[0], p)
                    pp = "certified" if pol == "certified0" else pol
                    tt = 0.0 if pol == "certified0" else tau
                    rows.append(lazy_rgd_p(A, p, Q0 @ Q0.T, eta, tt, pp, sh, tol, 5000, Pstar,
                                           np.random.default_rng(s), dt))
                if pol == "certified0" and dt != np.float64:
                    continue
                r = dict(problem=prob, dtype=dt.__name__, policy=pol, shift=sh, tau=tau if pol != "certified0" else 0.0,
                         final_dist=ms([x["final_dist"] for x in rows]), iters=ms([x["iters"] for x in rows]),
                         eig=ms([x["eig"] for x in rows]), time=ms([x["time"] for x in rows]),
                         divergences=int(sum(not x["conv"] for x in rows)))
                res.append(r)
                print(prob, dt.__name__, pol, sh, "it %.1f eig %.1f dist %.2e div %d" % (
                    r["iters"]["mean"], r["eig"]["mean"], r["final_dist"]["mean"], r["divergences"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "bench.json"), "w"), indent=1)


def exp_sens_alg(seeds):
    """Template Table 3: sensitivity of Algorithm 1 (certified + shift, B1, float64) to step size and tau."""
    res = []
    for c in [0.25, 0.5, 1.0, 1.5, 2.0]:
        rows = []
        for s in range(seeds):
            rng = np.random.default_rng(7000 + s)
            A, p, Pstar, eta = _bench_problem("B1", rng)
            Q0 = random_frame(rng, A.shape[0], p)
            rows.append(lazy_rgd_p(A, p, Q0 @ Q0.T, c * eta, 1e-10, "certified", True, 1e-10, 5000, Pstar,
                                   np.random.default_rng(s)))
            rows[-1]["always"] = lazy_rgd_p(A, p, Q0 @ Q0.T, c * eta, 0.0, "always", True, 1e-10, 5000, Pstar,
                                            np.random.default_rng(s))["eig"]
        res.append(dict(factor="step", value=c, iters=ms([x["iters"] for x in rows]), eig=ms([x["eig"] for x in rows]),
                        saved=ms([1 - x["eig"] / max(x["always"], 1) for x in rows]),
                        final_dist=ms([x["final_dist"] for x in rows]), divergences=int(sum(not x["conv"] for x in rows))))
        print("step", c, res[-1]["iters"]["mean"], res[-1]["saved"]["mean"], res[-1]["divergences"], flush=True)
    for tau in [1e-14, 1e-12, 1e-10, 1e-8, 1e-6]:
        rows = []
        for s in range(seeds):
            rng = np.random.default_rng(7000 + s)
            A, p, Pstar, eta = _bench_problem("B1", rng)
            Q0 = random_frame(rng, A.shape[0], p)
            r = lazy_rgd_p(A, p, Q0 @ Q0.T, eta, tau, "certified", True, 1e-10, 5000, Pstar, np.random.default_rng(s))
            r["always"] = r["iters"]
            rows.append(r)
        res.append(dict(factor="tau", value=tau, iters=ms([x["iters"] for x in rows]), eig=ms([x["eig"] for x in rows]),
                        saved=ms([1 - x["eig"] / max(x["iters"], 1) for x in rows]),
                        final_dist=ms([x["final_dist"] for x in rows]), divergences=int(sum(not x["conv"] for x in rows))))
        print("tau", tau, res[-1]["iters"]["mean"], res[-1]["saved"]["mean"], res[-1]["divergences"], flush=True)
    json.dump(res, open(os.path.join(OUT, "sens_alg.json"), "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("exp")
    ap.add_argument("--seeds", type=int, default=30)
    a = ap.parse_args()
    f = {"cert": exp_cert, "sharp": exp_sharp, "sens": exp_sens, "floor": exp_floor,
         "computed": exp_computed, "lazy_counts": exp_lazy_counts, "lazy_time": exp_lazy_time,
         "trace": lambda s: exp_trace(), "bench": exp_bench, "sens_alg": exp_sens_alg}
    if a.exp == "all":
        for k in ["sharp", "cert", "sens", "floor", "computed", "lazy_counts", "lazy_time"]:
            f[k](a.seeds if k != "lazy_time" else 5)
    else:
        f[a.exp](a.seeds)
