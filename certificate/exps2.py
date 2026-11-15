"""Additional experiments for the revised manuscript (referee comments):
baselines (McWeeny purification, periodic re-projection, frame-based QR retraction),
wall-clock timing with 10 seeds, the shift lemma, and the accuracy theorem.
Usage: python exps2.py <baselines|time2|shift|accuracy> [--seeds N]"""
import argparse, json, os, time
import numpy as np
from exps import (U64, OUT, cert_skip, dist_to_projectors, random_frame, ms, _mu_shift, _problem, gamma)


def run(A, p, P0, eta, tau, policy, shift, tol, maxit, Pstar, rng, monitor=True, hist=False):
    """Projector-form RGD with a re-projection policy; matrix products are reused where possible:
    at an exact projector the shift term vanishes (2 products per step), otherwise P^2 is needed
    (3 products), and the P^2 formed by the certified test is reused by the next step.
    policy: always | certified | apriori | purify (one McWeeny step per iteration) |
            cert_purify (certified test; on failure McWeeny steps, eigendecomposition as fallback) | periodic_k"""
    n = A.shape[0]
    w, V = np.linalg.eigh(P0); Q = V[:, -p:]; P = Q @ Q.T
    mu = _mu_shift(A, Q, rng) if shift else 0.0
    beta = 0.0; neig = gemm = npur = guard = 0; maxd = 0.0; H = []
    exact = True; PPc = None
    period = int(policy.split("_")[1]) if policy.startswith("periodic") else 0
    t0 = time.perf_counter()
    for k in range(maxit):
        PA = (A @ P).T; PAP = PA @ P; gemm += 2
        if exact:
            G = PA - PAP
        else:
            if PPc is None:
                PPc = P @ P; gemm += 1
            G = PA - PAP - mu * (P - PPc)
        E = eta * (G + G.T); eps = float(np.linalg.norm(E))
        if eps / eta < tol:
            break
        P = P + E; P = 0.5 * (P + P.T); exact = False; PPc = None
        do = False
        if policy == "always":
            do = True
        elif policy == "apriori":
            beta += eps ** 2 / np.sqrt(2); do = beta > tau
        elif period:
            do = (k + 1) % period == 0
        elif policy == "certified":
            PPc = P @ P; gemm += 1
            keep, g = cert_skip(P, PPc - P, eps, tau); guard += g; do = not keep
        elif policy == "purify":
            P2 = P @ P; P = 3 * P2 - 2 * (P2 @ P); P = 0.5 * (P + P.T); gemm += 2; npur += 1
        elif policy == "cert_purify":
            e0 = eps
            for it in range(4):
                PPc = P @ P; gemm += 1
                keep, g = cert_skip(P, PPc - P, e0, tau); guard += g; e0 = 0.0
                if keep:
                    break
                if it == 3 or g:
                    do = True; break
                P = 3 * PPc - 2 * (PPc @ P); P = 0.5 * (P + P.T); gemm += 1; npur += 1; PPc = None
        if do:
            w, V = np.linalg.eigh(P); Q = V[:, -p:]; P = Q @ Q.T; neig += 1; beta = 0.0; exact = True; PPc = None
            if shift:
                mu = _mu_shift(A, Q, rng)
        if monitor:
            maxd = max(maxd, dist_to_projectors(P)[0])
        if hist:
            w, V = np.linalg.eigh(P); Qk = V[:, -p:]
            H.append(float(np.linalg.norm(Qk @ Qk.T - Pstar)))
    el = time.perf_counter() - t0
    w, V = np.linalg.eigh(P); Q = V[:, -p:]
    out = dict(iters=k, eig=neig, gemm=gemm, purif=npur, guard=guard, max_dist=maxd, time=el,
               work=gemm + 4.5 * neig, final_err=float(np.linalg.norm(Q @ Q.T - Pstar) / np.sqrt(2)))
    if hist:
        out["hist"] = H
    return out


def run_frame(A, p, X0, eta, tol, maxit, Pstar):
    """Frame (Stiefel) representation with QR retraction, the usual implementation."""
    X = X0.copy(); t0 = time.perf_counter()
    for k in range(maxit):
        AX = A @ X
        D = AX - X @ (X.T @ AX)
        if np.sqrt(2) * np.linalg.norm(D) < tol:     # same scale as ||E||/eta for the projector step
            break
        X, _ = np.linalg.qr(X + eta * D)
    el = time.perf_counter() - t0
    return dict(iters=k, time=el, final_err=float(np.linalg.norm(X @ X.T - Pstar) / np.sqrt(2)))


POL = [("always", True), ("certified", True), ("apriori", True), ("periodic_10", True), ("purify", True),
       ("cert_purify", True), ("certified", False), ("apriori", False)]
SEEDS = {"syn_300_0.5": 10, "syn_300_0.1": 10}


def exp_baselines(seeds, problems=("digits_10", "digits_5", "syn_300_0.5"), tau=1e-10):
    res = []
    for prob in problems:
        for pol, sh in POL:
            rows = []
            for s in range(SEEDS.get(prob, seeds)):
                rng = np.random.default_rng(4242 + s)
                A, p, Pstar, t = _problem(prob, rng)
                Q0 = random_frame(rng, A.shape[0], p)
                rows.append(run(A, p, Q0 @ Q0.T, t, tau, pol, sh, 1e-10, 5000, Pstar, np.random.default_rng(s)))
            r = dict(problem=prob, policy=pol, shift=sh, seeds=len(rows), **{k: ms([x[k] for x in rows]) for k in rows[0]},
                     nonconv=int(sum(x["iters"] >= 4999 for x in rows)))
            res.append(r)
            print(prob, pol, "it %.1f eig %.1f gemm %.1f work %.1f maxd %.2e err %.2e nc %d guard %d" % (
                r["iters"]["mean"], r["eig"]["mean"], r["gemm"]["mean"], r["work"]["mean"], r["max_dist"]["max"],
                r["final_err"]["max"], r["nonconv"], r["guard"]["max"]), flush=True)
            json.dump(res, open(os.path.join(OUT, "baselines.json"), "w"), indent=1)


def exp_time2(seeds, tau=1e-10):
    res = []
    plan = [("syn_1000_0.5", ["always", "certified", "apriori", "periodic_5", "purify", "cert_purify", "frame_qr"], seeds),
            ("syn_1000_0.1", ["always", "certified", "periodic_5", "purify", "frame_qr"], max(seeds // 2, 1))]
    for prob, pols, ns in plan:
        for pol in pols:
            rows = []
            for s in range(ns):
                rng = np.random.default_rng(4242 + s)
                A, p, Pstar, t = _problem(prob, rng)
                Q0 = random_frame(rng, A.shape[0], p)
                if pol == "frame_qr":
                    rows.append(run_frame(A, p, Q0, t, 1e-10, 5000, Pstar))
                else:
                    rows.append(run(A, p, Q0 @ Q0.T, t, tau, pol, True, 1e-10, 5000, Pstar,
                                    np.random.default_rng(s), monitor=False))
            times = [x["time"] for x in rows]
            r = dict(problem=prob, policy=pol, **{k: ms([x[k] for x in rows]) for k in rows[0]},
                     time_median=float(np.median(times)), time_q1=float(np.percentile(times, 25)),
                     time_q3=float(np.percentile(times, 75)))
            res.append(r)
            print(prob, pol, "it %.1f time med %.2f err %.2e" % (r["iters"]["mean"], r["time_median"],
                                                                 r["final_err"]["max"]), flush=True)
            json.dump(res, open(os.path.join(OUT, "time2.json"), "w"), indent=1)


def exp_shift(seeds):
    """Lemma 5.3: growth factor of a normal defect without re-projection, observed vs predicted."""
    res = []
    for prob in ["syn_300_0.5", "digits_10"]:
        for sh in [False, True]:
            rows = []
            for s in range(seeds):
                rng = np.random.default_rng(4242 + s)
                A, p, Pstar, t = _problem(prob, rng)
                n = A.shape[0]
                w, V = np.linalg.eigh(A); U = V[:, ::-1]
                Q = U[:, :p] @ U[:, :p].T
                mu = _mu_shift(A, U[:, :p], np.random.default_rng(s)) if sh else 0.0
                th = w[::-1][:p] - mu; om = w[::-1][p:] - mu
                pred = max(np.max(np.abs(1 - t * (th[:, None] + th[None, :]))),
                           np.max(np.abs(1 + t * (om[:, None] + om[None, :]))))
                d = np.zeros(n); d[:p] = rng.standard_normal(p); d[p:] = rng.standard_normal(n - p)
                N = U @ np.diag(1e-9 * d / np.linalg.norm(d)) @ U.T
                P = Q + N; ds = []
                for k in range(40):
                    PA = (A @ P).T
                    G = PA - PA @ P - mu * (P - P @ P)
                    P = P + t * (G + G.T); P = 0.5 * (P + P.T)
                    ds.append(dist_to_projectors(P)[0])
                ds = np.array(ds); good = np.where(ds < 1e-5)[0]
                k1 = int(good[-1]) if len(good) else 10
                k0 = 5
                obs = (ds[k1] / ds[k0]) ** (1 / max(k1 - k0, 1))
                rows.append(dict(pred=float(pred), obs=float(obs), d40=ds[-1]))
            r = dict(problem=prob, shift=sh, pred=ms([x["pred"] for x in rows]), obs=ms([x["obs"] for x in rows]),
                     rel=ms([abs(x["obs"] / x["pred"] - 1) for x in rows]), d40=ms([x["d40"] for x in rows]))
            res.append(r)
            print(prob, sh, "pred %.6f obs %.6f" % (r["pred"]["mean"], r["obs"]["mean"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "shift.json"), "w"), indent=1)


def exp_accuracy(seeds, tau=1e-10):
    """Theorem 5.4: final error of the certified run against the bound q^k e0 + K tau/(1-q)."""
    res = []
    for prob in ["digits_5", "digits_10", "syn_300_0.5"]:
        rows = []
        for s in range(SEEDS.get(prob, seeds)):
            rng = np.random.default_rng(4242 + s)
            A, p, Pstar, t = _problem(prob, rng)
            Q0 = random_frame(rng, A.shape[0], p)
            a = run(A, p, Q0 @ Q0.T, t, 0.0, "always", True, 1e-10, 5000, Pstar, np.random.default_rng(s), hist=True)
            c = run(A, p, Q0 @ Q0.T, t, tau, "certified", True, 1e-10, 5000, Pstar, np.random.default_rng(s), hist=True)
            h = np.array(a["hist"]); m = len(h)
            i0, i1 = int(0.3 * m), int(0.7 * m)
            q = float((h[i1] / h[i0]) ** (1 / (i1 - i0)))
            w = np.linalg.eigvalsh(A); B = w[-1] - w[0]
            c1 = 1 + 2 * t * (3 + tau) * B
            K = np.sqrt(2) * c1 / (1 - c1 * tau)
            L = min(len(a["hist"]), len(c["hist"]))
            diff = float(np.max(np.abs(np.array(a["hist"][:L]) - np.array(c["hist"][:L]))))
            rows.append(dict(q=q, K=K, floor=K * tau / (1 - q), diff=diff,
                             err_always=a["final_err"], err_cert=c["final_err"]))
        r = dict(problem=prob, seeds=len(rows), **{k: ms([x[k] for x in rows]) for k in rows[0]},
                 max_diff=max(x["diff"] for x in rows), max_floor=max(x["floor"] for x in rows),
                 min_floor=min(x["floor"] for x in rows))
        res.append(r)
        print(prob, "q %.4f K %.2f floor %.2e diff %.2e" % (r["q"]["mean"], r["K"]["mean"], r["floor"]["mean"],
                                                            r["max_diff"]), flush=True)
    json.dump(res, open(os.path.join(OUT, "accuracy.json"), "w"), indent=1)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("exp"); ap.add_argument("--seeds", type=int, default=30)
    a = ap.parse_args()
    {"baselines": exp_baselines, "time2": exp_time2, "shift": exp_shift, "accuracy": exp_accuracy}[a.exp](a.seeds)
