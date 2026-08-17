"""
Main experiment script for Grassmann projector perturbation analysis.
Author: Pascal Kubwimana
Purpose: Verify quadratic vs. linear residual scaling under tangent and non-tangent perturbations.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# Create results directory if it doesn't exist
if not os.path.exists('../results'):
    os.makedirs('../results')

# -----------------------------
# 1. Generate random projector
# -----------------------------
def random_projector(n, p):
    """
    Construct a random orthogonal projector of rank p in R^n.
    """
    Q, _ = np.linalg.qr(np.random.randn(n, p))
    return Q @ Q.T

# -----------------------------
# 2. Tangent perturbation
# -----------------------------
def tangent_perturbation(P, p, n, epsilon):
    """
    Generate a tangent perturbation scaled to Frobenius norm epsilon.
    """
    # Use full n x n random matrix
    Delta = np.random.randn(n, n)
    E = P @ Delta @ (np.eye(n) - P) + (np.eye(n) - P) @ Delta.T @ P
    E = E * (epsilon / np.linalg.norm(E, 'fro'))
    return E

# -----------------------------
# 3. Non-tangent perturbation
# -----------------------------
def nontangent_perturbation(n, epsilon):
    """
    Generate a symmetric non-tangent perturbation scaled to Frobenius norm epsilon.
    """
    E = np.random.randn(n, n)
    E = (E + E.T) / 2  # enforce symmetry
    E = E * (epsilon / np.linalg.norm(E, 'fro'))
    return E

# -----------------------------
# 4. Residual and backward error
# -----------------------------
def residual(P_tilde):
    """
    Compute residual norm ||P̃² - P̃||_F.
    """
    return np.linalg.norm(P_tilde @ P_tilde - P_tilde, 'fro')

def backward_error(P_tilde, p):
    """
    Compute backward error relative to nearest rank-p projector.
    """
    eigvals, eigvecs = np.linalg.eigh(P_tilde)
    idx = np.argsort(eigvals)[::-1][:p]
    U = eigvecs[:, idx]
    P_star = U @ U.T
    return np.linalg.norm(P_tilde - P_star, 'fro')

# -----------------------------
# 5. Experiment loop
# -----------------------------
def experiment(n=100, p=10, epsilons=np.logspace(-6, -1, 6), trials=30):
    """
    Run perturbation experiments across multiple epsilon values.
    Returns a list of dictionaries with averaged results.
    """
    results = []
    for eps in epsilons:
        res_tangent, err_tangent = [], []
        res_nontangent, err_nontangent = [], []
        
        for _ in range(trials):
            P = random_projector(n, p)

            # Tangent perturbation
            E_t = tangent_perturbation(P, p, n, eps)
            P_tilde_t = P + E_t
            res_tangent.append(residual(P_tilde_t))
            err_tangent.append(backward_error(P_tilde_t, p))

            # Non-tangent perturbation
            E_n = nontangent_perturbation(n, eps)
            P_tilde_n = P + E_n
            res_nontangent.append(residual(P_tilde_n))
            err_nontangent.append(backward_error(P_tilde_n, p))

        results.append({
            "epsilon": eps,
            "tangent_residual": np.mean(res_tangent),
            "tangent_error": np.mean(err_tangent),
            "nontangent_residual": np.mean(res_nontangent),
            "nontangent_error": np.mean(err_nontangent),
            "tangent_ratio": np.mean(res_tangent) / (eps**2),
            "nontangent_ratio": np.mean(res_nontangent) / eps
        })
    return results

# -----------------------------
# 6. Rank dependence experiment
# -----------------------------
def rank_dependence_experiment(n=100, p_values=[1, 5, 10, 20, 40], epsilon=1e-6, trials=30):
    """
    Run experiments with varying rank p.
    """
    results = []
    
    for p in p_values:
        res_tangent, err_tangent = [], []
        
        for _ in range(trials):
            P = random_projector(n, p)
            E_t = tangent_perturbation(P, p, n, epsilon)
            P_tilde = P + E_t
            res_tangent.append(residual(P_tilde))
            err_tangent.append(backward_error(P_tilde, p))
        
        results.append({
            "p": p,
            "residual": np.mean(res_tangent),
            "backward_error": np.mean(err_tangent),
            "ratio": np.mean(err_tangent) / (epsilon**2)
        })
    
    return results

# -----------------------------
# 7. CIFAR-10 validation experiment
# -----------------------------
def cifar10_projector():
    """
    Load CIFAR-10 data and create PCA projector.
    Note: You need to have CIFAR-10 data downloaded.
    If not available, this returns a random projector as fallback.
    """
    try:
        from sklearn.decomposition import PCA
        from sklearn.datasets import fetch_openml
        
        print("Loading CIFAR-10 data...")
        # Load CIFAR-10 from OpenML (8,000 samples)
        X, y = fetch_openml('CIFAR_10', version=1, return_X_y=True, as_frame=False, parser='pandas')
        X = X[:8000]  # Use 8,000 samples
        X = X.reshape(-1, 32*32)  # Flatten images
        
        # Center the data
        X = X - X.mean(axis=0)
        
        # PCA
        pca = PCA(n_components=10)
        pca.fit(X)
        
        # Create projector
        Q = pca.components_.T
        P = Q @ Q.T
        
        print(f"CIFAR-10 projector created. Explained variance: {pca.explained_variance_ratio_.sum():.3f}")
        return P
    except:
        print("CIFAR-10 not available. Using random projector as fallback.")
        return random_projector(1024, 10)

def cifar10_experiment(epsilons=np.logspace(-6, -1, 4), trials=50):
    """
    Run experiment on CIFAR-10 derived projector.
    """
    P = cifar10_projector()
    n = P.shape[0]
    p = 10
    
    results = []
    for eps in epsilons:
        res_tangent, err_tangent = [], []
        res_nontangent, err_nontangent = [], []
        
        for _ in range(trials):
            # Tangent perturbation
            E_t = tangent_perturbation(P, p, n, eps)
            P_tilde_t = P + E_t
            res_tangent.append(residual(P_tilde_t))
            err_tangent.append(backward_error(P_tilde_t, p))

            # Non-tangent perturbation
            E_n = nontangent_perturbation(n, eps)
            P_tilde_n = P + E_n
            res_nontangent.append(residual(P_tilde_n))
            err_nontangent.append(backward_error(P_tilde_n, p))

        results.append({
            "epsilon": eps,
            "tangent_residual": np.mean(res_tangent),
            "tangent_error": np.mean(err_tangent),
            "nontangent_residual": np.mean(res_nontangent),
            "nontangent_error": np.mean(err_nontangent),
            "tangent_ratio": np.mean(res_tangent) / (eps**2)
        })
    return results

# -----------------------------
# 8. Plotting functions
# -----------------------------
def plot_residual_scaling(df, save_path="../results/residual_scaling_plot.png"):
    """
    Create log-log plot of residual scaling.
    """
    plt.figure(figsize=(8,6))
    plt.loglog(df["epsilon"], df["tangent_residual"], 'o-', 
               label="Tangent residual (O(ε²))", linewidth=2, markersize=8)
    plt.loglog(df["epsilon"], df["nontangent_residual"], 's-', 
               label="Non-tangent residual (O(ε))", linewidth=2, markersize=8)
    plt.loglog(df["epsilon"], df["epsilon"]**2, 'k--', 
               label="ε² reference", linewidth=1.5, alpha=0.7)
    plt.loglog(df["epsilon"], df["epsilon"], 'r--', 
               label="ε reference", linewidth=1.5, alpha=0.7)
    
    plt.xlabel("Perturbation norm ε", fontsize=12)
    plt.ylabel("Residual ∥P̃² - P̃∥_F", fontsize=12)
    plt.title("Quadratic vs Linear Scaling of Residuals", fontsize=14)
    plt.legend(loc='best', fontsize=11)
    plt.grid(True, which="both", ls="--", alpha=0.6)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Plot saved to {save_path}")
    plt.show()

def plot_rank_dependence(results, save_path="../results/rank_dependence.png"):
    """
    Plot backward error vs rank p.
    """
    df = pd.DataFrame(results)
    
    plt.figure(figsize=(8,6))
    plt.plot(df["p"], df["ratio"], 'o-', linewidth=2, markersize=8)
    plt.xlabel("Rank p", fontsize=12)
    plt.ylabel("Backward Error / ε²", fontsize=12)
    plt.title("Rank Dependence of Backward Error (ε = 1e-6)", fontsize=14)
    plt.grid(True, alpha=0.6)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300)
    print(f"Rank dependence plot saved to {save_path}")
    plt.show()

# -----------------------------
# 9. Main execution
# -----------------------------
if __name__ == "__main__":
    print("="*60)
    print("GRASSMANN PROJECTOR PERTURBATION ANALYSIS")
    print("="*60)
    
    # Run main experiment
    print("\n[1/3] Running main perturbation experiment...")
    print("n=100, p=10, ε ∈ [1e-6, 1e-1], trials=30")
    results = experiment()
    df = pd.DataFrame(results)
    
    # Save results
    csv_path = "../results/grassmann_experiment_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nResults saved to {csv_path}")
    
    # Display table
    print("\n" + "="*60)
    print("TABLE: Residuals and Backward Errors vs ε")
    print("="*60)
    print(df.to_string(index=False, float_format="%.3e"))
    
    # Plot residual scaling
    print("\n[2/3] Generating residual scaling plot...")
    plot_residual_scaling(df)
    
    # Run rank dependence
    print("\n[3/3] Running rank dependence experiment...")
    print("p ∈ [1, 5, 10, 20, 40], ε = 1e-6")
    rank_results = rank_dependence_experiment()
    rank_df = pd.DataFrame(rank_results)
    
    print("\n" + "="*60)
    print("RANK DEPENDENCE at ε = 1e-6")
    print("="*60)
    print(rank_df.to_string(index=False, float_format="%.3e"))
    
    # Plot rank dependence
    plot_rank_dependence(rank_results)
    
    # Optional: Run CIFAR-10 experiment
    print("\n" + "="*60)
    print("CIFAR-10 VALIDATION")
    print("="*60)
    print("Running CIFAR-10 experiment (this may take a moment)...")
    try:
        cifar_results = cifar10_experiment()
        cifar_df = pd.DataFrame(cifar_results)
        print("\nCIFAR-10 Results:")
        print(cifar_df.to_string(index=False, float_format="%.3e"))
        
        # Save CIFAR-10 results
        cifar_df.to_csv("../results/cifar10_results.csv", index=False)
        print("\nCIFAR-10 results saved to ../results/cifar10_results.csv")
    except Exception as e:
        print(f"CIFAR-10 experiment skipped: {e}")
    
    print("\n" + "="*60)
    print("EXPERIMENT COMPLETE!")
    print("All results and figures saved to 'results/' folder")
    print("="*60)
