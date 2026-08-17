"""
Complete experiment script for Grassmann projector perturbation analysis.
Produces both table and plot outputs.
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
# 6. Run experiment and produce outputs
# -----------------------------
if __name__ == "__main__":
    print("="*70)
    print("GRASSMANN PROJECTOR PERTURBATION ANALYSIS")
    print("n=100, p=10, trials=30")
    print("="*70)
    
    # Run experiment
    results = experiment()
    df = pd.DataFrame(results)
    
    # -----------------------------
    # 6a. Produce Table
    # -----------------------------
    print("\n" + "="*70)
    print("TABLE: Residuals and Backward Errors vs ε")
    print("="*70)
    print(df.to_string(index=False, float_format="%.3e"))
    
    # Save table to CSV
    csv_path = "../results/grassmann_experiment_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"\nTable saved to {csv_path}")
    
    # -----------------------------
    # 6b. Produce Plot
    # -----------------------------
    print("\n" + "="*70)
    print("GENERATING PLOT: Residual Scaling")
    print("="*70)
    
    plt.figure(figsize=(10, 7))
    
    # Plot data
    plt.loglog(df["epsilon"], df["tangent_residual"], 'o-', 
               color='blue', linewidth=2.5, markersize=10,
               label="Tangent residual (O(ε²))")
    plt.loglog(df["epsilon"], df["nontangent_residual"], 's-', 
               color='orange', linewidth=2.5, markersize=10,
               label="Non-tangent residual (O(ε))")
    
    # Reference lines
    plt.loglog(df["epsilon"], df["epsilon"]**2, 'k--', 
               linewidth=1.5, alpha=0.7, label="ε² reference")
    plt.loglog(df["epsilon"], df["epsilon"], 'r--', 
               linewidth=1.5, alpha=0.7, label="ε reference")
    
    # Labels and title
    plt.xlabel("Perturbation norm ε", fontsize=14)
    plt.ylabel("Residual ∥P̃² - P̃∥_F", fontsize=14)
    plt.title("Quadratic vs Linear Scaling of Residuals", fontsize=16)
    
    # Legend
    plt.legend(loc='lower right', fontsize=12)
    
    # Grid
    plt.grid(True, which="both", ls="--", alpha=0.3)
    
    # Add annotation with parameters
    plt.text(0.02, 0.98, f"n=100, p=10, trials=30", 
             transform=plt.gca().transAxes, fontsize=11,
             verticalalignment='top', 
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Save plot
    plot_path = "../results/residual_scaling_plot.png"
    plt.savefig(plot_path, dpi=300, bbox_inches='tight')
    print(f"Plot saved to {plot_path}")
    
    # Show plot
    plt.show()
    
    # -----------------------------
    # 6c. Print summary
    # -----------------------------
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"Tangent ratio (ε² constant): ~{df['tangent_ratio'].mean():.3f}")
    print(f"Non-tangent ratio (ε constant): ~{df['nontangent_ratio'].mean():.3f}")
    print(f"Tangent advantage at ε=1e-6: ~{df['nontangent_residual'].iloc[0] / df['tangent_residual'].iloc[0]:.0f}x")
    print("="*70)
    print("\nEXPERIMENT COMPLETE!")
