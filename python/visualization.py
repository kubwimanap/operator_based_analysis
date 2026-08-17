"""
Publication-quality visualization for Grassmann projector analysis.
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

# Set publication-quality style
plt.style.use('seaborn-v0_8-paper')
rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'legend.fontsize': 10,
    'figure.dpi': 300,
    'savefig.dpi': 300,
})

def plot_residual_scaling_publication(df, save_path="../results/figures/residual_scaling.pdf"):
    """
    Publication-quality residual scaling plot.
    """
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # Data
    ax.loglog(df["epsilon"], df["tangent_residual"], 'o-', 
              color='#2E86AB', linewidth=2.5, markersize=8,
              label="Tangent residual (O(ε²))")
    ax.loglog(df["epsilon"], df["nontangent_residual"], 's-', 
              color='#A23B72', linewidth=2.5, markersize=8,
              label="Non-tangent residual (O(ε))")
    
    # References
    ax.loglog(df["epsilon"], df["epsilon"]**2, 'k--', 
              linewidth=1.5, alpha=0.5, label="ε² reference")
    ax.loglog(df["epsilon"], df["epsilon"], 'r--', 
              linewidth=1.5, alpha=0.5, label="ε reference")
    
    # Labels and title
    ax.set_xlabel("Perturbation norm ε", fontsize=14)
    ax.set_ylabel("Residual $\\Vert \\widetilde{P}^2 - \\widetilde{P} \\Vert_F$", fontsize=14)
    ax.set_title("Quadratic vs Linear Scaling of Residuals", fontsize=16, pad=15)
    
    # Legend
    ax.legend(loc='lower right', frameon=True, fancybox=True, shadow=True)
    
    # Grid
    ax.grid(True, which="both", ls="--", alpha=0.3)
    
    # Add annotation
    ax.text(0.02, 0.98, f"n=100, p=10", 
            transform=ax.transAxes, fontsize=11,
            verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    plt.savefig(save_path, bbox_inches='tight')
    plt.show()
