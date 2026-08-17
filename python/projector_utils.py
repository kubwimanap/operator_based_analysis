"""
Projector perturbation utilities for Grassmann manifold analysis.
"""

import numpy as np

def random_projector(n, p):
    """Generate random rank-p orthogonal projector."""
    Q, _ = np.linalg.qr(np.random.randn(n, p))
    return Q @ Q.T

def tangent_perturbation(P, p, n, epsilon):
    """Generate tangent perturbation scaled to Frobenius norm epsilon."""
    Delta = np.random.randn(n, n)
    E = P @ Delta @ (np.eye(n) - P) + (np.eye(n) - P) @ Delta.T @ P
    E = E * (epsilon / np.linalg.norm(E, 'fro'))
    return E

def nontangent_perturbation(n, epsilon):
    """Generate symmetric non-tangent perturbation."""
    E = np.random.randn(n, n)
    E = (E + E.T) / 2
    E = E * (epsilon / np.linalg.norm(E, 'fro'))
    return E

def residual(P_tilde):
    """Compute residual norm ||P̃² - P̃||_F."""
    return np.linalg.norm(P_tilde @ P_tilde - P_tilde, 'fro')

def backward_error(P_tilde, p):
    """Compute backward error to nearest rank-p projector."""
    eigvals, eigvecs = np.linalg.eigh(P_tilde)
    idx = np.argsort(eigvals)[::-1][:p]
    U = eigvecs[:, idx]
    P_star = U @ U.T
    return np.linalg.norm(P_tilde - P_star, 'fro')
