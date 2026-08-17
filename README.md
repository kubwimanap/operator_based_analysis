# OPERATOR-BASED ANALYSIS OF NUMERICAL PROJECTION PERTURBATIONS IN GRASSMANN MANIFOLD OPTIMIZATION

**Authors:** Pascal Kubwimana, B. Prabhakar Reddy, Jason M. Mkenyeleye  
**Affiliation:** University of Dodoma  
**Contact:** kubwimanap4@gmail.com

---

## Overview

Numerical validation of quadratic residual scaling for tangent-structured perturbations on Grassmann manifolds.

### Key Results

- **Quadratic Identity:** P̃² - P̃ = E²
- **Sharp Bound:** ||P̃² - P̃||_F ≤ ε²/√2
- **Tangent scaling:** O(ε²)
- **Non-tangent scaling:** O(ε)

---

## Experimental Results

| ε | Tangent Residual | Non-tangent Residual | Tangent Ratio |
|---|---|---|---|
| 1.000e-06 | 2.361e-13 | 9.068e-07 | 0.236 |
| 1.000e-05 | 2.370e-11 | 9.060e-06 | 0.237 |
| 1.000e-04 | 2.370e-09 | 9.066e-05 | 0.237 |
| 1.000e-03 | 2.372e-07 | 9.076e-04 | 0.237 |
| 1.000e-02 | 2.370e-05 | 9.068e-03 | 0.237 |
| 1.000e-01 | 2.368e-03 | 9.066e-02 | 0.237 |

**Conclusion:** Tangent perturbations are ~3,800× more accurate at ε = 10⁻⁶.

---

## Rank Dependence (ε = 10⁻⁶)

| Rank p | Backward Error | Ratio |
|--------|---------------|-------|
| 1 | 7.100e-13 | 0.707 |
| 5 | 3.200e-13 | 0.325 |
| 10 | 2.400e-13 | 0.237 |
| 20 | 1.800e-13 | 0.178 |
| 40 | 1.500e-13 | 0.145 |

---

## Installation

```bash
git clone https://github.com/kubwimanap/operator_based_analysis.git
cd operator_based_analysis
pip install -r python/requirements.txt


operator_based_analysis/
├── python/
│   ├── experiments.py
│   ├── projector_utils.py
│   ├── visualization.py
│   └── requirements.txt
├── results/
│   └── grassmann_experiment_results.csv
├── Synthetic_matlab/
├── synthetic_dataset/
└── README.md


@article{kubwimana2026operator,
  title={Operator-Based Analysis of Numerical Projection Perturbations in Grassmann Manifold Optimization},
  author={Kubwimana, Pascal and Reddy, B. Prabhakar and Mkenyeleye, Jason M.},
  year={2026}
}

