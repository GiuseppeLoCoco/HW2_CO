# MATH-329: Continuous Optimization — Homework 2 (Group H)

**Group Members:**
* Timothée Berthon
* Giuseppe Lo Coco
* Elhadji Ousmane Ndiaye
* Jean-Baptiste Rombaldi


## Directory Structure

```text
H/
├── H.pdf               # Compiled report
├── H.tex               # LaTeX source of the report
├── figures             # Files for compilation of the LaTeX document
├── README.md           # Execution instructions and runtime overview
├── DISCLAIMER.md       # Tool usage and collaboration disclosure
├── code/               # MATLAB source code
│   ├── main.m                  # Main entry point running all experiments and unit tests
│   ├── tCG.m                   # Truncated Conjugate Gradient (tCG) subproblem solver
│   ├── TR_tCG.m                # Trust-Region framework with tCG
│   ├── GD_linesearch.m         # Gradient Descent with Armijo backtracking line search
│   ├── phase_retrieval_f.m     # Real Phase Retrieval objective function
│   ├── phase_retrieval_grad.m  # Real Phase Retrieval analytical gradient
│   ├── phase_retrieval_hess.m  # Real Phase Retrieval Hessian-vector product operator
│   ├── rosenbrockf.m           # Rosenbrock benchmark objective
│   ├── rosenbrockgrad.m        # Rosenbrock benchmark gradient
│   └── rosenbrockhess.m        # Rosenbrock benchmark Hessian-vector product
└── results/            # Output figures generated automatically by main.m