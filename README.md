# QSymDecomposition.jl

The package QSymDecomposition.jl provides tools to perform $q$-Symmetric decompositions of symmetric tensors. The function `qsym_decompose` implements Algorithm 1 from [Symmetric tensor decomposition on rational varieties](https://arxiv.org/abs/2606.25712) to Julia.
The `tests` folder of this repository contains Julia codes to reproduce the results in [Section 4. Effective decompositions of q-Symmetric tensors](https://arxiv.org/abs/2606.25712).

## Installation

To install the package, in Julia, enter the package manager by pressing `]` and then type

`add https://github.com/matteobechere/QSymDecomposition.jl`.

## Example usage, within Julia

```
using DynamicPolynomials
using QSymDecomposition

@polyvar X[1:9]
@polyvar Z[1:3]

q = [
    Z[1]^3,
    Z[1]^2*Z[2],
    Z[1]^2*Z[3],
    Z[1]*Z[2]^2,
    Z[1]*Z[3]^2,
    Z[2]^3,
    Z[2]^2*Z[3],
    Z[2]*Z[3]^2,
    Z[3]^3
]

p = (X[1] + X[2] + X[3] + X[4] + X[5] + X[6] + X[7] + X[8] + X[9])^3 +
    (X[1] - X[2] + X[3] + X[4] + X[5] - X[6] + X[7] - X[8] + X[9])^3 +
    (X[1] - X[2] - X[3] + X[4] + X[5] - X[6] - X[7] - X[8] - X[9])^3 +
    (8*X[1] + 4*X[2] + 4*X[3] + 2*X[4] + 2*X[5] + X[6] + X[7] + X[8] + X[9])^3 +
    (X[1] + 3*X[2] + X[3] + 9*X[4] + X[5] + 27*X[6] + 9*X[7] + 3*X[8] + X[9])^3 +
    (8*X[1] + 4*X[2] + 8*X[3] + 2*X[4] + 8*X[5] + X[6] + 2*X[7] + 4*X[8] + 8*X[9])^3 +
    (- 8*X[1] + 4*X[2] + 8*X[3] - 2*X[4] - 8*X[5] + X[6] + 2*X[7] + 4*X[8] + 8*X[9])^3

result = qsym_decompose(p, q, 3)
```

A $q$-Symmetric decomposition looks like $p_{\text{dec}}=\sum_{i=1}^r\omega_i(q_0(\xi_i)X_0+\cdots+q_n(\xi_i)X_n)^k$.

The function `qsym_decompose` outputs a Dict called `result`, containing
- `psi_polynomial` is the image of $p$ under the $\psi_q$ map;
- `weights` is a tuple containing the scalars $\omega_1,\ldots,\omega_r$;
- `points` is a matrix whose columns are the tuples $\xi_1,\ldots,\xi_r$;
- `qSympoints` is a matrix whose i-th column is the tuples $q_0(\xi_i),\ldots,q_n(\xi_i)$;
- `p_decomposition` is the expanded $p_{\text{dec}}$, i.e. $q$-Symmetric decomposition of $p$;
- `residual` is the (apolar) norm of the difference $p-p_{\text{dec}}$. A small `residual` ensures the decomposition was successful.
