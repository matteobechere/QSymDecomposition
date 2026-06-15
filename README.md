# QSymDecomposition

The package QSymDecomposition provides tools to perform $q$-Symmetric decompositions of symmetric tensors.

## Installation

To install the package, in Julia, enter the package manager by pressing `]` and then type
`add https://github.com/matteobechere/QSymDecomposition`.

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

