using DynamicPolynomials
using QSymDecomposition

@polyvar X[1:6]
@polyvar Z[1:3]


q = [
    Z[1]^2,
    Z[1]*Z[2],
    Z[1]*Z[3],
    Z[2]^2,
    Z[2]*Z[3],
    Z[3]^2
]

# The following polynomial is q-Symmetric. qsym_decompose should output an accurate decomposition

p = (4*X[1] - 2*X[2] + 2*X[3] + X[4] - X[5] + X[6])^4 -
    (X[1] + 2*X[2] - X[3] + 4*X[4] - 2*X[5] + X[6])^4 +
    (X[1] - X[2] + 2*X[3] + X[4] - 2*X[5] + 4*X[6])^4 +
    (4*X[1] + 2*X[2] - 2*X[3] + X[4] - X[5] + X[6])^4 -
    (X[1] + X[2] + X[3] + X[4] + X[5] + X[6])^4 +
    (X[1] - 2*X[2] - X[3] + 4*X[4] + 2*X[5] + X[6])^4 +
    (4*X[1] + 2*X[2] + 2*X[3] + X[4] + X[5] + X[6])^4

result = qsym_decompose(p, q, 2)
