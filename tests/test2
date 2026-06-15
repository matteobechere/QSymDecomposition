using DynamicPolynomials
using QSymDecomposition

@polyvar X[1:4]
@polyvar Z[1:2]

q = [
    (-Z[1]+Z[2])*(-2*Z[1]+Z[2])*(3*Z[1]+Z[2]),
    Z[1]^2*Z[2],
    Z[1]*Z[2]^2,
    Z[2]^3
]

# The following polynomial is q-Symmetric. qsym_decompose should output an accurate decomposition

p = X[1]^4 + 
    (21*X[1] + 4*X[2] + 2*X[3] + X[4])^4 / 194481 + 
    (100*X[1] + 9*X[2] + 3*X[3] + X[4])^4 / 100000000 + 
    (273*X[1] + 16*X[2] + 4*X[3] + X[4])^4 / 5554571841

result = qsym_decompose(p, q, 4)
