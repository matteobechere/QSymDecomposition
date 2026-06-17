using DynamicPolynomials
using LinearAlgebra
using Random
using Statistics
using Printf
using QSymDecomposition
using TensorDec

@polyvar X[1:6] # the amount of X variables must match the length of the tuple q (in this example, q has length 6 and there are 6 variables X[1],...,X[6]
@polyvar Z[1:4]

q = [
    Z[1]^2,
    Z[1]*Z[2],
    Z[1]*Z[3],
    Z[2]^2,
    Z[2]*Z[3],
    Z[3]^2
]

################################################################################
# Computes multinomial coefficients
################################################################################

function get_multinomial_coeff(mon, vars, total_deg)
    den = prod(factorial(degree(mon, v)) for v in vars)
    return factorial(total_deg) / den
end

################################################################################
# Builds a form from a Waring decomposition, given the points Xi and the weights w
################################################################################
function build_polynomial(Xi::AbstractMatrix, w::AbstractVector, d::Int, vars)
    p = zero(vars[1])
    for i in 1:size(Xi, 2)
        lin = sum(Xi[j,i]*vars[j] for j in 1:length(vars))
        p += w[i] * lin^d
    end
    return p
end

################################################################################
# Builds a form from a q-Sym decomposition, given the vector q, the evaluation points Xi, and the weights w
################################################################################
function map_decomp_through_q(w, Xi, q_vec, Z_vars, X_vars, k)
    n_points = size(Xi, 2)
    n_q = length(q_vec)

    q_eval_matrix = zeros(Float64, n_q, n_points)
    for i in 1:n_points
        pt = Xi[:, i]
        for j in 1:n_q
            q_eval_matrix[j, i] = q_vec[j](Z_vars => pt)
        end
    end

    return build_polynomial(q_eval_matrix, w, k, X_vars)
end


################################################################################
# Multi-trial experiment comparing direct decompositions and q-Symmetric decompositions
################################################################################
function run_qsym_experiment(; r, k, trials, threshold=1e-4)
    residuals = Float64[]
    residuals_direct = Float64[]

    for t in 1:trials
        println("\n--- Trial $t ---")
        Xi = randn(length(Z), r)
        for j in 1:r
            Xi[:, j] /= norm(Xi[:, j])
        end
        w  = randn(r)

        # forward map
        p_true = map_decomp_through_q(w, Xi, q, Z, X, k)

        try
            # decomposition
            p_rec = qsym_decompose(p_true, q, k)[:p_decomposition]

            # error
            res = apolar_norm(p_true - p_rec)
            push!(residuals, res)
            println(" Apolar residual: $res")

        catch e
            println("  Failed: $e")
            push!(residuals, NaN)
        end
        try
            # direct decomposition
            w_direct, Xi_direct = decompose(p_true) # specify the rank choice here for better accuracy (the rank MUST coincide with the choice of r at line 135)
            p_direct = build_polynomial(Xi_direct, w_direct, k, X)
            
            res_direct = apolar_norm(p_true - p_direct)
            push!(residuals_direct, res_direct)
            println(" Apolar direct residual: $res_direct")
        catch e
            println("  Failed: $e")
            push!(residuals_direct, NaN)
        end
    end

    clean = filter(!isnan, residuals)
    clean_direct = filter(!isnan, residuals_direct)

    println("\n--- Summary ---")
    println("Trials: ", trials)
    println("Valid runs: ", length(clean))
    println("Valid direct runs: ", length(clean_direct))
    println("Success threshold: ", threshold)

    n_success = count(r -> r < threshold, clean)
    n_success_direct = count(r -> r < threshold, clean_direct)


    println("Success rate: $n_success / $(length(clean))")
    @printf("Mean residual (all): %.5f\n", mean(clean))
    @printf("Mean residual (successes): %.5f\n", mean(filter(r -> r < threshold, clean)))
    @printf("Mean residual (failures): %.5f\n", mean(filter(r -> r >= threshold, clean)))

    println("Success rate direct method: $n_success_direct / $(length(clean_direct))")
    @printf("Mean residual (all): %.5f\n", mean(clean_direct))
    @printf("Mean residual (successes): %.5f\n", mean(filter(r -> r < threshold, clean_direct)))
    @printf("Mean residual (failures): %.5f\n", mean(filter(r -> r >= threshold, clean_direct)))
    return residuals
end

################################################################################
# Run the experiment. Change the parameters r, k and trials accordingly
################################################################################
results = run_qsym_experiment(r=6, k=4, trials=20);
