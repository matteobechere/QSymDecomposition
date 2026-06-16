export qsym_decompose

#######################################################################################################
# Function taking as input the polynomial p, the vector of polynomials q=[], the degree k and tolerance
#######################################################################################################
function qsym_decompose(p, q::AbstractVector, k::Int; Z_vars = unique(reduce(vcat, variables.(q))), tol=1e-8)
    X_vars = variables(p)

    ######################################################################################################################################
    # build_polynomial is to construct a polynomial starting from a matrix of nodes Xi, the weights w, the degree d and the variables vars
    # the output polynomial is sum_{j=1}^{size(Xi, 2)} w[j](sum_{i=1}^{size(Xi,1)} Xi[j,i]vars[i])^d
    ######################################################################################################################################
    function build_polynomial(Xi::AbstractMatrix, w::AbstractVector, d::Int, vars)
        n_terms = size(Xi, 2)
        p_reconst = zero(vars[1])
        for i in 1:n_terms
            linform = sum(Xi[j,i]*vars[j] for j in 1:length(vars))
            p_reconst += w[i] * linform^d
        end
        return p_reconst
    end
    ########################################################################################################################################################
    # Wq_matrix takes as input the variables, the degrees and the vector of polynomials q. It outputs the matrix associated to W_q w.r.t. the monomial bases
    # W_q is the substitution map s.t. X^alpha |-> q^alpha induced by X_i |-> q_i
    ########################################################################################################################################################
    function Wq_matrix(X_vars, Z_vars, k, h, q_vec)
        Bcod = reverse(monomials(Z_vars, h*k))
        Bdom = reverse(monomials(X_vars, k))
        cod_index = Dict(Bcod[i] => i for i in eachindex(Bcod))
        M = zeros(Float64, length(Bcod), length(Bdom))
        for (j, mon) in enumerate(Bdom)
            p_comp = prod(q_vec[i]^degree(mon, X_vars[i]) for i in 1:length(X_vars)) + zero(q_vec[1])
            for t in terms(p_comp)
                m = monomial(t)
                if haskey(cod_index, m)
                    M[cod_index[m], j] += coefficient(t)
                end
            end
        end
        return M, Bdom, Bcod
    end

    ################################################################################################################
    # get_multinomial_coeff computes multnomial coefficients given the total degree and a monomial in some variables
    ################################################################################################################
    function get_multinomial_coeff(mon, vars, total_deg)
        den = prod(factorial(degree(mon, v)) for v in vars)
        return factorial(total_deg) / den
    end

    ############################################################
    # q-Symmetric decomposition
    ############################################################
    println("Trying q-Symmetric decomposition...")

    h = maximum(maxdegree.(q))
    M, Bdom, Bcod = Wq_matrix(X_vars, Z_vars, k, h, q) # Bdom and Bcod are the bases of the domain and codomain respectively

    hk = Int(maxdegree(Bcod[1]))
    d_hk_inv = [1.0 / get_multinomial_coeff(m, Z_vars, hk) for m in Bcod] # diagonal matrix for normalization purposes (see definition of phi_q)

    k_deg = Int(maxdegree(Bdom[1]))
    d_k = [Float64(get_multinomial_coeff(m, X_vars, k_deg)) for m in Bdom] # diagonal matrix for normalization purposes (see definition of phi_q)

    A = M'
    A = d_k .* A
    A = A .* d_hk_inv'

    dom_index = Dict(Bdom[i] => i for i in eachindex(Bdom))
    v = zeros(Float64, length(Bdom))
    for t in terms(p)
        m = monomial(t)
        if haskey(dom_index, m)
            v[dom_index[m]] = coefficient(t)
        end
    end
    # check if the polynomial p is a q-Sym polynomial
    x = A \ v
    residual = norm(A * x - v)
    if residual > tol
        @warn "Polynomial is not a q-Symmetric polynomial. Residual: $residual"
    end

    # reconstruct psi(p)
    poly_z = zero(sum(Z_vars))
    for i in eachindex(x)
        if abs(x[i]) > tol # eliminates terms which are close to zero
            poly_z += x[i] * Bcod[i]
        end
    end

    wPsi, XiPsi = decompose(poly_z) # Waring decomposition of psi(p)
    n_points = size(XiPsi, 2)
    n_q = length(q)

    q_eval_matrix = zeros(eltype(XiPsi), n_q, n_points) # initializing the matrix containing the nodes of the q-Sym decomposition

    for i in 1:n_points
        pt = XiPsi[:, i]  # point in Z-space

        for j in 1:n_q
            q_eval_matrix[j, i] = q[j](Z_vars => pt)
        end
    end

    # now reconstruct the q-Sym decomposition in X variables
    p_dec = build_polynomial(q_eval_matrix, wPsi, k, X_vars)
    q_residual = apolar_norm(p - p_dec)

    return Dict(
        :weights => wPsi,
        :points => XiPsi,
        :qSympoints => q_eval_matrix,
        :psi_polynomial => poly_z,
        :p_decomposition => p_dec,
        :residual => q_residual
    )
end
