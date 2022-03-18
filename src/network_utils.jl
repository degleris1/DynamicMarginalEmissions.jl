
"""
    make_pfdf_matrix(A, β)

Convert incidence matrix `A` and edge susceptances `β` to the (two-sided) 
PFDF matrix `F`.

Note: the reference node is always set to 1.
"""
function make_pfdf_matrix(A, β)
    n, m = size(A)

    # Make susceptance matrices
    S_l = diagm(β) * A'
    S_b = A * S_l

    # Construct reactance matrix (psuedo-inverse)
    X_b = [
        zeros(n)'; 
        zeros(n-1) inv(S_b[2:n, 2:n])
    ]

    F = S_l * X_b

    return [F; -F]
end

"""
    Fn = make_pfdf_matrix(A, β, nodes)

Convert incidence matrix `A` and edge susceptances `β` to the (two-sided) 
the device-to-line PFDF matrix `Fn`.

Note: the reference node is always set to 1.
"""
function make_device_pfdf_matrix(A, β, nodes)
    F = make_pfdf_matrix(A, β)
    return F[:, nodes]
end