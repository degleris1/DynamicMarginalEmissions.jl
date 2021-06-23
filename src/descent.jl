# Code for training gradient based algorithms

function loss_and_grad(case_list, f̂q, f̂l, M, pmax, gmax, A, B)
    n, m = size(A)

    L = 0.0
    dfq, dfl = zeros(n), zeros(n)
    
    T = length(case_list)
    ∇L = zeros(kkt_dims(n, m))

    for case in case_list
        ĝ, opf, params = solve_pmp(f̂q, f̂l, case.d, pmax, gmax, A, B) 

        # Loss is averaged over cases AND generators, in order
        # to have a rough sense of error per generator
        L += (1/2) * norm(M*ĝ - case.g)^2 / (T*n)
        
        ∇L[1:n] .= M' * (M*ĝ - case.g)
        dfq_case, dfl_case = sensitivity_price(opf, ∇L, params...)
        dfq .+= dfq_case / T
        dfl .+= dfl_case / T
    end
    
    return L, dfq, dfl
end

function stochastic_loss_and_grad(sample, case_list, f̂q, f̂l, M, pmax, gmax, A, B)
    return loss_and_grad(case_list[sample], f̂q, f̂l, M, pmax, gmax, A, B)
end

function solve_pmp(f̂q, f̂l, d, pmax, gmax, A, B)
        params = (f̂q, f̂l, d, pmax, gmax, A, B)
        opf = PowerManagementProblem(params...)
        solve!(opf, () -> ECOS.Optimizer(verbose=false), verbose=false)
        ĝ = evaluate(opf.g)
        return ĝ, opf, params
end

function run_gradient_descent(
    f̂, train_cases, test_cases, pmax, gmax, A, B, M, 
    max_iter, batch_size, test_batch_size, α
)
    train_loss_hist = []
    test_loss_hist = []
    grad_hist = []

    @time for _ in 1:max_iter
        print(".")

        # Evaluate loss and gradient
        sample = rand(1:length(train_cases), batch_size)
        L, df = stochastic_loss_and_grad(sample, train_cases, f̂q, f̂l, M, pmax, gmax, A, B)
        
        # Compute test loss
        sample = rand(1:length(test_cases), test_batch_size)
        L_test, _ = stochastic_loss_and_grad(sample, test_cases, f̂q, f̂l, M, pmax, gmax, A, B)
        
        push!(train_loss_hist, L)
        push!(grad_hist, df)
        push!(test_loss_hist, L_test)
        
        # Take projected gradient step
        η = min(norm(df), α) / norm(df)
        f̂ .= max.(f̂ - η*df, 0)
    end

    println("\nCompleted $(max_iter) iterations.")
    return train_loss_hist, test_loss_hist, grad_hist, f̂
end