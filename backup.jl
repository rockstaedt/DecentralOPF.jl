

# ADMM Nodes
T = Dict(k => v for (k,v) in zip(N, ptdf[1, :]))

function subproblem(n, I_prev, lambda, mu)
    sub = Model(with_optimizer(Gurobi.Optimizer))
    set_silent(sub)
    Pn = map_pn[n]
    @variable(sub, 0 <= G[p=Pn] <= gmax[p]);
    @variable(sub, I);

    @expression(sub, penalty_term, 
        (I + sum(I_prev[n] for n in setdiff(N, [n])))^2
    );

    @objective(sub, Min, 
        sum(mc[p] * G[p] for p in Pn) 
        -sum(lambda*G) - mu * ptdf[n] * I + gamma/2 * penalty_term)

    @constraint(sub, I == sum(G[p] for p in Pn) - demand[n] )
    optimize!(sub)
    println(value.(penalty_term))
    return value.(I)
end


I_prev = Dict(k => v for (k,v) in zip(N, [-50. 90. -20.]))
# G_prev = Dict(k => v for (k,v) in zip(P, [80. 0. 120. 180.]))
lambda = [10.]
mu = [10.]
gamma = 0.3
for i in 1:100
    # i = 1
    I_new = Dict(k => v for (k,v) in zip(N, [0. 0. 0.]))
    println("Iteration $(i) with lambda $(lambda[end])")
    for n in N
        I_new[n] = subproblem(n, I_prev, lambda[end], mu[end])
        print("$(n): old  $(I_prev[n]) - new $(I_new[n]). ")
    end
    print("\n")
    I_prev = I_new
    lambda_new = lambda[end] - 0.5 * gamma * (sum(values(I_new)))
    mu_new = mu[end] - 0.5 * gamma * (40 - sum(ptdf[n] * I_new[n] for n in N))
    push!(lambda, lambda_new)
    push!(mu, mu_new)
end

I_new
lambda
mu
PTDF * [value(I_new[n]) for n in N]

