

using JuMP 
using Gurobi


P = ["pv", "gas", "wind", "coal"]
N = ["N1", "N2", "N3"]
L = ["L1", "L2", "L3"]

mc = Dict(k => v for (k,v) in zip(P, [0 50 0 30]))
gmax = Dict(k => v for (k,v) in zip(P, [80 150 120 300]))
map_np = Dict(k => v for (k,v) in zip(P, ["N1" "N1" "N2" "N3"]))
map_pn = Dict(k => v for (k,v) in zip(N, [["pv", "gas"],["wind"],["coal"]]))
demand =  Dict(k => v for (k,v) in zip(N, [150 30 200]))

slack = "N3"
PTDF = [
    -0.4 0.2 0;
    0.6 0.2 0;
    0.4 0.8 0]

    



fmax = [40 100 100]


# Dispatch 
dispatch = Model(with_optimizer(Gurobi.Optimizer))
@variable(dispatch, 0 <= G[p=P] <= gmax[p])
@variable(dispatch, I[n=N])


@objective(dispatch, Min, sum(G[p] * mc[p] for p in P));
@constraint(dispatch, EnergyBalance[n=N], demand[n] - sum(G[p] for p in map_np[n]) == I[n]);
@constraint(dispatch, Slack, sum(I) == 0);

optimize!(dispatch)
objective_value(dispatch)
value.(G)
dual.(Slack)

# OPF
@constraint(dispatch, Flow, -fmax' .<= PTDF * [I[n] for n in N]  .<= fmax');
optimize!(dispatch)
objective_value(dispatch)
value.(G)
dual.(Flow)
PTDF * [value(I[n]) for n in N]



# ADMM Generators
function subproblem(g, G_prev, lambda)
    sub = Model(with_optimizer(Gurobi.Optimizer))
    set_silent(sub)
    @variable(sub, 0 <= G <= gmax[g]);
    @expression(sub, penalty_term, 
        ((sum(values(demand)) - G - sum(G_prev[p] for p in setdiff(P, [g])))^2)
    );
    @objective(sub, Min, mc[g] * G - lambda * G + gamma/2 * penalty_term)
    optimize!(sub)
    return value.(G)
end


G_prev = Dict(k => v for (k,v) in zip(P, [0. 0. 90. 100.]))
# G_prev = Dict(k => v for (k,v) in zip(P, [80. 0. 120. 180.]))
lambda = [10.]
gamma = 0.2
for i in 1:100
    # i = 1
    G_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
    println("Iteration $(i) with lambda $(lambda[end])")
    for p in P
        G_new[p] = subproblem(p, G_prev, lambda[end])
        print("$(p): old  $(G_prev[p]) - new $(G_new[p]). ")
    end
    print("\n")
    G_prev = G_new
    lambda_new = lambda[end] + 0.5 * gamma * (sum(values(demand)) - sum(values(G_new)))
    push!(lambda, lambda_new)
end


# ADMM Nodes
ptdf = Dict(k => v for (k,v) in zip(N, PTDF[1, :]))

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
