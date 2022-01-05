

using JuMP 
using Gurobi
using DataFrames
using Mosek 
using MosekTools
include("utility.jl")

# %%
begin 
    plants = DataFrame(
        P = ["pv", "gas", "wind", "coal"],
        mc = [0, 50, 0, 30],
        gmax = [80, 150, 120, 300],
        node = ["N1", "N1", "N2", "N3"]
    );

    nodes = DataFrame(
        N = ["N1", "N2", "N3"],
        demand = [150, 30, 200],
        slack = [false, false, true]
    );

    lines = DataFrame(
        L = ["L1", "L2", "L3"],
        from = ["N2", "N3", "N2"],
        to = ["N1", "N1", "N3"],
        fmax = [40, 100, 100],
        susceptance = [1, 1, 2]

    );

    P = plants[:, :P]
    N = nodes[:, :N]
    L = lines[:, :L]

    mc = Dict(k => v for (k,v) in eachrow(plants[:, [:P, :mc]]))
    gmax = Dict(k => v for (k,v) in  eachrow(plants[:, [:P, :gmax]]))
    map_np = Dict(k => v for (k,v) in eachrow(plants[:, [:P, :node]]))
    map_pn = Dict(n => filter(row -> row.node == n, plants)[:, :P] for n in N)
    demand =  Dict(k => v for (k,v) in eachrow(nodes[:, [:N, :demand]]))

    ptdf = calculate_ptdf(nodes, lines)
end

# Dispatch 
dispatch = Model(with_optimizer(Gurobi.Optimizer))
@variable(dispatch, 0 <= G[p=P] <= gmax[p])
@variable(dispatch, I[n=N])


@objective(dispatch, Min, sum(G[p] * mc[p] for p in P));
@constraint(dispatch, EnergyBalance[n=N], sum(G[p] for p in map_pn[n]) - demand[n] == I[n]);
@constraint(dispatch, Slack, sum(I) == 0);

optimize!(dispatch)
objective_value(dispatch)
value.(G)

plants[!, :G_dispatch] = value.(G).data
lines[!, :F_dispatch] = ptdf * [value(I[n]) for n in N]

# OPF
@constraint(dispatch, Flow, -lines[:, :fmax] .<= ptdf * [I[n] for n in N]  .<= lines[:, :fmax]);
optimize!(dispatch)
objective_value(dispatch)

lines[!, :F_opf] = ptdf * [value(I[n]) for n in N]
plants[!, :G_opf] = value.(G).data
# lambda_ref = dual.(EnergyBalance)

# ADMM Generators
function subproblem(g, G_prev, lambda)
    sub = Model(with_optimizer(Gurobi.Optimizer))
    set_silent(sub)
    @variable(sub, 0 <= G <= gmax[g]);
    @expression(sub, penalty_term, 
        ( G + sum(G_prev[p] for p in setdiff(P, [g])) - sum(values(demand)))^2
    );
    @objective(sub, Min, mc[g] * G - lambda * G + gamma/2 * penalty_term)
    optimize!(sub)
    return value.(G)
end

function subproblem(g, G_prev, lambda)
    sub = Model(with_optimizer(Mosek.Optimizer))
    set_silent(sub)
    @variable(sub, 0 <= G <= gmax[g]);
    @variable(sub, T >= 0)

    @constraint(sub, penalty_term, 
        [T, sum(values(demand)) - G - sum(G_prev[p] for p in setdiff(P, [g]))] in SecondOrderCone()
    );
    @objective(sub, Min, mc[g] * G - lambda * G + gamma/2 * T);

    optimize!(sub)
    return value.(G)
end

begin
    G_prev = Dict(k => v for (k,v) in zip(P, [80. 0. 100. 100.]))
    G_prev = Dict(k => v for (k,v) in zip(P, [80. 0. 120. 180.]))
    lambda = [10.]
    gamma = 0.3
    rho = 0.1
    for i in 1:200
        G_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
        for p in P
            G_new[p] = subproblem(p, G_prev, lambda[end])
        end
        G_prev = G_new
        lambda_new = lambda[end] + rho * gamma * (sum(values(demand)) - sum(values(G_new)))
        push!(lambda, lambda_new)
    end
    result = copy(plants)
    result[!, :G_admm] = [G_prev[p] for p in P]
    println(result)
    println("Total Demand: ", sum(values(demand)))
    println("Total Generation: ", sum(values(G_prev)))
end

lambda
plants
G_prev



# ADMM OPF Generators
T = Dict(k => v for (k,v) in zip(N, ptdf[1, :]))
fmax = 40
function subproblem(g, G_prev, R_ref_prev, R_cref_prev, lambda, mu_ref, mu_cref)
    sub = Model(with_optimizer(Gurobi.Optimizer))
    set_silent(sub)
    @variable(sub, 0 <= G <= gmax[g]);
    @variable(sub, 0 <= R_ref);
    @variable(sub, 0 <= R_cref);
    @expression(sub, I[n=N], 
        demand[n] - sum(G_prev[p] for p in setdiff(map_pn[n], [g])) - (map_np[g] == n ? G : 0)
    );

    @expression(sub, penalty_term_eb, 
        ( G + sum(G_prev[p] for p in setdiff(P, [g])) - sum(values(demand)) - lambda/rho )^2
    );
    @expression(sub, penalty_term_flow_ref, 
        (sum(T[n] * I[n] for n in N) + R_ref - fmax - mu_ref/rho )^2
    );
    @expression(sub, penalty_term_flow_cref, 
        (R_cref - sum(T[n] * I[n] for n in N) - fmax - mu_cref/rho )^2
    );
    @expression(sub, penalty_term_cref, 
        (R_cref - R_cref_prev)^2
    );
    @expression(sub, penalty_term_ref, 
        (R_ref - R_ref_prev)^2
    );
    @expression(sub, penalty_term_g, 
        (G - G_prev[g])^2
);

    @objective(sub, Min, mc[g] * G 
        + rho/2 * penalty_term_eb 
        + rho/2 * penalty_term_flow_ref
        + rho/2 * penalty_term_flow_cref
        + 0.5 * penalty_term_g
        + 0.5 * penalty_term_ref
        + 0.5 * penalty_term_cref
        )

    optimize!(sub)
    return value.(G), value.(R_ref), value.(R_cref)
end

function net_injection(G)
    return Dict(n => demand[n] - sum(G[p] for p in map_pn[n]) for n in N)
end


begin 
    G_prev = Dict(k => v for (k,v) in zip(P, [50. 0. 120. 165.]))

    R_ref_avg = 0 
    R_cref_avg = 0
    # G_prev = Dict(k => v for (k,v) in zip(P, [80. 0. 120. 180.]))
    lambda = [10.]
    mu_ref = [1.]
    mu_cref = [-1.]
    gamma = 0.4
    rho = 0.8

    for i in 1:200
        # i = 1
        G_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
        R_ref_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
        R_cref_new = Dict(k => v for (k,v) in zip(P, [0. 0. 0. 0.]))
        for p in P
            G_new[p], R_ref_new[p], R_cref_new[p] = subproblem(
                p, G_prev, R_ref_avg, R_cref_avg, 
                lambda[end], mu_ref[end], mu_cref[end])
        end
        G_prev = G_new
        lambda_new = lambda[end] + rho * gamma * (sum(values(demand)) - sum(values(G_new)))
        push!(lambda, lambda_new)

        R_ref_avg = sum(R_ref_new[p] for p in P)/length(P)
        R_cref_avg = sum(R_cref_new[p] for p in P)/length(P)
        
        R_ref_prev = R_ref_new
        R_cref_prev = R_cref_new

        mu_ref_new = mu_ref[end] - gamma * rho * (sum(T[n] * net_injection(G_prev)[n] for n in N) + R_ref_avg - fmax)
        push!(mu_ref, mu_ref_new)
        mu_cref_new = mu_cref[end] - gamma * rho * (R_cref_avg - sum(T[n] * net_injection(G_prev)[n] for n in N) - fmax)
        push!(mu_cref, mu_cref_new)
    end
    result = copy(plants)
    result[!, :G_admm] = [G_new[p] for p in P]
    println(result)
    println("Total Demand: ", sum(values(demand)))
    println("Total Generation: ", sum(values(G_prev)))
end

mu_ref_new
mu_cref_new

R_cref_prev
R_ref_prev

obj = sum(mc[p]*G_prev[p] for p in P)

plants
G_prev
lambda