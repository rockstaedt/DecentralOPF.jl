using JuMP
using Gurobi

include("structures.jl")
include("plotting.jl")
include("penalty_terms.jl")


# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

node1 = Node("N1", [10, 250], false)
node2 = Node("N2", [50, 30], false)
node3 = Node("N3", [50, 150], true)

nodes = [node1, node2, node3]

line1 = Line("L1", node2, node1, 20, 1)
line2 = Line("L2", node3, node1, 30, 1)
line3 = Line("L3", node2, node3, 20, 2)

lines = [line1, line2, line3]

pv = Generator("pv", 3, 80, "yellow", node1)
wind = Generator("wind", 4, 120, "lightblue", node2)
coal = Generator("coal", 30, 300, "brown", node3)
gas = Generator("gas", 50, 120, "grey", node1)

generators = [pv, wind, coal, gas]

battery = Storage("battery", 1, 10, 20, "purple", node1)
ps = Storage("ps", 2, 15, 30, "blue", node1)

storages = [battery]

function optimize_subproblem(generator::Generator)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)

    # Maximum capacity of generator
    @variable(sub, 0 <= G[t=admm.T] <= generator.max_generation)

    # Slack variables
    @variable(sub, 0 <= R_ref[line=admm.L, t=admm.T])
    @variable(sub, 0 <= R_cref[line=admm.L, t=admm.T])


    prev_G = get_unit_results(generator, admm.iteration - 1)

    prev_node_results = get_node_id_to_result(admm.iteration - 1)

    # Injection term
    @expression(
        sub,
        injection[n=admm.N, t=admm.T],
        n == admm.node_to_id[generator.node] ? (
            G[t] + (prev_node_results[n].generation[t] - prev_G[t])
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        ) : (
            prev_node_results[n].generation[t]
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        )
    )

    add_penalty_terms!(sub)

    node_index = admm.node_to_id[generator.node]

    @objective(
        sub,
        Min,
        sum(G[t] * generator.marginal_costs
            + admm.lambdas[admm.iteration][t] * G[t]
            + admm.gamma/2 * sum(
                admm.ptdf[l, node_index] * (
                    admm.mues[admm.iteration][l, t]
                    - admm.rhos[admm.iteration][l, t]
                ) * G[t] for l in admm.L
            )
            + admm.gamma/2 * sub[:penalty_term_eb][t]
            # + sum(admm.mues[admm.iteration][:, t]) * G[t]
            + admm.gamma/2 * sub[:penalty_term_flow1][t]
            # - sum(admm.rhos[admm.iteration][:, t]) * G[t]
            + admm.gamma/2 * sub[:penalty_term_flow2][t]
            + admm.gamma/2 * sub[:penalty_term_R_ref][t]
            + admm.gamma/2 * sub[:penalty_term_R_cref][t]
            # + 1/2 * (G[t] - prev_G[t])^2
            for t in admm.T)
    )

    optimize!(sub)

    penalty_terms = PenaltyTerm(
        value.(sub[:penalty_term_eb]).data,
        value.(sub[:penalty_term_flow1]).data,
        value.(sub[:penalty_term_flow2]).data
    )

    result = ResultGenerator(
        generator,
        value.(G).data,
        penalty_terms,
        value.(R_ref).data,
        value.(R_cref).data
    )

    return result
end

function optimize_subproblem(storage::Storage)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)

    # Storage variables with maximum bounds
    @variable(sub, 0 <= G_S_d[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= G_S_c[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= storage_level[t=admm.T] <= storage.max_level)

    # Slack variables
    @variable(sub, 0 <= R_ref[line=admm.L, t=admm.T])
    @variable(sub, 0 <= R_cref[line=admm.L, t=admm.T])


    previous_S_d, previous_S_c = get_unit_results(storage, admm.iteration - 1)

    prev_node_results = get_node_id_to_result(admm.iteration - 1)

    # Injection term
    @expression(
        sub,
        injection[n=admm.N, t=admm.T],
        n == admm.node_to_id[storage.node] ? (
            prev_node_results[n].generation[t]
            + G_S_d[t] + (prev_node_results[n].discharge[t] - previous_S_d[t])
            - G_S_c[t] - (prev_node_results[n].charge[t] - previous_S_c[t])
            - admm.node_id_to_demand[n][t]
        ) : (
            prev_node_results[n].generation[t]
            + prev_node_results[n].discharge[t]
            - prev_node_results[n].charge[t]
            - admm.node_id_to_demand[n][t]
        )
    )

    add_penalty_terms!(sub)

    # Storage level constraint
    @constraint(
        sub,
        StorageBalance[t=admm.T],
        storage_level[t] == (
            (t == 1 ? 0 : storage_level[t-1]) + G_S_c[t] - G_S_d[t]
        )
    )

    node_index = admm.node_to_id[storage.node]

    @objective(
        sub,
        Min,
        sum(storage.marginal_costs * (G_S_d[t] + G_S_c[t])
            + admm.lambdas[admm.iteration][t] * (G_S_d[t] - G_S_c[t])
            + admm.gamma/2 * sum(
                admm.ptdf[l, node_index] * (
                    admm.mues[admm.iteration][l, t]
                    - admm.rhos[admm.iteration][l, t]
                ) * (G_S_d[t] - G_S_c[t]) for l in admm.L
            ) 
            + admm.gamma/2 * sub[:penalty_term_eb][t]
            # + sum(admm.mues[admm.iteration][:, t]) * (G_S_d[t] - G_S_c[t])
            + admm.gamma/2 * sub[:penalty_term_flow1][t]
            # - sum(admm.rhos[admm.iteration][:, t]) * (G_S_d[t] - G_S_c[t])
            + admm.gamma/2 * sub[:penalty_term_flow2][t]
            + admm.gamma/2 * sub[:penalty_term_R_ref][t]
            + admm.gamma/2 * sub[:penalty_term_R_cref][t]
            # + 1/2 * (G_S_d[t] - previous_S_d[t])^2
            # + 1/2 * (G_S_c[t] - previous_S_c[t])^2
            for t in admm.T)
    )

    optimize!(sub)

    penalty_terms = PenaltyTerm(
        value.(sub[:penalty_term_eb]).data,
        value.(sub[:penalty_term_flow1]).data,
        value.(sub[:penalty_term_flow2]).data
    )

    result = ResultStorage(
        storage,
        value.(G_S_d).data,
        value.(G_S_c).data,
        value.(storage_level).data,
        penalty_terms,
        value.(R_ref).data,
        value.(R_cref).data
    )

    return result
end

function update_lambda()
    lambdas = admm.lambdas[admm.iteration] + admm.gamma * (
        admm.results[admm.iteration].generation
        + admm.results[admm.iteration].discharge
        - admm.results[admm.iteration].charge
        - admm.total_demand
    )
    push!(admm.lambdas, lambdas)
end

function update_mue()
    injection = admm.results[admm.iteration].injection
    avg_R_ref, _ = get_slack_results(admm.iteration)
    mues = (
        admm.mues[admm.iteration]
        + admm.gamma * (admm.ptdf * injection + avg_R_ref .- admm.L_max)
    )
    push!(admm.mues, mues)
end

function update_rho()
    injection = admm.results[admm.iteration].injection
    _, avg_R_cref = get_slack_results(admm.iteration)
    rhos = (
        admm.rhos[admm.iteration]
        + admm.gamma * (avg_R_cref- admm.ptdf * injection .- admm.L_max)
    )
    push!(admm.rhos, rhos)
end

function update_duals()
    update_lambda()
    update_mue()
    update_rho()
end

function check_convergence()
    eps = 10^(-5)
    # if admm.iteration != 1
    #     # Iteration counter is not increased yet. Hence, the new lamba is
    #     # compared to the lambda of the current iteration.
    #     lambda_convergence = (
    #         abs.(admm.lambdas[end] - admm.lambdas[admm.iteration]) .< eps
    #     )
    #     mue_convergence = (
    #         abs.(admm.mues[end] - admm.mues[admm.iteration]) .< eps
    #     )
    #     rho_convergence = (
    #         abs.(admm.rhos[end] - admm.rhos[admm.iteration]) .< eps
    #     )
    #     convergence = hcat(lambda_convergence, mue_convergence)
    #     convergence = hcat(convergence, rho_convergence)
    #     admm.converged = all(convergence)
    # end

    if admm.converged
        println("Converged")
    else
        println("Not converged")
        admm.iteration += 1
    end
end

function calculate_iteration()
    println("\nIteration: $(admm.iteration)")
    println("###############")
    print_duals(admm.iteration)

    if !isnothing(admm.storages)
        unit_to_result = Dict(zip(
            vcat(generators, storages),
            vcat(
                optimize_subproblem.(generators),
                optimize_subproblem.(storages)
            )
        ))
    else
        unit_to_result = Dict(zip(generators, optimize_subproblem.(generators)))
    end

    result = Result(unit_to_result)

    push!(admm.results, result)

end

function print_results(print_generators::Bool,
                       print_storages::Bool,
                       print_penalty::Bool,
                       iteration::Int)
    unit_to_result = admm.results[iteration].unit_to_result

    for (unit, result) in unit_to_result
        if typeof(unit) == Generator && print_generators
            print("$(result.generator.name): ")
            println(result.generation)
        elseif typeof(unit) == Storage && print_storages
            println("$(result.storage.name): ")
            println("\tdischarge: $(result.discharge)")
            println("\tcharge: $(result.charge)")
            println("\tlevel: $(result.level)")
        end
    end
    # TODO: Split between nodes and lines
    if print_penalty
        println("Penalty Sum per Node:")
        node_to_result = admm.results[iteration].node_to_result
        for (node, node_result) in node_to_result
            println("\t$(node.name): $(node_result.penalty_term)")
        end
        println("Penalty Sum Total: $(admm.results[iteration].penalty_term)")
    end
end

function print_duals(iteration::Int)
    println(
        "\nlambdas: $(admm.lambdas[iteration])"
        * "\nmues: $(admm.mues[iteration])"
        * "\nrhos: $(admm.rhos[iteration])\n"
    )
end

begin
    admm = ADMM(0.0001, nodes, generators, storages, lines)

    for i in 1:2000
        calculate_iteration()
        
        println("Generation Results: ")
        print_results(true, true, false, admm.iteration)
        
        update_duals()
        
        check_convergence()
        if admm.converged
            break
        end
    end
end

function get_nodal_price(iteration::Int)
    nodal_price = zeros(length(admm.N), length(admm.T))
    for t in admm.T
        nodal_price[:, t] = admm.lambdas[iteration][t] .+ sum((admm.mues[iteration][l, t] + admm.rhos[iteration][l, t])  * admm.ptdf[l, :] for l in admm.L)
    end
    return nodal_price
end


plot_mues(2)
plot_rhos(1)
plot_lambdas()