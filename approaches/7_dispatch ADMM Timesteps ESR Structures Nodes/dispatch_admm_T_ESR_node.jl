using JuMP
using Gurobi
using PlotlyJS

include("data_node.jl")


# Trick to avoid multiple license printing
gurobi_env = Gurobi.Env()

node1 = Node("N1", [70, 120], false)
node2 = Node("N2", [30, 120], false)
node3 = Node("N3", [20, 200], true)

nodes = [node1, node2, node3]

line1 = Line("L1", node2, node1, 40, 1)
line2 = Line("L2", node3, node1, 100, 1)
line3 = Line("L3", node2, node3, 100, 2)

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

    sum_G, sum_S_d, sum_S_c = getSumOfIterationResults(admm.iteration - 1)

    previous_G = getIterationResults(generator, admm.iteration - 1)
    
    # Penalty term
    @expression(
        sub,
        penalty_term[t=admm.T], 
        (G[t] + sum_G[t] - previous_G[t] + sum_S_d[t]
            - sum_S_c[t] - admm.total_demand_t[t])^2
    )

    node_id = admm.node_to_id[generator.node]
    @objective(
        sub,
        Min,
        sum(G[t] * generator.marginal_costs
            + admm.lambdas_n_t[admm.iteration][node_id, t] * G[t]
            + admm.gamma/2 * penalty_term[t]
            for t in admm.T)
    )

    optimize!(sub)

    result = ResultGenerator(
        generator,
        value.(G).data,
        value.(penalty_term).data
    )

    return result
end

function optimize_subproblem(storage::Storage)
    sub = Model(with_optimizer(Gurobi.Optimizer, gurobi_env))
    set_silent(sub)

    # Maximum capacity of storage
    @variable(sub, 0 <= G_S_d[t=admm.T] <= storage.max_power)
    @variable(sub, 0 <= G_S_c[t=admm.T] <= storage.max_power)
    
    # Maximum storage level
    @variable(sub, 0 <= storage_level[t=admm.T] <= storage.max_level)

    sum_G, sum_S_d, sum_S_c = getSumOfIterationResults(admm.iteration - 1)

    previous_S_d, previous_S_c = getIterationResults(storage, admm.iteration-1)
    
    @expression(
        sub,
        penalty_term[t=admm.T], 
        (sum_G[t] + G_S_d[t] + (sum_S_d[t] - previous_S_d[t])
            - G_S_c[t] - (sum_S_c[t] - previous_S_c[t])
            - admm.total_demand_t[t])^2
    )

    node_id = admm.node_to_id[storage.node]
    @objective(
        sub,
        Min,
        sum(storage.marginal_costs * (G_S_d[t] + G_S_c[t])
            + admm.lambdas_n_t[admm.iteration][node_id, t] * (
                G_S_d[t] - G_S_c[t]
            ) + admm.gamma/2 * penalty_term[t]
            for t in admm.T)
    )

    @constraint(
        sub,
        StorageBalance[t=admm.T],
        storage_level[t] == (
            (t == 1 ? 0 : storage_level[t-1]) + G_S_c[t] - G_S_d[t]
        )
    )

    optimize!(sub)

    result = ResultStorage(
        storage,
        value.(G_S_d).data,
        value.(G_S_c).data,
        value.(storage_level).data,
        value.(penalty_term).data
    )

    return result
end

function update_lambda()
    lambdas = zeros(Float64, length(admm.N), length(admm.T))
    for node_id in 1:length(admm.N)
        sum_G_t, sum_S_d_t, sum_S_c_t = getSumOfIterationResults(admm.iteration)
        node_lambdas = (admm.lambdas_n_t[admm.iteration][node_id, :] 
            + admm.gamma * (
                sum_G_t + sum_S_d_t - sum_S_c_t - admm.total_demand_t
            )
        )
        lambdas[node_id, :] = node_lambdas
    end
    push!(admm.lambdas_n_t, lambdas)
end

function check_convergence()
    eps = 10^(-5)
    if admm.iteration != 1
        # Iteration counter is not increased yet. Hence, the new lamba is
        # compared to the lambda of the current iteration.
        admm.converged = all(
            abs.(admm.lambdas_n_t[end] - admm.lambdas_n_t[admm.iteration])
            .< eps
        )
    end
    if admm.converged
        println("Converged")
    else
        println("Not converged")
        admm.iteration += 1
    end
end

function calculate_iteration()
    println("\nIteration: $(admm.iteration) with "
        * "lambdas: $(admm.lambdas_n_t[admm.iteration])")

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
    if print_penalty
        println("Penalty Sum per Node:")
        node_to_sum_penalty_t = admm.results[iteration].node_to_sum_penalty_t
        for (node, sum_penalty_t) in node_to_sum_penalty_t
            println("\t$(node.name): $(sum_penalty_t)")
        end
        println("Penalty Sum Total: $(admm.results[iteration].sum_penalty_t)")
    end
end

function plot_lambdas(node::Node)
    lambdas = []
    for lambda_n_t in admm.lambdas_n_t
        push!(lambdas, lambda_n_t[admm.node_to_id[node], :]) 
    end
    node_lambda_matrix = hcat(lambdas...)
    traces = [
        scatter(x=1:admm.iteration, y=node_lambda_matrix[1, :], name="t=1")
    ]
    for t in admm.T
        if t != 1
            push!(
                traces,
                scatter(
                    x=1:admm.iteration,
                    y=node_lambda_matrix[t, :],
                    name="t=$(t)"
                )
            )
        end
    end
    plot(traces, Layout(title="Lambda for node $(node.name)"))
end

function plot_sum_penalty()
    penalty_terms = []
    if admm.converged
        iterations = 1:admm.iteration
    else
        iterations = 1:admm.iteration-1
    end
    for i in iterations
        push!(penalty_terms, admm.results[i].sum_penalty_t)
    end
    penalty_matrix = hcat(penalty_terms...)
    traces = [scatter(x=1:admm.iteration, y=penalty_matrix[1, :], name="t=1")]
    for t in admm.T
        if t != 1
            push!(
                traces,
                scatter(
                    x=1:admm.iteration,
                    y=penalty_matrix[t, :],
                    name="t=$(t)"
                )
            )
        end
    end
    plot(traces)
end

function plot_generation()
    generation = Array{Float64}(undef, 0, length(admm.T))
    charge = Array{Float64}(undef, 0, length(admm.T))
    labels = Array{String}(undef, 1,0)
    colors = Array{String}(undef, 1,0)
    colors_charge = Array{String}(undef, 1,0)
    for (unit, result) in admm.results[end].unit_to_result
        if typeof(unit) == Generator
            generation = vcat(generation, result.generation')
        elseif typeof(unit) == Storage
            generation = vcat(generation, result.discharge')
            charge = vcat(charge, -result.charge')
            colors_charge = hcat(colors_charge, unit.plot_color)
        end
        labels = hcat(labels, unit.name)
        colors = hcat(colors, unit.plot_color)
    end
    areaplot(
        generation',
        label=labels,
        color=colors,
        legend=:bottomright,
        xlabel="Timesteps",
        ylabel="MW",
        width=0
    )
    demand_t = zeros(length(admm.T))
    for node in admm.nodes
        demand_t += node.demand
    end
    plot!(demand_t, color=:black, width=1, label = "Demand", linestyle=:dash)
    areaplot!(
        charge',
        label="",
        color=colors_charge,
        width=0
    )
    hline!([0], color=:black, width=2, label="")
end

begin
    admm = ADMM(0.008, nodes, generators, storages, lines)

    for i in 1:200
        calculate_iteration()
        
        println("Generation Results: ")
        print_results(true, true, true, admm.iteration)
        
        update_lambda()
        print("Updated Lambdas: ")
        println(admm.lambdas_n_t[admm.iteration + 1])
        
        check_convergence()
        if admm.converged
            break
        end
    end
end

plot_lambdas(node1)
plot_sum_penalty()
using Plots
plot_generation()