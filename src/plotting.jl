using PlotlyJS

function plot_lambdas()
    lambda_matrix = hcat(admm.lambdas...)
    traces = [scatter(
        x=1:admm.iteration,
        y=lambda_matrix[t, :],
        name="t=$(t)")
        for t in admm.T]
    PlotlyJS.plot(traces)
end

# TODO: Adapt function
# -> auf Leitung anpassen

# function plot_mues(node::Node)
#     mues = []
#     for mue_n_t in admm.mues
#         push!(mues, mue_n_t[admm.node_to_id[node], :]) 
#     end
#     node_mue_matrix = hcat(mues...)
#     traces = [
#         scatter(x=1:admm.iteration, y=node_mue_matrix[1, :], name="t=1")
#     ]
#     for t in admm.T
#         if t != 1
#             push!(
#                 traces,
#                 scatter(
#                     x=1:admm.iteration,
#                     y=node_mue_matrix[t, :],
#                     name="t=$(t)"
#                 )
#             )
#         end
#     end
#     plot(traces, Layout(title="Mue for line $(node.name)"))
# end

# TODO: Adapt function

# function plot_sum_penalty()
#     penalty_terms = []
#     if admm.converged
#         iterations = 1:admm.iteration
#     else
#         iterations = 1:admm.iteration-1
#     end
#     for i in iterations
#         push!(penalty_terms, admm.results[i].sum_penalty_t)
#     end
#     penalty_matrix = hcat(penalty_terms...)
#     traces = [scatter(x=1:admm.iteration, y=penalty_matrix[1, :], name="t=1")]
#     for t in admm.T
#         if t != 1
#             push!(
#                 traces,
#                 scatter(
#                     x=1:admm.iteration,
#                     y=penalty_matrix[t, :],
#                     name="t=$(t)"
#                 )
#             )
#         end
#     end
#     plot(traces)
# end

# TODO: Adapt function

# function plot_generation()
#     generation = Array{Float64}(undef, 0, length(admm.T))
#     charge = Array{Float64}(undef, 0, length(admm.T))
#     labels = Array{String}(undef, 1,0)
#     colors = Array{String}(undef, 1,0)
#     colors_charge = Array{String}(undef, 1,0)
#     for (unit, result) in admm.results[end].unit_to_result
#         if typeof(unit) == Generator
#             generation = vcat(generation, result.generation')
#         elseif typeof(unit) == Storage
#             generation = vcat(generation, result.discharge')
#             charge = vcat(charge, -result.charge')
#             colors_charge = hcat(colors_charge, unit.plot_color)
#         end
#         labels = hcat(labels, unit.name)
#         colors = hcat(colors, unit.plot_color)
#     end
#     areaplot(
#         generation',
#         label=labels,
#         color=colors,
#         legend=:bottomright,
#         xlabel="Timesteps",
#         ylabel="MW",
#         width=0
#     )
#     demand_t = zeros(length(admm.T))
#     for node in admm.nodes
#         demand_t += node.demand
#     end
#     plot!(demand_t, color=:black, width=1, label = "Demand", linestyle=:dash)
#     areaplot!(
#         charge',
#         label="",
#         color=colors_charge,
#         width=0
#     )
#     hline!([0], color=:black, width=2, label="")
# end