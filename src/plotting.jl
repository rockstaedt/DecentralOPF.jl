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

function plot_mues(t)
    traces = [scatter(
        x=1:admm.iteration,
        y=[admm.mues[i][l, t] for i in 1:admm.iteration],
        name="l=$(l)")
        for l in admm.L]
    PlotlyJS.plot(traces)
end

function plot_rhos(t)
    traces = [scatter(
        x=1:admm.iteration,
        y=[admm.rhos[i][l, t] for i in 1:admm.iteration],
        name="l=$(l)")
        for l in admm.L]
    PlotlyJS.plot(traces)
end

function plot_generation(generator::Generator, t::Int)
    traces = [scatter(
        x=1:admm.iteration-1,
        y=[admm.results[i].unit_to_result[generator].generation[t] for i in 1:admm.iteration-1],
        name=generator.name
    )]
    PlotlyJS.plot(traces)
end

function plot_line_utilization(t::Int)
    traces = [
        scatter(
            x=1:admm.iteration-1,
            y=[admm.results[i].line_utilization[l, t] for i in 1:admm.iteration-1],
            name="Utilization for line $(l)"
        ) for l in admm.L
    ]
    PlotlyJS.plot(traces)
end

function plot_slack_variables(line_id::Int)
    traces = [
        scatter(
            x=1:admm.iteration-1,
            y=[admm.results[i].avg_R_ref[line_id] for i in 1:admm.iteration-1],
            name="R_ref for $(line_id)"
        ),
        scatter(
            x=1:admm.iteration-1,
            y=[admm.results[i].avg_R_cref[line_id] for i in 1:admm.iteration-1],
            name="R_cref for $(line_id)"
        ),
        ]
    PlotlyJS.plot(traces)
end

function plot_injection(t::Int)
    traces = [
        scatter(
            x=1:admm.iteration-1,
            y=[admm.results[i].injection[n] for i in 1:admm.iteration-1],
            name="I for node $(n)"
        ) for n in admm.N
    ]
    PlotlyJS.plot(traces)
end

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