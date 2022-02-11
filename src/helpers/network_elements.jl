function update(node_result::ResultNode,
    unit_result::Union{ResultStorage, ResultGenerator})
    if typeof(unit_result) == ResultStorage
        node_result.discharge += unit_result.discharge
        node_result.charge += unit_result.charge
    elseif typeof(unit_result) == ResultGenerator
        node_result.generation += unit_result.generation
    end
    node_result.penalty_term = sum_up(
        node_result.penalty_term,
        unit_result.penalty_term
    )
    return node_result
end

function get_nodal_price(iteration::Int)
    nodal_price = zeros(length(admm.N), length(admm.T))
    for t in admm.T
        nodal_price[:, t] = (
            admm.lambdas[iteration][t] .+ sum((admm.mues[iteration][l, t]
            + admm.rhos[iteration][l, t])  * admm.ptdf[l, :] for l in admm.L)
        )
    end
    return nodal_price
end