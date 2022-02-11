mutable struct ADMM
    iteration::Int
    gamma::Float64
    lambdas::Vector{Vector{Float64}}
    mues::Vector{Matrix{Float64}}
    rhos::Vector{Matrix{Float64}}
    T::Vector{Int}
    N::Vector{Int}
    L::Vector{Int}
    nodes::Union{Vector{Node}, Nothing}
    generators::Vector{Generator}
    storages::Union{Vector{Storage}, Nothing}
    lines::Union{Vector{Line}, Nothing}
    results::Vector{Result}
    convergence::Convergence
    ptdf::Matrix{Float64}
    total_demand::Vector{Float64}
    node_id_to_demand::Dict{Int, Vector{Int}}
    node_to_id::Dict{Node, Int}
    node_to_units::Dict{Node, Vector{Union{Generator, Storage}}}
    L_max::Vector{Float64} 

    function ADMM(gamma::Float64,
                  nodes::Union{Vector{Node}, Nothing},
                  generators::Vector{Generator},
                  storages::Union{Vector{Storage}, Nothing},
                  lines::Union{Vector{Line}, Nothing})
        admm = new()
        admm.iteration = 1
        admm.gamma = gamma
        admm.T = collect(1:length(nodes[1].demand))
        admm.N = collect(1:length(nodes))
        admm.L = collect(1:length(lines))
        admm.lambdas = [zeros(Float64, length(admm.T))]
        admm.mues = [zeros(Float64, length(admm.L), length(admm.T))]
        admm.rhos = [zeros(Float64, length(admm.L), length(admm.T))]
        admm.nodes = nodes
        admm.generators = generators
        admm.storages = storages
        admm.lines = lines
        admm.L_max = [line.max_capacity for line in admm.lines]
        admm.results = []
        admm.convergence = Convergence()
        admm.ptdf = calculate_ptdf(admm.nodes, admm.lines)
        admm.total_demand = zeros(length(admm.T))
        admm.node_id_to_demand = Dict()
        admm.node_to_id = Dict()
        admm.node_to_units = Dict()
        if !isnothing(admm.nodes)
            for (id, node) in enumerate(admm.nodes)
                admm.total_demand += node.demand
                admm.node_id_to_demand[id] = node.demand
                admm.node_to_id[node] = id
            end
        end
        for unit in vcat(admm.generators, admm.storages)
            if haskey(admm.node_to_units, unit.node)
                push!(admm.node_to_units[unit.node], unit)
            else
                admm.node_to_units[unit.node] = [unit]
            end
        end
        return admm
    end
end