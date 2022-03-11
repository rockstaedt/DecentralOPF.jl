include("imports.jl")

include("cases/three_node.jl")

admm = ADMM(0.3, nodes, generators, storages, lines)

run!(admm)

np = get_nodal_price(admm.iteration)