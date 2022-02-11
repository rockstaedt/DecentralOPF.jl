mutable struct Node
    name::String
    demand::Vector{Int}
    slack::Bool
end

mutable struct Generator
    name::String
    marginal_costs::Int
    max_generation::Int
    plot_color::String
    node::Node
end

mutable struct Storage
    name::String
    marginal_costs::Int
    max_power::Int
    max_level::Int
    plot_color::String
    node::Node
end

mutable struct Line
    name::String
    from::Node
    to::Node
    max_capacity::Int
    susceptance::Int
end