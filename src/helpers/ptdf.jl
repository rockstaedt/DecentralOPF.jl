function calculate_ptdf(nodes::Vector{Node}, lines::Vector{Line})

    N = collect(1:length(nodes))
    SLACK = nothing

    incidence = zeros(Int, (length(lines), length(nodes)))
    susceptances = []

    for (l, line) in enumerate(lines)
        push!(susceptances, line.susceptance)
        for (n, node) in enumerate(nodes)
            if isnothing(SLACK)
                if node.slack
                    SLACK = n
                end
            end
            if line.from == node
                incidence[l, n] = 1
            elseif line.to == node
                incidence[l, n] = -1
            else
                incidence[l, n] = 0
            end
        end
    end
    
    B = Diagonal([susceptances...])
    # Line suceptance Matrix
    Bl = B*incidence
    # Nodes suceptance Matrix
    Bn = (incidence'*B)*incidence
    

    B_inv = zeros(length(N), length(N))

    B_inv[setdiff(N, SLACK), setdiff(N, SLACK)] = (
        inv(Bn[setdiff(N, SLACK), setdiff(N, SLACK)])
    )
    PTDF = Bl*B_inv
    return PTDF
end