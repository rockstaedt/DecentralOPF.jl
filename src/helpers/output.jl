function export_results(admm::ADMM, filename::String)
    parent_dir = "results/"

    df_duals = get_dual_df(admm)
    CSV.write(parent_dir * filename * "_duals.csv", df_duals)

    df_generators = get_generator_df(admm)
    CSV.write(parent_dir * filename * "_generators.csv", df_generators)

    df_storages = get_storage_df(admm)
    CSV.write(parent_dir * filename * "_storages.csv", df_storages)
end

function get_dual_df(admm::ADMM)
    iterations = []
    duals = []
    timesteps = []
    lines = []
    values = []
    for (dual_name, dual_values) in [["lambda", admm.lambdas], ["rho", admm.rhos], ["mue", admm.mues]]
        for i in 1:admm.iteration
            for t in admm.T
                if dual_name == "lambda"
                    push!(iterations, i)
                    push!(duals, dual_name)
                    push!(timesteps, t)
                    push!(lines, "")
                    push!(values, dual_values[i][t])
                else
                    for l in admm.L
                        push!(iterations, i)
                        push!(duals, dual_name)
                        push!(timesteps, t)
                        push!(lines, l)
                        push!(values, dual_values[i][l,t])
                    end
                end
            end
        end
    end

    return DataFrame(iteration=iterations, dual=duals, timestep=timesteps, line=lines, value=values)
end

function get_generator_df(admm::ADMM)
    iterations = []
    generators = []
    timesteps = []
    generations = []
    for generator in admm.generators
        for i in 1:admm.iteration
            generator_result = admm.results[i].unit_to_result[generator]
            for t in admm.T
                push!(iterations, i)
                push!(generators, generator.name)
                push!(timesteps, t)
                push!(generations, generator_result.generation[t])
            end
        end
    end

    return DataFrame(iteration=iterations, generator=generators, timestep=timesteps, generation=generations)
end

function get_storage_df(admm::ADMM)
    iterations = []
    storages = []
    timesteps = []
    charges = []
    discharges = []
    for storage in admm.storages
        for i in 1:admm.iteration
            storage_result = admm.results[i].unit_to_result[storage]
            for t in admm.T
                push!(iterations, i)
                push!(storages, storage.name)
                push!(timesteps, t)
                push!(charges, storage_result.charge[t])
                push!(discharges, storage_result.discharge[t])
            end
        end
    end

    return DataFrame(iteration=iterations, storage=storages, timestep=timesteps, charge=charges, discharge=discharges)
end