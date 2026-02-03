"""
Library containing functions for topological data analysis
"""
module TDA
using Ripserer
using DataFrames
using CSV
using DelimitedFiles
using PersistenceDiagrams
using Plots

export  main,
        run_PH,
        compute_PH_directory,
        RipsererPD_to_array,
        array_to_ripsererPD,
        compute_PI,
        get_PD0_max,
        get_PD1_max,
        get_PD0_max2,
        get_DPD0_max,
        plot_PD



"""
    run_PH(df)
Computes persistence diagrams in dimension 0 and 1 from locations

### Inputs
- `cells`: dataframe with columns `x` and `y` for (x,y) coordinates

### Outputs
- `PD0`: persistence diagram in dimension 0
- `PD1`: persistence diagram in dimension 1
"""
function run_PH(df)

    # convert to Ripser input
    P = [tuple(df[i, :x], df[i, :y]) for i = 1:size(df,1)]

    # ripser (cohomology)
    PD = ripserer(P)
    PD0 = RipsererPD_to_array(PD[1])
    PD1 = RipsererPD_to_array(PD[2])

    return PD0, PD1
end

# compute PH from an input directory. Saves persistence diagrams (dim 0 and 1) as arrays. 
# returns a dictionary with all persistence diagrams
function compute_PH_directory(input_dir, output_dir)

    # output directories for PD0 and PD 1
    output_dir_PD0 = output_dir * "PD0/"
    output_dir_PD1 = output_dir * "PD1/"

    # create output directory if it doesn't exist
    isdir(output_dir) || mkdir(output_dir)
    isdir(output_dir_PD0) || mkdir(output_dir_PD0)
    isdir(output_dir_PD1) || mkdir(output_dir_PD1)

    # load all CSV files
    csv_files = [item for item in walkdir(input_dir)][1][3:end][1]
    csv_files = [f for f in csv_files if f != ".DS_Store"]

    PD0_dict = Dict()
    PD1_dict = Dict()

    for file in csv_files
        filename = split(file, ".")[1]

        df = CSV.read(input_dir * file, DataFrame)
        if size(df,1) != 0
            PD0, PD1 = run_PH(df)
        
            # save
            writedlm(output_dir_PD0 * filename * ".csv", PD0, ",")
            writedlm(output_dir_PD1 * filename * ".csv", PD1, ",")
            
            # combine to dictionary
            PD0_dict[filename] = PD0
            PD1_dict[filename] = PD1
        else
            PD0 = nothing
            PD1 = nothing
            
            writedlm(output_dir_PD0 * filename * ".csv", zeros(), ",")
            writedlm(output_dir_PD1 * filename * ".csv", zeros(), ",")
            
            # combine to dictionary
            PD0_dict_[filename] = PD0
            PD1_dict[filename] = PD1
        end
    end
    return PD0_dict, PD1_dict
end

# convert to 2D arrays
function RipsererPD_to_array(PD)
    n = size(PD,1)
    PD_array = zeros(n, 2)
    for i = 1:n
        PD_array[i,:] .= PD[i]
    end
    return PD_array
end

array_to_ripsererPD(PD_array) = PersistenceDiagram([(PD_array[i,1], PD_array[i,2]) for i = 1:size(PD_array,1)])

function compute_PI(PH_dict; sigma = 50, size = 20)
    
    PI = PersistenceImage([PH_dict[k] for k in keys(PH_dict)], sigma=sigma, size = size)
    PH_PI = Dict()
    for i in keys(PH_dict)
        PH_PI[i] = PI(PH_dict[i])
    end
    return PH_PI
end

# get maximum values of persistence diagrams (for plotting purposes)
get_PD0_max(PD_dict) = maximum([maximum(PD_dict[i][1:end-1,:]) for (i,v) in PD_dict if v != reshape(Array([0.0]), 1, 1) ])
get_PD1_max(PD_dict) = maximum([maximum(PD_dict[i]) for (i,v) in PD_dict if v != reshape(Array([0.0]), 1, 1) ])
get_PD0_max2(PD_dict) = maximum([sort(hcat(PD_dict[i]...), dims = 1)[end-1] for (i,v) in PD_dict if v != reshape(Array([0.0]), 1, 1) ])


get_DPD0_max(PD_dict) = maximum([sort(hcat(PD_dict[i]...), dims = 2)[end-1] for (i,v) in PD_dict if v != reshape(Array([0.0]), 1, 1) ])



function plot_PD(barcode; 
    highlight = [], highlight_color = :deeppink2, cutoff = nothing, inf_coord = nothing,
    pd_min = nothing, pd_max = nothing, threshold_lw = 2, diagonal_lw = 2, inf_markerstrokewidth = 5,
    kwargs...)
    points = barcode


    if size(barcode,1) == 0
        # plot diagonal line
        p = plot([0, 1], [0, 1], 
        labels ="", 
        linewidth = diagonal_lw,
        framestyle = :box,
        xlims = (0,1),
        ylims = (0,1),
        aspect_ratio = :equal,
        color = "grey"; 
        kwargs...)
        return p
    end
    # find index of points with death parameter == death
    idx = findall(x -> x == Inf, points[:,2])

    # plot points with death < Inf
    idx2 = [i for i in 1:size(points,1) if i ∉ idx]
    p = scatter(points[idx2,1], points[idx2,2]; kwargs..., color = "grey", labels = "", alpha = 0.5)

    # find max death value
    max_death = maximum(points[idx2, 2])

    # plot points with death parameter == Inf
    if inf_coord == nothing
        death_plot = ones(size(idx,1)) * max_death
    else
        death_plot = ones(size(idx, 1)) * inf_coord
    end
    scatter!(p, points[idx,1], death_plot, marker = :xcross, 
            markerstrokewidth = inf_markerstrokewidth,
            aspect_ratio = :equal, legend=:bottomright, labels="", color ="red"; kwargs...)

    # plot diagonal line
    if pd_max == nothing
        
        min_birth = minimum(points[:,1]) * 0.8
        max_death = max_death * 1.1
        plot!(p, [min_birth, max_death], [min_birth, max_death], 
            labels ="", 
            linewidth = diagonal_lw,
            framestyle = :box,
            xlims = (min_birth, max_death),
            ylims = (min_birth, max_death),
            color = "grey"; 
            kwargs...)
    else
        max_death = pd_max
        min_birth = pd_min
        plot!(p, [min_birth, max_death], [min_birth, max_death], 
            labels ="", 
            linewidth = diagonal_lw,
            framestyle = :box,
            xlims = (min_birth, max_death),
            ylims = (min_birth, max_death),
            color = "grey"; 
            
            kwargs...)
    end

    # if highlight is provided, color specific points with the given color
    if highlight != []
        scatter!(points[highlight,1], points[highlight,2]; kwargs..., color = highlight_color, labels = "")
    end

    # plot the cutoff (dotted line) if provided
    if cutoff != nothing
        f(x) = x + cutoff
        plot!(f, linestyle = :dash, c = "black", label = "", linewidth = threshold_lw)
    end

    return p
end

end