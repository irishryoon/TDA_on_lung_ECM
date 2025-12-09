"""
Script that computes persistent homology features (persistence diagrams, persistence images)
"""
module compute_persistence


include("TDA.jl")
using .TDA
using DelimitedFiles
using Plots
using ArgParse


function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table! s begin
        "input"
            help = "a positional argument. Input directory."
            required = true

        "output"
            help = "a positional argument. Output directory."
            required = true
    end

    return parse_args(s)
end


function main()
    ##########################################################################
    # Set directories
    ##########################################################################
    parsed_args = parse_commandline()


    input_dir = parsed_args["input"]
    output_dir = parsed_args["output"]

    input_dir = input_dir * "/"
    output_dir = output_dir * "/"

    if isdir(input_dir) == false
        throw("The input directory does not exist")
    end 

    # create output directory if it doesn't exist
    isdir(output_dir) || mkdir(output_dir)
    
    # prompt to enter input directory
    #print("\nEnter input directory: \n")  
    # input_dir = readline() * "/"

    # # check if input directory exists
    # if isdir(input_dir) == false
    #     throw("The input directory does not exist")
    # end 

    # # prompt to enter output directory
    # print("\nEnter output directory: \n")
    # output_dir = readline() * "/"

    # # create output directory if it doesn't exist
    # isdir(output_dir) || mkdir(output_dir)

    ##########################################################################
    # compute persistence diagrams 
    ##########################################################################

    println("\nComputing persistence diagrams... \n")
    PD_output = output_dir * "PD/"
    PD0_dict, PD1_dict = compute_PH_directory(input_dir, PD_output)

    ##########################################################################
    # compute persistence images
    ##########################################################################
    println("\nComputing persistence images... \n")
    # convert array to Ripserer PD
    PH0 = Dict(k => array_to_ripsererPD(v) for (k,v) in PD0_dict if v != nothing)
    PH1 = Dict(k => array_to_ripsererPD(v) for (k,v) in PD1_dict if v != nothing)

    # compute PI
    PI0 = compute_PI(PH0);
    PI1 = compute_PI(PH1);

    # save PersistenceImages as arrays
    output_dir_PI = output_dir * "PI/" 
    output_dir_PI0 = output_dir_PI * "PI0/"
    output_dir_PI1 = output_dir_PI * "PI1/"

    # create output directory if it doesn't exist
    isdir(output_dir_PI) || mkdir(output_dir_PI)
    isdir(output_dir_PI0) || mkdir(output_dir_PI0)
    isdir(output_dir_PI1) || mkdir(output_dir_PI1)

    for (k,v) in PI0
        writedlm(output_dir_PI0 * k * ".csv", v, ",")
    end

    for (k,v) in PI1
        writedlm(output_dir_PI1 * k * ".csv", v, ",")
    end
    
    ##########################################################################
    # save plots of persistence diagrams and persistence images 
    ##########################################################################
    println("\nSaving plots of persistence diagrams and persistence images... \n")
    output_PD_figures = output_dir * "PD_figures/"
    output_PD0_figures = output_PD_figures * "PD0/"
    output_PD1_figures = output_PD_figures * "PD1/"

    output_PI_figures = output_dir * "PI_figures/"
    output_PI0_figures = output_PI_figures * "PI0/"
    output_PI1_figures = output_PI_figures * "PI1/"

    # create output directory if it doesn't exist
    isdir(output_PD_figures) || mkdir(output_PD_figures)
    isdir(output_PD0_figures) || mkdir(output_PD0_figures)
    isdir(output_PD1_figures) || mkdir(output_PD1_figures)

    isdir(output_PI_figures) || mkdir(output_PI_figures)
    isdir(output_PI0_figures) || mkdir(output_PI0_figures)
    isdir(output_PI1_figures) || mkdir(output_PI1_figures)

    # get maximum values (for plotting purposes)
    PD0_dict = Dict(k =>v for (k,v) in PD0_dict if v != nothing)
    PD1_dict = Dict(k => v for (k,v) in PD1_dict if v != nothing)
    max0 = get_PD0_max(PD0_dict)
    max1 = get_PD1_max(PD1_dict)

    ### save figures
    # dim 0
    for k in keys(PD0_dict)
        p = histogram(PD0_dict[k][:,2], xlims = (0, max0), label = "")
        savefig(output_PD0_figures * k * ".png")
    end

    # dim 1
    for k in keys(PD1_dict)
        p = plot_PD(PD1_dict[k], pd_min = 0, pd_max = max1, frame = :box)
        savefig(output_PD1_figures * k * ".png")
    end

    ### save persistence images
    # dim 0
    for k in keys(PI0)
        p = heatmap(PI0[k], rightmargin=7Plots.mm, size = (500, 400))
        savefig(output_PI0_figures * k * ".png")
    end

    # dim 1
    for k in keys(PI1)
        p = heatmap(PI1[k], rightmargin=7Plots.mm, size = (500, 400))
        savefig(output_PI1_figures * k * ".png")
    end

end

main()


end