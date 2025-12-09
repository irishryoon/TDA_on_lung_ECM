"""
This script computes Dowker persistence in dimension 0 and 1. 

This should be run as a script. It should not imported as a library.

@author Iris Yoon
irishryoon@gmail.com
"""

module compute_Dowker_script
include("dowker_persistence.jl")
include("Eirene_var.jl")
include("TDA.jl")

using .Dowker
using .Eirene_var
using .TDA
using DelimitedFiles
using DataFrames
using Plots



function main()
    ##########################################################################
    # Set directories
    ##########################################################################

    # prompt to enter input directory
    print("\nEnter input directory 1: \n")  
    input_dir1 = readline() * "/"

    print("\nEnter input directory 2: \n")  
    input_dir2 = readline() * "/"

    # check if input directory exists
    if isdir(input_dir1) == false
        throw("The input directory does not exist")
    end 
    if isdir(input_dir2) == false
        throw("The input directory does not exist")
    end 

    # prompt to enter output directory
    print("\nEnter output directory: \n")
    output_dir = readline() * "/"

    # create output directory if it doesn't exist
    isdir(output_dir) || mkdir(output_dir)

    # output directory for persistence diagrams
    output_Dowker_PD = output_dir * "Dowker_PD/"
    output_Dowker_PD0 = output_dir * "Dowker_PD/Dowker_PD0/"
    output_Dowker_PD1 = output_dir * "Dowker_PD/Dowker_PD1/"

    isdir(output_Dowker_PD) || mkdir(output_Dowker_PD)
    isdir(output_Dowker_PD0) || mkdir(output_Dowker_PD0)
    isdir(output_Dowker_PD1) || mkdir(output_Dowker_PD1)

    ##########################################################################
    # compute Dowker persistence diagrams 
    ##########################################################################
    println("\nComputing Dowker persistence homology... \n")
    # use one of the input directories to get all files 
    csv_files = [item for item in walkdir(input_dir1)][1][3:end][1]
    DPD0_dict = Dict()
    DPD1_dict = Dict()

    for (idx, file) in enumerate(csv_files)
        if file[end-3:end] == ".csv"
    
            filename = split(file, ".")[1]
            # load the two point clouds 
            pointcloud1 = readdlm(input_dir1 * file, ',')
            pointcloud2 = readdlm(input_dir2 * file, ',')
    
            # compute Dowker PermutedDimsArray
            W_barcode0, W_barcode1, _ = compute_Dowker(pointcloud1, pointcloud2)
    
            # save Dowker persistence diagrams
            if W_barcode0 == nothing
                writedlm(output_Dowker_PD0 * file, zeros(), ",") 
                DPD0_dict[filename] = reshape(Array([0.0]), 1, 1)
            else
                writedlm(output_Dowker_PD0 * file, W_barcode0, ",") 
                DPD0_dict[filename] = W_barcode0 
            end
            
            if W_barcode1 == nothing
                writedlm(output_Dowker_PD1 * file, zeros(), ",")
                DPD1_dict[filename] = reshape(Array([0.0]), 1, 1)
            else
                writedlm(output_Dowker_PD1 * file, W_barcode1, ",")
                DPD1_dict[filename] = W_barcode1
            end
        end
    end

    
    ##########################################################################
    # compute Dowker persistence images
    ##########################################################################
    println("\nComputing Dowker persistence images... \n")
    # convert to Ripser PD
    PH0 = Dict(k => array_to_ripsererPD(v) for (k,v) in DPD0_dict if v != reshape(Array([0.0]), 1, 1))
    PH1 = Dict(k => array_to_ripsererPD(v) for (k,v) in DPD1_dict if v != reshape(Array([0.0]), 1, 1))

    # compute PI
    PI0 = compute_PI(PH0)
    PI1 = compute_PI(PH1);


    ### save PersistenceImages as arrays
    output_dir_PI = output_dir * "Dowker_PI/" 
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
    output_PD_figures = output_dir * "Dowker_PD_figures/"
    output_PD0_figures = output_PD_figures * "PD0/"
    output_PD1_figures = output_PD_figures * "PD1/"

    output_PI_figures = output_dir * "Dowker_PI_figures/"
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
    DPD0_dict = Dict(k =>v for (k,v) in DPD0_dict if v != nothing)
    DPD1_dict = Dict(k => v for (k,v) in DPD1_dict if v != nothing)
    # compute maximum PD values (for plotting)
    PD0_max = get_DPD0_max(DPD0_dict) # not a typo. This is because Dowker PD0 doesn't just end with one connected component
    PD1_max = get_PD1_max(DPD1_dict);


    ### save figures
    # dim 0
    for k in keys(DPD0_dict)
        p = TDA.plot_PD(DPD0_dict[k], pd_min = 0, pd_max = PD0_max * 1.01, frame = :box)
        savefig(output_PD0_figures * k * ".png")
    end

    # dim 1
    for k in keys(DPD1_dict)
        p = TDA.plot_PD(DPD1_dict[k], pd_min = 0, pd_max = PD1_max * 1.01, frame = :box)
        savefig(output_PD1_figures * k * ".png")
    end

    ### save persistence images
    # dim 0
    for k in keys(PI0)
        p = heatmap(PI0[k],  rightmargin=7Plots.mm, size = (500, 400))
        savefig(output_PI0_figures * k * ".png")
    end

    # dim 1
    for k in keys(PI1)
        p = heatmap(PI1[k],  rightmargin=7Plots.mm, size = (500, 400))
        savefig(output_PI1_figures * k * ".png")
    end

end

main()
end