# Code for creating point cloud from ECM image.
# This is separated from ECM_TDA.jl because using "Images" clashes with "Ripserer"
module preprocess

using Random
using Distributions
using Images
using Plots
using CSV
using DataFrames
export n_ECM_samples,
        sample_ECM_points,
       find_island_points,
       get_ECM_mean_pixels,
       sample_points_from_ECM_directory

#-----------------------------------------------------------------------------------
# functions for sampling points from ECM
#-----------------------------------------------------------------------------------
"""
    n_ECM_samples()
Determines the number of points to sample from mean ECM pixel value. Creates a piecewise-linear function with three pieces,
with f(low) = min_sample and f(high) = max_sample
The argument `x` is mean ECM pixel value, scaled to be [0,1], where x close to 1 means there is higher ECM content.
"""
function n_ECM_samples(x, low, high, min_sample, max_sample)
    if x <= low
        return min_sample
    elseif x <= high
        slope = (max_sample - min_sample) / (high - low)
        return slope * (x - low) + min_sample
    else
        return max_sample
    end
end

"""
    sample_ECM_points(image; <keyword arguments>)
Sample points from an ECM image. 
1. For each pixel (with pixel value p in [0,1]), sample points according to the binomial distribution with p^2 as probability.  
2. Remove "island points" 
3. Downsample
"""

"""
    sample_ECM_points(image; <keyword arguments>)
Sample points from an ECM image. 
1. For each pixel (with pixel value p in [0,1]), sample points according to the binomial distribution with p^2 as probability.  
2. Remove "island points" 
3. Downsample
"""
function sample_ECM_points(image; 
    invert = true,
    vicinity = 100,
    n_points = 5,
    n_samples = 2000)
    
    
    if invert == true
        image = 1 .- image
    end
    
    ### 1. Sample points according to binomial distribution
    # for each pixel location with pixel value p, sample the point according to Binomial(1, p^2)
    sampled_img = image.^2 .> rand(Uniform(0,1), size(image))
    
    # get index of sampled points
    inds = Tuple.(findall(!iszero, sampled_img))
    points = hcat(first.(inds), last.(inds))
    
    ### 2. Remove island points 
    image_size = size(image,1)
    island_idx = find_island_idx(points, sampled_img, image_size; vicinity = vicinity, n_points = n_points) 
    
    ### 3. Downsample
    include_idx = setdiff(1:size(points,1), island_idx)

    if length(include_idx) <= n_samples
        sampled_idx = include_idx
    else 
        sampled_idx = sample(include_idx, n_samples, replace = false)
    end
    sampled_points = points[sampled_idx,:]

    # for consistency
    points = points[:,[2,1]] # needed for consistency when plotted 
    sampled_points = sampled_points[:,[2,1]] # needed for consistency when plotted 
    
    return sampled_points, points, sampled_img, island_idx
end

"""
    find_island_idx(points; <keyword arguments>)
Given sampled points from an ECM image, identify points that are "islands". That is, identify points that have less than a specified number of points in its vicinity. 
"""
function find_island_idx(points, sampled, subregion_size; vicinity = 50, n_points = 5)
    island_idx = []

    for i = 1:size(points,1)
        # look around neighborhood of size vicinity
        vicinity_half = Int(vicinity/2)
        x, y = points[i,:]

        xmin = maximum([x-vicinity_half, 1])
        xmax = minimum([x+vicinity_half, subregion_size])

        ymin = maximum([y-vicinity_half, 1])
        ymax = minimum([y+vicinity_half, subregion_size])
        sub = sampled[xmin:xmax, ymin:ymax]
        if sum(sub) < n_points
            push!(island_idx, i)
        end
    end
    return island_idx
end

function get_ECM_mean_pixels(directory)
    files =  [item for item in walkdir(directory)][1][3]
    ECM_mean_pixels = []
    for filename in files
        if filename != ".DS_Store"
            img = Array(Images.load(directory * filename))
            push!(ECM_mean_pixels, mean(Float64.(img)))
        end
    end
    return ECM_mean_pixels
end

function sample_points_from_ECM_directory(ECM_directory, min_sample, max_sample, low, high; 
    c_ECM =  "#259ea1",
     sample_plot_directory = "",
     sample_CSV_directory = "")
    files = [item for item in walkdir(ECM_directory)][1][3]
    files = [f for f in files if f != ".DS_Store"]
    for filename in files
        img = Array(Images.load(ECM_directory * filename))
        figure_file = split(filename,".")[1] * ".pdf"

        # compute (inverted) mean pixel value of image
        img_mean_inv = 1- mean(Float64.(img))

        # compute number of points to sample
        n_sample = Int64(round(n_ECM_samples(img_mean_inv, low, high, min_sample, max_sample)))

        # sample points
        resampled, points, sampled, island_idx = sample_ECM_points(img, vicinity = 100, n_points = 5, n_samples = n_sample)

        # plot the results 
        p = scatter(resampled[:,1], resampled[:,2], yflip = :true, c = c_ECM, label = "", frame = :box, ticks = [], size = (500, 500))
        savefig(sample_plot_directory * figure_file)

        # save sampled points to CSV
        csv_file = split(filename, ".")[1] * ".csv"
        df = DataFrame(resampled, [:x, :y])
        CSV.write(sample_CSV_directory * csv_file, df)
    end
end

end