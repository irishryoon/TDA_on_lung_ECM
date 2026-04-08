# Topological analysis of lung adenocarcinoma

This repository contains code for paper <a href="https://www.biorxiv.org/content/10.1101/2024.01.05.574362">"Deciphering the diversity and sequence of extracellular matrix and cellular spatial patterns in lung adenocarcinoma using topological data analysis"</a>. 

# Install
Download <a href="https://julialang.org/downloads/">Julia</a>. We recommend downloading Julia v.1.10.2. For detailed instructions on installing Julia and the necessary packages, see [Julia_install_instructions.pdf][Julia_install_instructions.pdf\]

# Tutorial: Quick-start on computing topological features

## 1. Computing persistent homology features

* The code computes topological features from point cloud data. 
* The input file is a csv file with the xy- coordinates of the point clouds. See `compute_topological_features/test_input` for example input files. 
	* The point clouds can represent cell locations. They can also be points sampled from ECM. 
* One can compute topological features from multiple point clouds that are extracted from distinct images. In such cases, ensure the following:
	* The original images from which the point clouds are extracted must be of the same size.
	* Each point cloud must be saved in a separate CSV file. See `compute_topological_features/test_input` for example. 

### Using Julia REPL
* To compute the features, start Julia. Make sure you are in the root directory of this git repository. Run the following:

```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/compute_persistence.jl")
```

The script will prompt you to input the directory of the input CSV files and your desired output directory.

For example, one can use `tutorial/PH_inputs` and `tutorial/PH_outputs` when prompted to provide the input and output directories. 

### From command line
Alternatively, one can use the following script. 

```
julia src/compute_persistence_script.jl "PATH_TO_INPUT_DIRECTORY" "PATH_TO_OUTPUT_DIRECTORY"
```

For example, 

`julia src/compute_persistence.jl "tutorial/PH_inputs" "tutorial/PH_outputs"`


### Outputs
The outputs are the following:
* persistence diagrams (dimensions 0 and 1)
	* `PD` contains the persistence diagrams, saved as an array.
	* `PD_figures` contains plots of the persistence diagrams. 
* persistence images (dimensions 0 and 1)
	* `PI` contains the persistence images, saved as an array.
	* `PI_figures` contains plots of the persistence images. 

## 2. Computing Dowker persistent homology features
* The code computes Dowker persistece diagrams and Dowker persistence images from two point clouds.
* The input directory must have two subdirectories, one for each cell type. Let's say `cell_type1` and `cell_type2` are the names of the subdirectories. Here, cell type can refer to cancer cells, leukocytes, or ECM. 
	* In the subdirectories `cell_type1` and `cell_type2`, there must be csv files with the same name. This indicates that the csv files in the two subdirectories represent the point cloud from the same image. 
	* As before, each csv file is the xy-coordinates of the point clouds. 
* Note: Dowker persistence diagram computation takes longer than computing the usual persistence diagram. We therefore usually compute Dowker persistence diagrams on a subsample of the point cloud. Because of the random nature of subsampling, the subsampled cells will lead to slightly different Dowker persistence diagrams. Any analysis that uses Dowker persistence diagrams, such as PCA or UMAP, will also depend on the subsampling. 

### Using Julia REPL
To compute the Dowker features, start Julia. Make sure you are in the root directory of this git repository. Run the following:

```
using Pkg
Pkg.activate(".")
Pkg.instantiate()
include("src/compute_Dowker_persistence.jl")
```

The script will prompt you to input the directories of the two cell types and the output directory.

For example, one can use `tutorial/DowkerPH_inputs/cell_type1`, `tutorial/DowkerPH_inputs/cell_type2`, and `tutorial/DowkerPH_outputs` when prompted to provide the input and output directories. 

### Using command line
Alternatively, one can use the following script:
```
julia src/compute_Dowker_persistence_script.jl "PATH_TO_INPUT_DIRECTORY_1" "PATH_TO_INPUT_DIRECTORY_2" "PATH_TO_OUTPUT_DIRECTORY"
```

For example,
`julia src/compute_Dowker_persistence_script.jl "tutorial/DowkerPH_inputs/cell_type1" "tutorial/DowkerPH_inputs/cell_type2" "tutorial/DowkerPH_outputs"
`


# Analysis
* See the "Analysis" directory for topological feature computation and analysis on a small example dataset. 
* `src/ECM_TDA.jl` contains code for computing topological (PH and Dowker PH features), along with other helper functions that are relevant to the analysis
* To use the script, run the following.
```
include("src/ECM_TDA.jl")
using .ECM_TDA
```
