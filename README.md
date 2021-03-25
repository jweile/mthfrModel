# Shifting landscapes of human MTHFR missense variant effects
_**A collection of scripts accompanying the manuscript**_

## Dependencies
* R version 3.6.3 or higher
* [OpenPyMol](https://github.com/schrodinger/pymol-open-source/releases) 1.8 or higher

Many of the scripts in this repository depend on the R packages below. Please install them before using.
* yogitools : Install via devtools from [jweile/yogitools](https://github.com/jweile/yogitools)
* hgvsParseR : Install via devtoolsfrom [VariantEffect/yogitools](https://github.com/VariantEffect/yogitools)
* yogiroc : Install via devtoolsfrom [jweile/yogiroc](https://github.com/jweile/yogiroc)
* maveLLR : Install via devtoolsfrom [jweile/maveLLR](https://github.com/jweile/maveLLR)
* pbmcapply : Install via CRAN
* beeswarm : Install via CRAN
* optimization : Install via CRAN
* mclust : Install via CRAN
* ncdf4 : Install via CRAN

## File descriptions
* **map_data/** : A folder containing the experimentally measured fitness values in the two genetic backgrounds and four folinate concentrations. W12, W25, W100 and W200 contain the WT background maps for 12.5, 25, 100 and 200ug/ml folinate, while V12, V25, V100 and V200 contain the Ala222Val background maps for the equivalent folinate concentrations.
* **reference_data/**: A folder containing reference datasets, encompassing MTHFR structural features (derived from the structural alignment generated below), curated genotype/phenotype information (in the files named `age_of_onset`); a copy of the MTHFR gnomAD database entry; results from a previous low-throughput study of MTHFR variants by Goyette and Rozen 2000; and in-silico variant effect predictions from PolyPhen-2, SIFT, and PROVEAN.
* **results/**: This folder is where the scripts will write their results to upon execution.
* **fitModels.R** : This is the foundational script which will run the linear model fits. It also generates lists of confident folinate-responsive residues and genetic interactions with the A222V background. It requires at least 8 CPU cores and takes roughly 2h to run on a contemporary machine.
* **MTHFR_structural_alignment.pml**: An OpenPyMol script that performs a structural alignment between Human and E.coli MTHFR stuctures to determine putative binding site information for folate, FAD and SAH.
* **compareGoyette.R**: This script compares the variant effect map data against previously published low-throughput data from Goyette and Rozen 2000.
* **findSuppressors.R**: This script uses the output of `fitModels.R` to determine variants which significatly suppress the phenotype of A222V.
* **hypercomp.R**: This script runs statistical tests to establish significant hypercomplementation in MTHFR's regulatory domain, and near the SAM binding site in particular.
* **prsNrs.R**: This script evaluates the performance of the map with respect to distinguishing pathogenic variants from gnomAD and 1000 genomes variants, based on the reference data sets described above.
* **molecular_dynamics.R**: This script analyzes molecular dynamics simulation data of the region surrounding Trp165. It analyzes the distance between Trp165 and FAD, as well as FAD and the active site center. It also performs a clustering analysis of the Trp165-FAD binding mode and generates a Markov model of the interaction. The underlying data itself surpasses github storage limits and is not provided here, but is available upon request. 

