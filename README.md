# GEN-tools

This project contains the scripts used in my master thesis:  
*“Integrated analysis of compositional multiomics time series data for the discovery of microbial
interactions and metabolic networks”*

The thesis was carried out in the working group of [Systems Biology of Microbial Communities](https://www.ufz.de/index.php?en=38437) at the Helmholtz Center of Environmental Research (UFZ), Leipzig. It was submitted on July 19th, 2019 by [Martin Graf](https://git.ufz.de/grafma). The [GEN-tools](https://git.ufz.de/UMBSysBio/magraf) project on GitLab is frozen and will further be maintained here on GitHub.

 ## Purpose of the project
Microbial Interference Networks (MIN) analysis has been shown to be a useful framework for the reconstruction of complex interactions between species. With the increasing popularity of metatranscriptomic sequencing, approaches established in MIN are now applied to metatranscriptomic data in Gene Expression Network (GEN) analysis. 

In this project, I applied four conceptually different approaches of MIN, namely [CoNet](https://doi.org/10.1371/journal.pcbi.1002606) and the [sparCC](https://doi.org/10.1371/journal.pcbi.1002687), [gLasso]( https://doi.org/10.1093/biostatistics/kxm045) and [MB](http://pages.stat.wisc.edu/~shao/stat992/lasso-cons1.pdf) method from the [SpiecEasi](https://github.com/zdk123/SpiecEasi) package, on two sets of metatranscriptomic time series data. The scripts cover all aspects from the data preparation, preprocessing and filtering, initial data export and threshold calculation for CoNet, up to the network estimation and results export for the visualization and network analysis in Cytoscape. 


## Scripts
* [TableMerge.R]( TableMerge.R) is the script that is used to generate the initial metatranscriptomic count table. It is designed to handle a *.biom* file (or multiple .txt/.csv table files). The input can hold **multiple records from the same samples**. It the ```Biome2Table()``` function offers the possibility to identify the records of each sample by a **regular expression** and join them as one column (sample). (This function became necessary because the initial MT sequence files were too large to perform the aligned on the full sample and thus, some files were split into multiple records).

* [PreProcessing.R](PreProcessing.R) contains the preprocessing pipeline and export of the preprocessed count table inputs for all GEN tools. Preprocessing includes optional **zero replacement**, **rank-bases filtering** without omitting the counts **(preserves the closure)**, and **table splitting** (by experiment). 

* [CoNet_thresholds.R](CoNet_thresholds.R) generates the threshold files for CoNet. The threshold limits are set, based on the number of pairwise associations to consider.

* [GEN-tools.R]( GEN-tools.R) covers the network estimation with sparCC, gLasso, and MB (all from the [SpiecEasi](https://github.com/zdk123/SpiecEasi) package). Note that the parameters, that are used for the algorithms are optimized for a specific dataset in terms of accuracy and computational expenses. A description of how to adjust the parameters for a certain dataset is given [here](https://github.com/zdk123/SpiecEasi/blob/master/README.md). For visualization and further analysis of the networks, the .txt output files of this script can be imported into [Cytoscape]( https://cytoscape.org) via ```File->Import->Networks from file…```.

* [tools.R](tools.R) contains all functions that were used in this project. The functions are designed for a general usage and are transferable (or adaptable) to other projects. 


* [CoNet.bat](CoNet.bat) a batch script to easily run all four computational steps of CoNet on a **Windows** node outside of Cytoscape. Note that CoNet requires at least *Java 7*. I found that *Java 10* did not work as well and used an old version of *Java 8* (```jdk1.8.0_162```) instead. 
