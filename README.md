# repDilPCR

<div align="justify">
repDilPCR is a software tool to analyze qPCR data. It has been inspired by the efficient dilution-replicate design for real-time PCR assays by Hui and Feng (Kwokyin Hui & Zhong-Ping Feng (2013) Efficient experimental design and analysis of real-time PCR assays, Channels, 7:3, 160-170, DOI: https://doi.org/10.4161/chan.24024) and is the first tool to enable the analysis of experiments performed according to this design. The statistical and the graphical functions of the program can also be used with preprocessed data obtained by more conventional assay designs and evaluation methods.

## Key features
* Ability to use multiple reference genes
* Imputation of missing Cq values (only for reference genes)
* Statistical hypothesis testing with guided selection of appropriate statistical tests
* Preparation of publication-ready plots
* High level of automatization of the whole analysis
* Open-source
* Fast and easy to use
  * single-click analysis of a whole dataset
  * ~1 min from uploading raw Cq values to getting publishable plots
* Multiple customization options
* Possibility to export tabular data at each intermediate step to analyze with third-party programs

## Introduction
In a qPCR experiment, it is of key importance to determine the efficiency of the PCR reaction for each amplicon and primer pair for correct evaluation and interpretation of the data. Different approaches to determine efficiency have been developed, from the classical calibration curve-based method to sophisticated methods that rely on fitting linear or non-linear models on individual amplification curves. Occupying the middle ground between these two extremes is the dilution-replicate experimental design of Hui and Feng that has remained somehow overlooked, most probably due to the lack up to now of a dedicated software tool to apply the method. This is a multiple linear regression-based approach with a number of advantages. It requires fewer reactions than the traditional approach with calibration curves produced by a separate set of dilutions of a standard sample. In the dilution-replicate design, standard curves are determined from so-called dilution-replicates of experimental reactions that serve both to control technical variance and to determine efficiency. Like this, all samples contribute to the efficiency estimate, thus precision increases with the number of samples on a plate. Furthermore, the traditional approach requires that the linear dynamic range of the independent standard curve covers all sample Cq values which sometimes leads to the necessity to repeat experiments using different dilutions. In contrast, with the dilution-replicate design it is guaranteed that the sample Cq values will be within range.
repDilPCR utilizes the described dilution-replicate analytical method and extends it by adding the possibility to use multiple reference genes. It also offers capabilities for performing statistical tests and plotting publication-ready graphs. The program has been designed with the philosophy to automate and speed up analysis of qPCR data (typically less than one minute from raw Cq values to publication-ready plots) and to help users with little knowledge of statistics to select and perform the appropriate statistical tests, at least in the case of one-factor experimental designs. At the same time, the program allows experienced users to export intermediate data and perform more sophisticated analyses with external statistical software, e.g. if two-way ANOVA is necessary.