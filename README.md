# repDilPCR

<div align="justify">
repDilPCR is a software tool to analyze qPCR data. It has been inspired by the efficient dilution-replicate design for real-time PCR assays by Hui and Feng (Kwokyin Hui & Zhong-Ping Feng. Efficient experimental design and analysis of real-time PCR assays. Channels 2013, 7:3, 160-170, https://doi.org/10.4161/chan.24024) and is the first tool to enable the analysis of experiments performed according to this design. The statistical and the graphical functions of the program can also be used with precalcocessed data obtained by more conventional assay designs and evaluation methods.

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

repDilPCR utilizes the described dilution-replicate analytical method (Kwokyin Hui & Zhong-Ping Feng. Efficient experimental design and analysis of real-time PCR assays. Channels 2013, 7:3, 160-170, https://doi.org/10.4161/chan.24024) and extends it by adding the possibility to use multiple reference genes. It also offers capabilities for performing statistical tests and plotting publication-ready graphs. The program has been designed with the philosophy to automate and speed up analysis of qPCR data (typically less than one minute from raw Cq values to publication-ready plots) and to help users with little knowledge of statistics to select and perform the appropriate statistical tests, at least in the case of one-factor experimental designs. At the same time, the program allows experienced users to export intermediate data and perform more sophisticated analyses with external statistical software, e.g. if two-way ANOVA is necessary.

Although the primary goal of the program is to enable analysis of qPCR data via the dilution-replicate approach, the statistical and plotting functions can also be used with precalcocessed data (relative expression values) obtained by usual assay designs and evaluation methods.

## Installation
The program can be installed on a local computer or on a server (see below). Alternatively, users can freely access a working installation hosted on a server at the University Hospital in Ulm, Germany (http://not-yet-available). This service is anonymous, does not require registration and complies with common standards for protection of user data: raw data uploaded by the user are processed on the server and used to generate results that can be downloaded by the user; after the user closes the session by closing the browser window all uploaded data and processed results are automatically deleted from the server.

#### Prerequisites
* A working installation of R (version 3.6.0 or more recent) on a computer with a Linux or Windows operating system. (Theoretically MacOS should be possible, too, but I haven't had the chance to test whether it works.) The Rstudio integrated development environment is recommended for convenient use of the script but not required.
* The following R packages have to be installed: `car`, `gridExtra`, `tidyverse`, `mice`, `PMCMRplus`, `ggbeeswarm` and `ggsignif` (needed for both the ordinary R script and the Shiny app), as well as `shiny`, `shinycssloaders` and `shinyalert` (needed for the Shiny app only). It's possible that `PMCMRplus` will initially fail to install on a Debian or Ubuntu Linux system. The solution is to first install the GNU Multiple Precision Arithmetic Library (e.g. `gmp-6.2.1.tar.lz`) from https://gmplib.org/, as well as the package `libmpfr-dev` (`sudo apt install libmpfr-dev`).

#### Installation on a local computer
Download the zip archive of all files in the repository by clicking on "Code" and then on "Download ZIP" on the GitHub page of the repDilPCR project or by following this direct download link: https://github.com/deyanyosifov/repDilPCR/archive/refs/heads/main.zip.  (Note: the possibility for download is currently available to invited users only. It will be be made available to the general public after publishing the code. For now, if you are interested to get a copy of the program, write to me at deyan.yosifov@uniklinik-ulm.de.) Unzip the archive, this action will create  a new directory named `repDilPCR-main` in the current directory. You can rename the new directory to `repDilPCR` or whatever other name you choose and place it in a convenient place on your computer.  For the purpose of this manual, we will assume that your installation is located in the directory `repDilPCR` in your home folder on a Linux machine, i.e. `~/repDilPCR`. If your situation is different, just replace the `~/repDilPCR` part in the further instructions with the actual path to your installation.

#### Installation on a server
This option is only possible on a server running Linux. Apart from the prerequisites stated above, you will need to install the Shiny Server. (It can be downloaded from https://www.rstudio.com/products/shiny/download-server/, detailed installation instructions are available at https://docs.rstudio.com/shiny-server/#install-shiny.) Installation of repDilPCR on a server is similar to installing on a local computer but the `repDilPCR` directory will have to be placed in `/srv/shiny-server/`. The `shiny` user must have read and write access to `/srv/shiny-server/repDilPCR` and its contents.

## Usage
The program can be used both as an ordinary R script on a local computer or as a Shiny app (either on a local computer or on a server) accessed through a web browser.

### Preparation of the data
This preparatory step is the same no matter whether you intend to use the R script or the Shiny app. Preparing the data will take you some time, especially the first time when you do it, but don't worry, what might seem like wasted time is actually very well invested time because it will bring you considerable time savings later on by enabling automated analysis by repDilPCR.

The input data have to be arranged in a CSV file following a specific format which is different depending on whether one wants to feed into the program raw Cq values obtained from an experiment performed according to the dilution-replicate approach or already calculated relative expression values. Exemplary CSV files for these two use cases are provided in the files `Test_data.csv` and `Test_data_precalc.csv`, respectively. Note for international users: a point must be used as the decimal separator and a comma as the field separator (to separate values in each row). Cq values can usually easily be exported by the software controlling your qPCR machine. An important thing that you have to pay attention to is that you select a common threshold for all amplicons in your experiment before exporting the Cq values. Depending on the manufacturer of the machine and the respective software, they might be referred to as Ct ("cycle threshold") or Cp ("crossing point") values but these different names stand for the same thing. Here, we adhere to the MIQE guidelines and the respective terminology (Cq = quantification cycle). If you have not prepared a CSV file before, don't panic. This is easily done in any commonly available spreadsheet software like Microsoft Excel or LibreOffice Calc. You just have to select the CSV file format when saving and to pay attention to the regional settings (point as the decimal separator and comma as the field separator).

**In the case of using raw Cq values,** the CSV file needs to have the following layout. The first row contains column titles. The first three columns have predetermined names that must not be changed. The first column is called "Replicates" and it should contain the names of the samples with a suffix that identifies the biological replicate. The suffix consists of an underscore ("_") plus additional numbers and/or letters. For example, if you have an experimental group called "HeLa" and your experiment contained biological triplicates, they could be entered in the column as "HeLa_1", "HeLa_2" and "HeLa_3". The same convention can be used if your experiment compares groups of patients and each patient provided a single sample. In this case, each subject would be considered a biological replicate of the respective group and assigned an individual suffix number, e.g. if you compared the expression of a gene in 3 groups of 10 subjects each - healthy, leukemia and lymphoma, your labels in the "Replicates" column might be "Healthy_1", "Healthy_2", ..., "Healthy_10", "Leukemia_1", "Leukemia_2", ..., "Lymphoma_1", "Lymphoma_2" and so on (or you could use subject name initials, instead of numbers, as long as they don't occur more than once). If your experiment did not include replicates, the single sample can be entered as either "HeLa_1" or just "HeLa" in the table, it doesn't matter. The underscore character is reserved for the purpose of identifying replicates and should not occur elsewhere in the names of the samples. The second column is called "Pairs" and can contain optional information about grouping of samples in pairs. This is useful if the experiment consists of samples or sample groups that should be compared and statistically evaluated in pairs, e.g. if you have a number of cell lines in each of which you have measured the expression of your gene of interest in two conditions (control and treated) and you are interested in comparing these two conditions for each cell line but don't care about comparisons among cell lines. Then, all replicates from both samples of a given cell line should be assigned the same number in the column "Pairs". For example, all replicates (biological as well as dilution replicates) of "HeLa control" and "HeLa treated" should receive a number 1 in  the "Pairs" column, all replicates of "Jurkat control" and "Jurkat treated" should receive a number 2 and so on. In this way, the program will know which samples constitute a pair. If your experiment does not contain pairs or you do not wish to use this functionality, leave the column empty but do not delete its title. The third column is called "Dilution" and contains the dilution factors according to the dilution-replicate design. For example, if you had 5-fold serial dilutions in your experiment, you can use as factors the numbers 1, 5 and 25 (the actual numbers do not matter, only the ratios among them are important, so you could just as well use 2, 10 and 50). If we return to the example given above, this means that for each biological replicate "HeLa_1", "HeLa_2" and "HeLa_3" you would have three dilution replicates, accordingly each biological replicate should be listed three times in the "Replicates" column and the respective rows in the "Dilution" column should be filled with the numbers 1, 5 and 25 where 1 corresponds to the most concentrated dilution replicate and 25 corresponds to the most diluted one. The following columns should contain the Cq values for the assessed genes, first the reference genes (RG) and then the genes of interest (GOI). The titles of these columns should be the names of the respective genes/amplicons. Here is a quick look at a part of the `Test_data.csv` file as an example of the required format:

```
Replicates,Pairs,Dilution,RG1,RG2,RG3,GOI1,GOI2
CellLine1_1,1,1,22.91,24.69,21.85,23.06,22.31
CellLine1_1,1,5,25.67,27.97,23.99,25.68,24.25
CellLine1_1,1,25,28.08,30.71,26.66,28.42,27.28
CellLine1_2,1,1,23.70,25.33,21.70,23.66,21.55
CellLine1_2,1,5,26.28,27.96,23.63,26.13,24.92
CellLine1_2,1,25,27.74,30.05,25.90,28.02,26.48
CellLine1_3,1,1,23.60,24.67,21.82,23.21,21.55
CellLine1_3,1,5,26.14,28.69,24.29,25.82,24.17
CellLine1_3,1,25,27.42,30.23,25.94,27.64,26.70
..............................................
CellLine1_9,1,1,22.86,24.96,21.44,22.04,21.34
CellLine1_9,1,5,25.45,27.97,23.81,25.49,25.08
CellLine1_9,1,25,27.48,30.03,26.40,28.16,26.29
CellLine1 (trt)_1,1,1,23.27,24.66,21.50,25.15,23.48
CellLine1 (trt)_1,1,5,25.60,27.38,24.32,26.98,25.74
CellLine1 (trt)_1,1,25,28.09,30.45,26.09,28.95,28.12
....................................................
CellLine1 (trt)_9,1,1,22.35,25.30,21.46,24.59,22.49
CellLine1 (trt)_9,1,5,25.58,27.69,24.27,26.67,25.71
CellLine1 (trt)_9,1,25,27.57,29.78,26.16,29.12,27.54
CellLine2_1,2,1,23.41,22.73,22.28,23.43,21.03
CellLine2_1,2,5,25.72,26.29,24.86,26.31,24.36
CellLine2_1,2,25,28.06,27.88,26.61,27.54,26.05
CellLine2_2,2,1,22.99,23.16,21.27,23.16,21.49
.............................................
```

In this example fragment, three experimental groups are present ("CellLine1", "CellLine1 (trt)" and "CellLine2"), the first two of which form a pair (baseline vs. treated), as the same number (1) is present in the second column of all rows occupied by replicates of these two samples ("CellLine2" is part of another pair, as its replicates have the number 2 in the second column.) Each experimental group has 9 biological replicates and each of them has 3 dilution replicates (denoted by the dilution factors 1, 5 and 25 in the third column). Note that we don't call them "technical replicates" as the concentration of the template in them is not identical, nevertheless they simultaneously fulfil two roles: producing calibration curves and controlling technical variance. The example fragment above includes three reference genes (labelled "RG1", "RG2" and "RG3") and two genes of interest ("GOI1" and "GOI2"). The order in which you enter your samples in the table is important because this will also be the order in which your samples will be displayed on the graphs, i.e. the sample that the program finds in the first row will become the leftmost sample on the graphs and so on.

**In the case of using relative expression values,** the CSV file that has to be prepared has a simpler layout. Again, the first row contains column titles but now only the first two columns are obligatory and with predetermined names that must not be changed: "Replicates" and "Pairs". Their specification is the same as in the case when raw Cq values are used (see above). The next columns should contain the relative expression levels (linearly scaled) of the evaluated genes of interest in each biological replicate. Accordingly, the titles of these columns would be the respective gene/amplicon names. Here is a quick look at a part of the `Test_data_precalc.csv` file as an example of the required format:

```
Replicates,Pairs,GOI1,GOI2
CellLine1_1,1,1.58,1.42
CellLine1_2,1,1.33,1.70
CellLine1_3,1,1.81,1.98
.......................
CellLine1_9,1,1.91,1.56
CellLine1 (trt)_1,1,0.58,0.60
.............................
CellLine1 (trt)_9,1,0.63,0.78
CellLine2_1,2,1.10,1.72
CellLine2_2,2,1.19,1.27
.......................
```

In this example fragment, three experimental groups are present ("CellLine1", "CellLine1 (trt)" and "CellLine2"), the first two of which form a pair (baseline vs. treated), as the same number (1) is present in the second column of all rows occupied by replicates of these two samples ("CellLine2" is part of another pair, as its replicates have the number 2 in the second column.) Each experimental group has 9 biological replicates. The example fragment above includes two genes of interest ("GOI1" and "GOI2"). The order in which you enter your samples in the table is important because this will also be the order in which your samples will be displayed on the graphs, i.e. the sample that the program finds in the first row will become the leftmost sample on the graphs and so on.

### Usage of the Shiny app
The Shiny app can be used via any modern web browser. Users have two options:
* access a publicly available Shiny server with repDilPCR installed on it, for example the installation hosted at the University Hospital of Ulm (http://not-yet-available)
* if repDilPCR is installed on their local computer (see above for how-to), they can start the R environment and issue the following commands:
  ```
  library(shiny)
  runApp("~/repDilPCR/app.R", launch.browser = TRUE)
  ```
  replacing the `~/repDilPCR` part with the actual path to their installation, if deviating. This should launch the program and automatically start a new browser    window or tab to access it.


### Usage of the R script


