# repDilPCR

<div align="justify">
repDilPCR is a software tool to analyze qPCR data. It has been inspired by the efficient dilution-replicate design for real-time PCR assays by Hui and Feng (Kwokyin Hui & Zhong-Ping Feng. Efficient experimental design and analysis of real-time PCR assays. Channels 2013, 7:3, 160-170, https://doi.org/10.4161/chan.24024) and is the first tool to enable the analysis of experiments performed according to this design. The statistical and the graphical functions of the program can also be used with preprocessed data obtained by more conventional assay designs and evaluation methods.

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

Although the primary goal of the program is to enable analysis of qPCR data via the dilution-replicate approach, the statistical and plotting functions can also be used with preprocessed data (relative expression values) obtained by usual assay designs and evaluation methods.

## Installation
The program can be installed on a local computer or on a server (see below). Alternatively, users can freely access a working installation hosted on a server at the University Hospital in Ulm, Germany (http://not-yet-available). This service is anonymous, does not require registration and complies with common standards for protection of user data: raw data uploaded by the user are processed on the server and used to generate results that can be downloaded by the user; after the user closes the session by closing the browser window all uploaded data and processed results are automatically deleted from the server.

#### Prerequisites
* A working installation of R (version 3.6.0 or more recent) on a computer with a Linux or Windows operating system. (Theoretically MacOS should be possible, too, but I haven't had the chance to test whether it works.) The Rstudio integrated development environment is recommended for convenient use of the script but not required.
* The following R packages have to be installed: `car`, `gridExtra`, `tidyverse`, `mice`, `PMCMRplus`, `ggbeeswarm` and `ggsignif` (needed for both the ordinary R script and the Shiny app), as well as `shiny`, `shinycssloaders` and `shinyalert` (needed for the Shiny app only). It's possible that `PMCMRplus` will initially fail to install on a Debian or Ubuntu Linux system. The solution is to first install the GNU Multiple Precision Arithmetic Library (e.g. `gmp-6.2.1.tar.lz`) from https://gmplib.org/, as well as the package `libmpfr-dev` (`sudo apt install libmpfr-dev`).

#### Installation on a local computer
Download the zip archive of all files in the repository by clicking on "Code" and then on "Download ZIP" on the GitHub page of the repDilPCR project or by following this direct download link: https://github.com/deyanyosifov/repDilPCR/archive/refs/heads/main.zip.  (Note: the possibility for download is currently available to invited users only. It will be be made available to the general public after publishing the code. For now, if you are interested to get a copy of the program, write to me at deyan.yosifov@uniklinik-ulm.de.) Unzip the archive, this action will create  a new directory named `repDilPCR-main` in the current directory. You can rename the new directory to `repDilPCR` or whatever other name you choose and place it in a convenient place on your computer.  For the purpose of this manual, we will assume that your installation is located in the directory `repDilPCR` in your home folder on a Linux machine, i.e. `~/repDilPCR`. If your situation is different, just replace the `~/repDilPCR` part in the further instructions with the actual path to your installation.

#### Installation on a server
This option is only possible on a server running Linux. Apart from the prerequisites stated above, you will need to install the Shiny Server. (It can be downloaded from a("here", href = "https://www.rstudio.com/products/shiny/download-server/"), detailed installation instructions are available a("here", href = "https://docs.rstudio.com/shiny-server/#install-shiny").) Installation of repDilPCR on a server is similar to installing on a local computer but the `repDilPCR` directory will have to be placed in `/srv/shiny-server/`. The `shiny` user must have read and write access to `/srv/shiny-server/repDilPCR` and its contents.

## Usage
The program can be used both as an ordinary R script on a local computer or as a Shiny app (either on a local computer or on a server) accessed through a web browser.

### Usage of the R script
#### Preparation of the data
The input data have to be arranged in a CSV file following a specific format. It is strongly recommended to perform each analysis in a dedicated directory. The recommended 

the  called `repDilPCR` somewhere on your computer. Download the files `repDilPCR_lib.R` and `repDilPCR.R` and place them into the `repDilPCR` directory. Optionally, you can download the `Test_data.csv` and `Test_data_prepr.csv` files, too. These files contain example data and can be used to test the functions of the program and as templates for the proper formatting of your own data. The first file (`Test_data.csv`) contains raw Cq values obtained from an experiment performed according to the dilution-replicate approach and can accordingly be used as template for the formatting of data obtained from such experiments. The second file (`Test_data_prepr.csv`) contains preprocessed data (relative expression values, linear scale) and can accordingly be used as template for the formatting of such data independent of the experimental design (dilution-replicate or classical, with or without calibration curves).
