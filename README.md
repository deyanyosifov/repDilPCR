# repDilPCR

<div align="justify">
repDilPCR is a software tool to analyze qPCR data. It has been inspired by the efficient dilution-replicate design for real-time PCR assays by Hui and Feng (Kwokyin Hui & Zhong-Ping Feng. Efficient experimental design and analysis of real-time PCR assays. Channels 2013, 7:3, 160-170, https://doi.org/10.4161/chan.24024) and is the first tool to enable the analysis of experiments performed according to this design. The statistical and the graphical functions of the program can also be used with preprocessed data obtained by more conventional assay designs and evaluation methods.

## Key features
* Fast and easy to use
  * single-click analysis of a whole dataset
  * ~1 min from uploading raw Cq values to getting publishable plots
* Ability to use multiple reference genes
* Imputation of missing Cq values (only for reference genes)
* Statistical hypothesis testing with guided selection of appropriate statistical tests
* Preparation of publication-ready plots
* High level of automatization of the whole analysis
* Open-source
* Multiple customization options
* Possibility to export tabular data at each intermediate step to analyze with third-party programs

## Installation
The program can be installed on a local computer or on a server (see below). Alternatively, users can freely access and use a working installation hosted on a server at the German Cancer Research Center in Heidelberg, Germany (https://repdilpcr.eu). This service is anonymous, does not require registration and complies with common standards for protection of user data: raw data uploaded by the user are processed on the server and used to generate results that can be downloaded by the user; after the user closes the session by closing the browser window all uploaded data and processed results are automatically deleted from the server (**Warning: if you use a local installation of repDilPCR, do not store your data in the folder where repDilPCR is installed!**).

#### Prerequisites for installation
* A working installation of R (version 3.6.0 or more recent) on a computer with a Linux or Windows operating system. (Theoretically MacOS should be possible, too, but I haven't had the chance to test whether it works.) The RStudio integrated development environment is recommended for convenient use of the script but not required.
* The following R packages have to be installed: `car`, `gridExtra`, `tidyverse`, `mice`, `PMCMRplus`, `scales`, `RColorBrewer`, `ggbeeswarm` and `ggsignif` (needed for both the ordinary R script and the Shiny app), as well as `shiny`, `shinycssloaders` and `shinyalert` (needed for the Shiny app only). It's possible that `PMCMRplus` will initially fail to install on a Linux system. The solution is to first install the GNU Multiple Precision Arithmetic Library (e.g. `gmp-6.2.1.tar.lz`) from https://gmplib.org/, as well as the GNU Multiple Precision Floating-Point Reliable Library (`sudo apt install libmpfr-dev` on a Debian-based distribution or `sudo yum install mpfr-devel` on a RedHat-based distribution).

#### Installation on a local computer
Download the zip archive of all files in the repository by clicking on "Code" and then on "Download ZIP" on the GitHub page of the repDilPCR project or by following this direct download link: https://github.com/deyanyosifov/repDilPCR/archive/refs/heads/main.zip.  Unzip the archive, this action will create  a new directory named `repDilPCR-main` in the current directory. You can rename the new directory to `repDilPCR` or whatever other name you choose and place it in a convenient place on your computer.  For the purpose of this manual, we will assume that your installation is located in the directory `repDilPCR` in your home folder on a Linux machine, i.e. `~/repDilPCR`. If your situation is different, just replace the `~/repDilPCR` part in the further instructions with the actual path to your installation.

#### Installation on a server
This option is only possible on a server running Linux. Apart from the prerequisites stated above, you will need to install the Shiny Server. (It can be downloaded from https://www.rstudio.com/products/shiny/download-server/, detailed installation instructions are available at https://docs.rstudio.com/shiny-server/#install-shiny.) Installation of repDilPCR on a server is similar to installing on a local computer but the `repDilPCR` directory will have to be placed in `/srv/shiny-server/`. The `shiny` user must have read and write access to `/srv/shiny-server/repDilPCR` and its contents.

#### Browser compatibility
repDilPCR has been confirmed to work with the following combinations of browsers and operating systems:
| OS      | Version      | Chrome | Firefox | Microsoft Edge | Safari |
| ------- | ------------ | ------ | ------- | -------------- | ------ |
| Linux   | Ubuntu 20.04 | 107.0  | 106.0   | n/a            | n/a    |
| Windows | 10           | 107.0  | 107.0   | 107.0          | n/a    |
| MacOS   | Monterey     | 107.0  | 105.0   | 106.0          | 15.6   |

## Usage
A step-by-step manual and detailed information about the algorithms of the program will be published soon.

## Frequently asked questions

*Is repDilPCR only usable with raw Cq data obtained from experiments according to the dilution-replicate design?*

No, you can also feed precalculated relative quantities into the program. Of course, in this case you will only be able to use the statistical and graphical functions of repDilPCR.

*repDilPCR seems nice and I would like to try it but I have doubts about the dilution-replicate method. Does it give comparable results to the traditional methods?*

In our lab, we have directly compared the dilution-replicate approach with two of the most widely used traditional methods, the standard curve method and LinRegPCR (Ruijter et al., Amplification efficiency: linking baseline and bias in the analysis of quantitative PCR data. Nucleic Acids Research 2009, 37:6, e45, https://doi.org/10.1093/nar/gkp045), using the same samples. Generally, the three methods yielded very similar results but the standard curve method and repDilPCR outperformed LinRegPCR in an assay for miRNAs as the amplification curves in this assay had lower plateau and shorter exponential part which was not sufficient for LinRegPCR to reliably determine the so called window of linearity that this method relies upon for calculating reaction efficiences.

*Do I need technical replicates with the dilution-replicate method?*

The dilution-replicate method is based on so called dilution replicates. They serve both as a means to produce standard curves and to control technical variance. However, they are not called "technical replicates" as the concentration of the template in them is not identical. An estimation for the technical repeatability in an assay can be obtained from the standard curves plot and the goodness of fit (R<sup>2</sup> and adjusted R<sup>2</sup>). High values of adjusted R<sup>2</sup> (close to 1) speak for low technical variance in the assay as a whole. Single samples with high technical variance can be identified on the plot as one or more of their data points would lie away from the respective  regression lines. Such samples should be regarded with caution and excluded from the analysis if they deviate too much.

*Are standard curves with only 3 dilution points reliable enough to estimate the efficiency of the reaction?*

Although each individual standard curve would be based on only 3 dilutions, the multiple linear regression with parallel slopes would be based on the dilution points of all samples. As the slope is constrained, it effectively acts like an extra point in each sample. The degrees of freedom for the fit are given by the formula [(dilution points − 1) × (number of samples) − 1]. Thus, the more the samples, the better the precision. See the original publication of the method for more information (Kwokyin Hui & Zhong-Ping Feng. Efficient experimental design and analysis of real-time PCR assays. Channels 2013, 7:3, 160-170, https://doi.org/10.4161/chan.24024).

*Do the error bars on plots reflect technical variance, biological variance or both?*

Technical variance in a qPCR experiment is usually much smaller than variance between biological replicates and repDilPCR simply ignores it when calculating the variance in samples or experimental groups. You will see error bars only if you had biological replicates in your experiment (see "Preparation of the data" for information on how to name biological replicates so that the program will recognize them as such). You can get an impression of the technical variance by looking at the standard curves plots and the statistics below them (see the previous question). If you are specifically interested in determining technical variance with the method, you should perform replicates of the whole dilution series for each sample, which replicates you can name according to the convention for biological replicates (see "Preparation of the data").

*What do the error bars show? I see they are sometimes not symmetric, is this an error?*

For plots in logarithmic scale, the error bars denote the standard deviations and are symmetric. For plots in linear scale, the error bars show the 95%-confidence intervals. Confidence intervals are calculated from standard deviations in logarithmic scale and then converted to the linear scale. It is perfectly normal if they look asymmetric around the mean value.

*Why so complicated? Can't error bars depict standard deviation also on plots in linear scale?*

qPCR data are not normally distributed in linear scale. In this case, showing confidence intervals is more informative.

*If data are not normally distributed in linear scale, then why does repDilPCR use parametric statistical tests that assume normal distribution of the data?*

All statistical tests in repDilPCR are performed on log-transformed data, even when for display purposes you choose to prepare a plot in linear scale.

*Can repDilPCR perform two-way ANOVA?*

No. A leading concept in the design of the program was to make it simple (but sound) and offering high level of automation. Experiments with higher number of explanatory variables are too versatile and writing an algorithm that would be able to select an appropriate statistical test in each possible case would be too complicated. However, this should not stop you from using the dilution-replicate design and repDilPCR for the initial processing of the data as  repDilPCR allows downloading of all intermediate and final results in the open and widely compatible CSV format, so further processing with specialized statistical programs or popular spreadsheet software will not be a problem.

*Significance bars on my plots overlap and p-values are unreadable. What can I do?*

If there are a lot of experimental groups and a lot of the comparisons are significantly different, repDilPCR's algorithm may fail to prevent overlapping of significance bars. You can influence the algorithm by changing the spacing factor under `Distance between significance bars on plots` in the Shiny app (the equivalent variable in the R script is called `sp.f`). Increasing it will increase the distance between significance bars. Conversely, if the distances between significance bars are too large and they are wasting space on plots, you can try decreasing the spacing factor.

*Where does the name repDilPCR come from?*

Of course, the name and the logo of the program are allusion to the fact that developing the program was sponsored by the reptilians with the aim to subject mankind by enabling wider usage of the qPCR method, which has already been used by them to spread belief in the existence of a pandemic caused by the non-existing virus SARS-CoV-2. More sinister things are about to come. Watch out and wear your tin foil hat!