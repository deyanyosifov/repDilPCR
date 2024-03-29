## Title: repDilPCR - an R Script to Analyze qPCR Data by the Dilution-replicate Method
## File name: repDilPCR_CLI.R
## Version: 1.1.3
## Date: 2024-03-09
## Author: Deyan Yordanov Yosifov
## Maintainer: Deyan Yordanov Yosifov <deyan.yosifov@uniklinik-ulm.de>
## Copyright: University Hospital Ulm, Germany, 2021
## License: GNU General Public License v3 (GPL-3), https://www.gnu.org/licenses/gpl-3.0.html

## Variables
# Mandatory
input.table <- c("~/repDilPCR/Test_data.csv") ## Full path to the Cq-table
RG <- 3 ## Number of reference genes

# Optional
impute <- TRUE # Enter TRUE or FALSE to respectively enable or disable imputation of missing Cq values of reference genes. Do not set to TRUE if your experiment does not contain replicates.
ref.sample <- "default" # Reference sample, in which gene expression will be regarded as 1 (100%) on linear scale and 0 on log2-scale, respectively. If "default", this will be the first sample in the table, resp. the leftmost sample on the plots. Change to the name of another sample (without a trailing _ and replicate number) to make it the reference sample. If "" (empty), results will be shown in their original form, without forcing any particular sample to be 1 (100%) or 0.
statistics <- TRUE # Enter TRUE or FALSE to respectively enable or disable tests for statistical significance between samples or sample groups
test.type <- "parametric" # Type of statistical tests - "parametric" or "non-parametric". Default: "parametric"
posthoc <- "all to one" # Comparisons to test for statistically significant differences. Possible values: "all to one", "all pairs", "selected pairs". Default: "all to one".
alpha <- 0.05 # Significance level (alpha). Default value: 0.05
font.size <- 9 # Size of text on plots. Default value: 9
sign.repr <- "values" # Choose whether to display statistical significance with numeric p-values or with asterisks (significance levels). Possible values: "values" and "asterisks". Default: "values".
sp.f <- 1.5 # Spacing factor influencing the distance between significance bars on plots. In most cases, repDilPCR will succeed to distribute significance bars so that they will not overlap. If your significance bars overlap (which can be the case if you compare a lot of experimental groups), you can try increasing this factor. Conversely, if the distances between significance bars are too big and they are wasting space on plots, you can try decreasing the factor. The default value is 1.5.
colour.scheme <- "default" # Choose a colour scheme. Change the default colour scheme (scale_fill_hue()) by replacing the word "default" with the name of a ColorBrewer palette, e.g. Set1, Set2, Set3, Dark2 or Paired, or with the word "Grayscale". Default value: "default".
plot.format <- "PDF" # Choose format and adjust settings of graphical output. Possible formats are "PDF" (default), "PNG", "both" or "none".
png.size <- c(190,134) # Width and height of PNG plots in mm. The default values are 190 and 134 and fit well to the A5 page size format in landscape orientation.
png.dpi <- 96 # Resolution of PNG plots in dpi (default value: 96).


#----------------- Do not edit below this line unless you know what you are doing! ---------------
## Load packages
# source("~/Shiny/repDilPCR/repDilPCR_lib.R")
# source(paste0(getwd(),"/repDilPCR_script.R"))
source("repDilPCR_lib.R")


## Prepare data
if (length(str_extract_all(read_delim(input.table, ","), "\\d;\\d")[[1]]) > 5) {
  qPCR <- read.csv(input.table, sep = ";", dec = ",", stringsAsFactors = TRUE)
} else {  
  qPCR <- read.csv(input.table, sep = ",", dec = ".", stringsAsFactors = TRUE)
}
if (colnames(qPCR)[3] == "Dilution") {
  inp.data <- rd.preprocess(qPCR, RG)
  qPCR <- inp.data$qPCR
  ref.genes <- inp.data$ref.genes
  all.genes <- inp.data$all.genes
  GOIs <- inp.data$GOIs
} else {
    inp.data <- rd.preprocess.2(qPCR)
    qPCR <- inp.data$qPCR
    GOIs <- inp.data$GOIs
}

if (colnames(qPCR)[3] == "Dilution") {
  ## Imputation of missing Cq values of reference genes (optional step)
  if (impute == TRUE) {
    qPCR <- rd.impute(qPCR, ref.genes)
    write.csv(qPCR, file=paste0(gsub(".csv", "", input.table), "_with_imputed_missing_values.csv"))
  }

  ## Make a fictituous reference gene by averaging the Cq values of the real reference genes at each sample/dilution/replicate combination
  qPCR <- rd.ref(qPCR, RG)

  ## Multiple linear regression for standard curves
  model.list <- rd.mlr.model(qPCR, all.genes)

  ## Efficiency table
  eff.df <- rd.eff(model.list, all.genes)
  write.csv(eff.df, file=paste0(gsub(".csv", "", input.table), "_efficiency_table.csv"))

  ## Prepare colour scales
  colour_scale_regr <- rd.col.scale.regr(qPCR = qPCR, colour.scheme = colour.scheme, posthoc = posthoc)
  colour_scale <- rd.col.scale(qPCR = qPCR, colour.scheme = colour.scheme, posthoc = posthoc)
  
  ## Plot multiple regressions with separate regression curves but common slope
  stand.curves.output <- rd.plot.mlr(qPCR = qPCR, model.list = model.list, eff.df = eff.df, all.genes = all.genes, font.size = font.size, colour_scale_regr = colour_scale_regr)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_standard_curves.pdf"), width = 210, height = 297, units = "mm", stand.curves.output$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.1(ggplot.object = stand.curves.output$stand.curves, fname = input.table, type.name = "standard_curves", png.size = png.size, png.dpi = png.dpi)
  }
  rm(stand.curves.output)

  ## Multiple linear regression for Cq-Cq plots
  Cq.list <- rd.Cq.Cq(qPCR = qPCR, GOIs = GOIs)

  ## Plot Cq-Cq plots
  Cq.plots.output <- rd.plot.Cq.Cq(qPCR = qPCR, Cq.list = Cq.list, GOIs = GOIs, font.size = font.size, colour_scale_regr = colour_scale_regr)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_Cq-Cq_plots.pdf"), width = 210, height =  297, units = "mm", Cq.plots.output$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.1(ggplot.object = Cq.plots.output$Cq.plots, fname = input.table, type.name = "Cq-Cq_plots", png.size = png.size, png.dpi = png.dpi)
  }
  rm(Cq.plots.output)

  ## Relative quantities
  rel.q.results <- rd.rel.quant(qPCR = qPCR, Cq.list = Cq.list, eff.df = eff.df, GOIs = GOIs)
} else {
  ## Prepare colour scales
  colour_scale <- rd.col.scale(qPCR = qPCR, colour.scheme = colour.scheme, posthoc = posthoc)
  ## Relative quantities
  rel.q.results <- rd.rel.quant.2(qPCR = qPCR, GOIs = GOIs)
}

## Logarithmic tables
rel.q.results.log <- rd.log(rel.q.detailed = rel.q.results$rel.q.detailed, rel.q.df = rel.q.results$rel.q.df)

# Normalize relative quantities according to expression levels in a selected sample (if any)
rel.q.norm.results <- c(rel.q.results, rel.q.results.log)
rel.q.norm.results <- rd.normalize(rel.q.detailed = rel.q.norm.results$rel.q.detailed, rel.q.detailed.log = rel.q.norm.results$rel.q.detailed.log, rel.q.df = rel.q.norm.results$rel.q.df, rel.q.log = rel.q.norm.results$rel.q.log, rel.q.mean = rel.q.norm.results$rel.q.mean, rel.q.mean.log = rel.q.norm.results$rel.q.mean.log, ref.sample = ref.sample, GOIs = GOIs)
noref.warn <- c("A valid name of a sample to be used as baseline reference was not provided! Calculated relative quantities are not normalized to a particular sample.\n")
missingref.warn <- c("The chosen reference sample had missing data for one or more genes of interest! The respective plots will be empty and the result tables will not contain data for these genes. Choose a different reference sample or none to process and display all data.\n")
nonorm.warn <- c("The chosen reference sample had missing data for one or more genes of interest! Statistical tests have not been performed! Please choose a different reference sample or none.\n")


## Calculate confidence intervals
rel.q.confint <- rd.confint(rel.q.mean = rel.q.norm.results$rel.q.mean, rel.q.mean.log = rel.q.norm.results$rel.q.mean.log, p = alpha)

## Statistical tests
statistics.results <- rel.q.norm.results
statistics.results$rel.q.mean <- rel.q.confint$rel.q.mean
statistics.results$rel.q.mean.log <- rel.q.confint$rel.q.mean.log
statistics.results <- rd.statistics(rel.q.df = statistics.results$rel.q.df, rel.q.log = statistics.results$rel.q.log, rel.q.mean = statistics.results$rel.q.mean, rel.q.mean.log = statistics.results$rel.q.mean.log, statistics = statistics, test.type = test.type, posthoc = posthoc, ref.sample = statistics.results$ref.sample, nonorm = statistics.results$nonorm, p = alpha, sp.f = sp.f)
nostatref.warn <- c("A reference group for the Dunnett post-hoc test has not been chosen. Please check your input.\n")

## Plot relative expression
p1.results <- rd.plot.p1(rel.q.df = statistics.results$rel.q.df, rel.q.mean = statistics.results$rel.q.mean, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
if (plot.format == "both" | plot.format == "PDF") {
  ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_expression_dotplot.pdf"), width = 210, height =  297, units = "mm", p1.results$ml)
}
if (plot.format == "both" | plot.format == "PNG") {
  png.plot.2(ggplot.object = p1.results$p1, fname = input.table, type.name = "relative_expression_dotplot", png.size = png.size, png.dpi = png.dpi)
}
rm(p1.results)

if (test.type == "parametric") {
  p2.results <- rd.plot.p2(rel.q.mean = statistics.results$rel.q.mean, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_expression_dotplot_CI.pdf"), width = 210, height =  297, units = "mm", p2.results$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.2(ggplot.object = p2.results$p2, fname = input.table, type.name = "relative_expression_dotplot_CI", png.size = png.size, png.dpi = png.dpi)
  }
  p3.results <- rd.plot.p3(rel.q.mean = statistics.results$rel.q.mean, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_expression_bar_graph.pdf"), width = 210, height =  297, units = "mm", p3.results$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.2(ggplot.object = p3.results$p3, fname = input.table, type.name = "relative_expression_bar_graph", png.size = png.size, png.dpi = png.dpi)
  }
  rm(p2.results)
  rm(p3.results)
}

if (test.type == "non-parametric") {
  p2n.results <- rd.plot.p2n(rel.q.df = statistics.results$rel.q.df, rel.q.mean = statistics.results$rel.q.mean, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_expression_boxplot.pdf"), width = 210, height =  297, units = "mm", p2n.results$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.2(ggplot.object = p2n.results$p2n, fname = input.table, type.name = "relative_expression_boxplot", png.size = png.size, png.dpi = png.dpi)
  }
  rm(p2n.results)
}

p4.results <- rd.plot.p4(rel.q.log = statistics.results$rel.q.log, rel.q.mean = statistics.results$rel.q.mean, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
if (plot.format == "both" | plot.format == "PDF") {
  ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_log_expression_dotplot.pdf"), width = 210, height =  297, units = "mm", p4.results$ml)
}
if (plot.format == "both" | plot.format == "PNG") {
  png.plot.2(ggplot.object = p4.results$p4, fname = input.table, type.name = "relative_log_expression_dotplot", png.size = png.size, png.dpi = png.dpi)
}
rm(p4.results)

if (test.type == "parametric") {
  p5.results <- rd.plot.p5(rel.q.mean.log = statistics.results$rel.q.mean.log, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_log_expression_dotplot_SD.pdf"), width = 210, height =  297, units = "mm", p5.results$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.2(ggplot.object = p5.results$p5, fname = input.table, type.name = "relative_log_expression_dotplot_SD", png.size = png.size, png.dpi = png.dpi)
  }
  p6.results <- rd.plot.p6(rel.q.mean.log = statistics.results$rel.q.mean.log, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_log_expression_bar_graph.pdf"), width = 210, height =  297, units = "mm", p6.results$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.2(ggplot.object = p6.results$p6, fname = input.table, type.name = "relative_log_expression_bar_graph", png.size = png.size, png.dpi = png.dpi)
  }
  rm(p5.results)
  rm(p6.results)
}

if (test.type == "non-parametric") {
  p5n.results <- rd.plot.p5n(rel.q.log = statistics.results$rel.q.log, rel.q.mean = statistics.results$rel.q.mean, res.posthoc = statistics.results$res.posthoc, ref.sample = statistics.results$ref.sample, GOIs = GOIs, statistics = statistics.results$statistics, posthoc = posthoc, sign.repr = sign.repr, p = alpha, stat.test = statistics.results$stat.test, font.size = font.size, nonorm = rel.q.norm.results$nonorm, colour_scale = colour_scale)
  if (plot.format == "both" | plot.format == "PDF") {
    ggplot2::ggsave(paste0(gsub(".csv", "", input.table), "_relative_log_expression_boxplot.pdf"), width = 210, height =  297, units = "mm", p5n.results$ml)
  }
  if (plot.format == "both" | plot.format == "PNG") {
    png.plot.2(ggplot.object = p5n.results$p5n, fname = input.table, type.name = "relative_log_expression_boxplot", png.size = png.size, png.dpi = png.dpi)
  }
  rm(p5n.results)
}

## Save tables
save.tables <- rd.save.tables(input.table = input.table, rel.q.detailed = rel.q.norm.results$rel.q.detailed, rel.q.detailed.log = rel.q.norm.results$rel.q.detailed.log, rel.q.mean = statistics.results$rel.q.mean, rel.q.mean.log = statistics.results$rel.q.mean.log, p = alpha)

## Print warning messages if any
if (exists("statistics.results")) {
  rd.warnings <- rd.warn(ref.sample = statistics.results$ref.sample, rel.q.mean = statistics.results$rel.q.mean, noref.warn = noref.warn, statistics = statistics.results$statistics, posthoc = posthoc, nostatref.warn = nostatref.warn, frw = statistics.results$frw, few.repl.warn = statistics.results$few.repl.warn, rel.q.mean.log = rel.q.results.log$rel.q.mean.log, missingref.warn = missingref.warn, nonorm = rel.q.norm.results$nonorm, nonorm.warn = nonorm.warn, sel.pairs = statistics.results$sel.pairs, sel.pairs.warn = statistics.results$sel.pairs.warn)
}
rd.warnings.2 <- rd.warn.2(csv.wrong.format = inp.data$csv.wrong.format)
