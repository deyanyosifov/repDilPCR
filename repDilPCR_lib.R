## Title: A Library of Functions Used by both the repDilPCR R Script and the repDilPCR Shiny App
## File name: repDilPCR_lib.R
## Version: 1.0.3
## Date: 2021-08-09
## Author: Deyan Yordanov Yosifov
## Maintainer: Deyan Yordanov Yosifov <deyan.yosifov@uniklinik-ulm.de>
## Copyright: University Hospital Ulm, Germany, 2021
## License: GNU General Public License v3 (GPL-3), https://www.gnu.org/licenses/gpl-3.0.html

## Variables
# Mandatory
# input.table <- c("/path/to/file.csv") ## Full path to the Cq-table
# RG <- 3 ## Number of reference genes
#
# # Optional
# impute <- TRUE # Enter TRUE or FALSE to respectively enable or disable imputation of missing Cq-values of reference genes
# ref.sample <- "default" # Reference sample, in which gene expression will be regarded as 1 (100%) on linear scale and 0 on log2-scale, respectively. If "default", this will be the first sample in the table, resp. the leftmost sample on the plots. Change to the name of another sample (without a trailing _ and replicate number) to make it the reference sample. If "" (empty), results will be shown in their original form, without forcing any particular sample to be 1 (100%) or 0.
# plot.format <- "both" # Choose format and adjust settings of graphical output. Possible formats are "PDF" (default), "PNG", "both" or "none".
# png.size <- c(190,134) # Width and height of PNG plots in mm. The default values are 210 and 148 and correspond to the A5 page size format in landscape orientation.
# png.dpi <- 96 # Resolution of PNG plots in dpi (default value: 96).
# statistics <- TRUE # Enter TRUE or FALSE to respectively enable or disable tests for statistical significance between samples or sample groups
# test.type <- "parametric" # Type of statistical tests - "parametric" or "non-parametric". Default: "parametric"
# posthoc <- "all pairs" # Comparisons to test for statistically significant differences. Possible values: "all to one", "all pairs", "selected pairs". Default: "all to one".
# p <- 0.05 # Significance level (alpha). Default value: 0.05
# sign.repr <- "values" # Choose whether to display statistical significance with numeric p-values or with asterisks (significance levels). Possible values: "values" and "asterisks". Default: "values".
# sp.f <- 1.5 # Spacing factor influencing the distance between significance bars on plots. In most cases, repDilPCR will succeed to distribute significance bars so that they will not overlap. If your significance bars overlap (which can be the case if you compare a lot of experimental groups), you can try increasing this factor. Conversely, if the distances between significance bars are too big and they are wasting space on plots, you can try decreasing the factor. The default value is 1.5.
#
## Load packages
library(car)
library(gridExtra)
library(tidyverse)
library(mice)
library(PMCMRplus)
library(ggbeeswarm)
library(ggsignif)

## Convenience functions for formatting numbers and plotting to PNG files
get.rid.zeros <- function (x) ifelse(x > .Machine$double.eps,identity(x),x <- .Machine$double.eps)

png.plot.1 <- function(ggplot.object, fname, type.name, png.size, png.dpi) {
  path <- gsub(".csv", "", fname)
  fs <- paste0(path, "_", type.name, "_legend.png")
  w <- unlist(ggplot.object[[1]]$widths[[3]])
  if (length(w) < 3) {
    w <- as.numeric(w[1])
  } else {
    w <- as.numeric(w[3:length(w)])
  }
  h <- unlist(ggplot.object[[1]]$heights[[3]])
  if (length(h) < 3) {
    h<- as.numeric(h[1])
  } else {
    h <- as.numeric(h[3:length(h)])
  }
  # png(file = fs, width = as.numeric(ggplot.object[[1]]$widths[[3]]), height = as.numeric(ggplot.object[[1]]$heights[[3]]), units = "cm", bg = "white", res=png.dpi, type = "cairo", antialias = "default")
  png(file = fs, width = sum(w), height = sum(h), units = "cm", bg = "white", res=png.dpi, type = "cairo", antialias = "default")
  plot(ggplot.object[[1]])
  dev.off()
  for (i in 2:length(ggplot.object)) {
    fs2 <- paste0(path, "_", type.name, "_", ggplot.object[[i]]$labels$title, ".png")
    png(file = fs2, width = png.size[1], height = png.size[2], units = "mm", bg = "white", res=png.dpi, type = "cairo", antialias = "default")
    plot(ggplot.object[[i]])
    dev.off()
    fs <- c(fs, fs2)
  }
  return(fs)
}


png.plot.2 <- function(ggplot.object, fname, type.name, png.size, png.dpi) {
  path <- gsub(".csv", "", fname)
  fs <- c()
  for (i in 1:length(ggplot.object)) {
    fs2 <- paste0(path, "_", type.name, "_", ggplot.object[[i]]$labels$title, ".png")
    png(file = fs2, width = png.size[1], height = png.size[2], units = "mm", bg = "white", res=png.dpi, type = "cairo", antialias = "default")
    plot(ggplot.object[[i]])
    dev.off()
    fs <- c(fs, fs2)
  }
  return(fs)
}


## Prepare data
# qPCR <- read.csv(input.table, sep = ",", dec = ".", stringsAsFactors = TRUE) ## input data file name (insert it between the quotation marks)
rd.preprocess <- function(qPCR, RG) {
  num.genes <- ncol(qPCR) - 3
  col.names <- colnames(qPCR)
  gene.names <- col.names[4:ncol(qPCR)]
  ref.genes <- gene.names[1:RG]
  GOIs <- gene.names[(RG+1):num.genes]
  all.genes <- c(GOIs, "NF", ref.genes)
  qPCR$Samples <- as.factor(gsub("_.*","",qPCR$Replicates))
  qPCR$Replicates <- factor(qPCR$Replicates, levels = unique(qPCR$Replicates))
  qPCR$Samples <- factor(qPCR$Samples, levels = unique(qPCR$Samples))
  preprocess.output <- list("qPCR" = qPCR, "ref.genes" = ref.genes, "all.genes" = all.genes, "GOIs" = GOIs)
  return(preprocess.output)
}

rd.preprocess.2 <- function(qPCR) {
  rownames(qPCR) <- qPCR$Replicates
  # qPCR <- qPCR[,2:ncol(qPCR)]
  GOIs <- colnames(qPCR)[3:ncol(qPCR)]
  qPCR$Samples <- as.factor(gsub("_.*","",qPCR$Replicates))
  qPCR$Replicates <- factor(qPCR$Replicates, levels = unique(qPCR$Replicates))
  qPCR$Samples <- factor(qPCR$Samples, levels = unique(qPCR$Samples))
  preprocess.output <- list("qPCR" = qPCR, "GOIs" = GOIs)
  return(preprocess.output)
}

## Imputation of missing Cq values of reference genes (optional step)
rd.impute <- function(qPCR, ref.genes) {
set.seed(40075017)
imputed_data <- mice::mice(qPCR[,c(ref.genes, "Samples", "Dilution")], m = 20, maxit = 20, method = 'midastouch', print = F)
for (i in ref.genes) {
  a <- rowMeans(imputed_data$imp[[i]])
  qPCR[names(a),i] <- a
}
return(qPCR)
# write.csv(qPCR, file=paste0(gsub(".csv", "", input.table), "_with_imputed_missing_values.csv"))
}

## Make a fictituous reference gene by averaging the Cq values of the real reference genes at each sample/dilution/replicate combination
rd.ref <- function(qPCR, RG) {
  qPCR$NF <- apply(as.data.frame(qPCR[,seq(RG)+3], stringsAsFactors = FALSE),1, mean)
  return(qPCR)
}


## Multiple linear regression for standard curves
rd.mlr.model <- function (qPCR, all.genes) {
  model.list <- vector("list", length(all.genes))
  names(model.list) <- all.genes
  for (i in 1:length(all.genes)) {
    a <- all.genes[i]
    b <- "log(1/Dilution) + Replicates"
    fml <- paste(a, b, sep="~")
    model.list[[i]] <- lm(fml, data = qPCR)
  }
  return(model.list)
}

## Efficiency table
rd.eff <- function (model.list, all.genes) {
  efficiencies <- list()
  efficiencies.left.conf.int <- list()
  efficiencies.right.conf.int <- list()
  for (i in 1:length(all.genes)) {
    efficiencies[[i]] <- exp(1)^-(1/model.list[[i]]$coefficients[2])
    efficiencies.left.conf.int[[i]] <- exp(1)^-(1/confint(model.list[[i]])[2,1])
    efficiencies.right.conf.int[[i]] <- exp(1)^-(1/confint(model.list[[i]])[2,2])
  }
  eff <- unlist(efficiencies)
  eff.left.conf.int <- unlist(efficiencies.left.conf.int)
  eff.right.conf.int <- unlist(efficiencies.right.conf.int)
  eff.df <- as.data.frame(t(rbind(eff, eff.left.conf.int, eff.right.conf.int)), stringsAsFactors = FALSE)
  rownames(eff.df) <- all.genes
  colnames(eff.df) <- c("Efficiency", "Lower bound of the 95% confidence interval", "Upper bound of the 95% confidence interval")
  return(eff.df)
}
# write.csv(eff.df, file=paste0(gsub(".csv", "", input.table), "_efficiency_table.csv"))


## Plot multiple regressions with separate regression curves but common slope
# Variant for plotting from a script
rd.plot.mlr <- function(qPCR, model.list, eff.df, all.genes, font.size) {
  pred.data <- list()
  for (i in 1:length(all.genes)) {
    pred.data[[i]] <- predict(model.list[[i]])
  }
  max.length <- nrow(qPCR)
  d <- length(pred.data)
  for (i in 1:d) {
    pred.data[[d+1]] <- c(rep(NA, max.length))
    names(pred.data[[d+1]]) <- rownames(qPCR)
    for (j in names(pred.data[[i]])) {
      pred.data[[d+1]][j] <- pred.data[[i]][j]
    }
    pred.data[[i]] <- pred.data[[d+1]]
    pred.data[[d+1]] <- NULL
  }
  pred.df <- data.frame(pred.data, stringsAsFactors = FALSE)
  pred.names <- paste0(all.genes,"_pred")
  colnames(pred.df) <- pred.names

  stand.curves <- list()

  legend <- ggplot2::ggplot(data=cbind(qPCR, pred.df), ggplot2::aes_string(x = "log(1/Dilution)", y = "NF", colour="Replicates", fill="Replicates")) +
    ggplot2::geom_point(size=2) +
    ggplot2::geom_line(ggplot2::aes_(y = as.name("NF_pred"))) +
    ggplot2::theme_bw() +
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=4), fill=ggplot2::guide_legend(ncol=4)) +
    ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
    ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  stand.curves[[1]] <- ggplot2::ggplot_build(legend)
  legend.gtable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(legend))
  leg <- which(sapply(legend.gtable$grobs, function(x) x$name) == "guide-box")
  legend.grob <- legend.gtable$grobs[[leg]]
  stand.curves[[1]] <- legend.grob

  for (i in 1:length(all.genes)) {
    stand.curves[[i+1]] <- ggplot2::ggplot(data=cbind(qPCR, pred.df), ggplot2::aes_string(x = "log(1/Dilution)", y = all.genes[i], colour="Replicates", fill="Replicates")) +
      ggplot2::geom_point(size=2) +
      ggplot2::geom_line(ggplot2::aes_(y = as.name(pred.names[i]))) +
      ggplot2::xlim(floor(min(log(1/qPCR$Dilution))), ceiling(max(log(1/qPCR$Dilution)))) +
      ggplot2::ylim(floor(min(pred.df[,i], na.rm = TRUE)), ceiling(max(pred.df[,i], na.rm = TRUE))) +
      ggplot2::xlab("ln(1/dilution)") +
      ggplot2::ylab("Cq") +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::ggtitle(all.genes[i]) +
      ggplot2::geom_text(x=-0.75, y=floor(max(pred.df[,i]-1, na.rm = TRUE)), label = paste0("Eff = ", round(eff.df[i,1], 2), " (", round(eff.df[i,2], 2), "-", round(eff.df[i,3], 2), ")"), color="gray20", size=4) +
      ggplot2::labs(caption = gsub("^.*(Residual standard.*$)", "\\1", paste(capture.output(summary(model.list[[i]])), collapse = "\n"))) +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1))
  }
  ml <- gridExtra::marrangeGrob(grobs = stand.curves, nrow = 2, ncol = 1)
  stand.curves.output <- list("stand.curves" = stand.curves, "ml" = ml)
  return(stand.curves.output)
}

# Variant for plotting from a Shiny app
rd.plot.mlr.s <- function(qPCR, model.list, eff.df, all.genes, font.size) {
  pred.data <- list()
  for (i in 1:length(all.genes)) {
    pred.data[[i]] <- predict(model.list[[i]])
  }
  max.length <- nrow(qPCR)
  d <- length(pred.data)
  for (i in 1:d) {
    pred.data[[d+1]] <- c(rep(NA, max.length))
    names(pred.data[[d+1]]) <- rownames(qPCR)
    for (j in names(pred.data[[i]])) {
      pred.data[[d+1]][j] <- pred.data[[i]][j]
    }
    pred.data[[i]] <- pred.data[[d+1]]
    pred.data[[d+1]] <- NULL
  }
  pred.df <- data.frame(pred.data, stringsAsFactors = FALSE)
  pred.names <- paste0(all.genes,"_pred")
  colnames(pred.df) <- pred.names

  stand.curves <- list()

  for (i in 1:length(all.genes)) {
    stand.curves[[i]] <- ggplot2::ggplot(data=cbind(qPCR, pred.df), ggplot2::aes_string(x = "log(1/Dilution)", y = all.genes[i], colour="Replicates", fill="Replicates")) +
      ggplot2::geom_point(size=2) +
      ggplot2::geom_line(ggplot2::aes_(y = as.name(pred.names[i]))) +
      ggplot2::xlim(floor(min(log(1/qPCR$Dilution))), ceiling(max(log(1/qPCR$Dilution)))) +
      ggplot2::ylim(floor(min(pred.df[,i], na.rm = TRUE)), ceiling(max(pred.df[,i], na.rm = TRUE))) +
      ggplot2::xlab("ln(1/dilution)") +
      ggplot2::ylab("Cq") +
      ggplot2::ggtitle(all.genes[i]) +
      ggplot2::geom_text(x=-0.75, y=floor(max(pred.df[,i]-1, na.rm = TRUE)), label = paste0("Eff = ", round(eff.df[i,1], 2), " (", round(eff.df[i,2], 2), "-", round(eff.df[i,3], 2), ")"), color="gray20", size=5) +
      ggplot2::labs(caption = gsub("^.*(Residual standard.*$)", "\\1", paste(capture.output(summary(model.list[[i]])), collapse = "\n"))) +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=4), fill=ggplot2::guide_legend(ncol=4)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="bottom") +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  }
  return(stand.curves)
}

## Multiple linear regression for Cq-Cq plots
rd.Cq.Cq <- function(qPCR, GOIs) {
  Cq.list <- vector("list", length(GOIs))
  names(Cq.list) <- GOIs
  for (i in 1:length(GOIs)) {
    a <- GOIs[i]
    b <- ("NF + Replicates")
    fml <- paste(a, b, sep="~")
    Cq.list[[i]] <- lm(fml, data = qPCR)
  }
  return(Cq.list)
}

## Plot Cq-Cq plots
# Variant for plotting from a script
rd.plot.Cq.Cq <- function(qPCR, Cq.list, GOIs, font.size) {
  pred.Cq <- list()
  for (i in 1:length(Cq.list)) {
    pred.Cq[[i]] <- predict(Cq.list[[i]])
  }
  max.length <- nrow(qPCR)
  d <- length(pred.Cq)
  for (i in 1:d) {
    pred.Cq[[d+1]] <- c(rep(NA, max.length))
    names(pred.Cq[[d+1]]) <- rownames(qPCR)
    for (j in names(pred.Cq[[i]])) {
      pred.Cq[[d+1]][j] <- pred.Cq[[i]][j]
    }
    pred.Cq[[i]] <- pred.Cq[[d+1]]
    pred.Cq[[d+1]] <- NULL
  }
  pred.Cq.df <- data.frame(pred.Cq, stringsAsFactors = FALSE)
  pred.Cq.names <- paste0(GOIs,"_pred")
  colnames(pred.Cq.df) <- pred.Cq.names

  Cq.plots <- list()

  legend <- ggplot2::ggplot(data=cbind(qPCR, pred.Cq.df), ggplot2::aes_string(x = "NF", y = GOIs[1], colour="Replicates", fill="Replicates")) +
    ggplot2::geom_point(size=2) +
    ggplot2::geom_line(ggplot2::aes_(y = as.name(pred.Cq.names[1]))) +
    ggplot2::theme_bw() +
    ggplot2::guides(colour=ggplot2::guide_legend(ncol=4), fill=ggplot2::guide_legend(ncol=4)) +
    ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
    ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  legend.gtable <- ggplot2::ggplot_gtable(ggplot2::ggplot_build(legend))
  leg <- which(sapply(legend.gtable$grobs, function(x) x$name) == "guide-box")
  legend.grob <- legend.gtable$grobs[[leg]]
  Cq.plots[[1]] <- legend.grob

  for (i in 1:length(GOIs)) {
    Cq.plots[[i+1]] <- ggplot2::ggplot(data=cbind(qPCR, pred.Cq.df), ggplot2::aes_string(x = "NF", y = GOIs[i], colour="Replicates", fill="Replicates")) +
      ggplot2::geom_point(size=2) +
      ggplot2::geom_line(ggplot2::aes_(y = as.name(pred.Cq.names[i]))) +
      ggplot2::xlim(floor(min(qPCR$NF)), ceiling(max(qPCR$NF))) +
      ggplot2::ylim(floor(min(pred.Cq.df[,i])), ceiling(max(pred.Cq.df[,i]))) +
      # ggplot2::coord_fixed(ratio=1) +
      ggplot2::xlab(c("NF (Cq)")) +
      ggplot2::ylab(paste(GOIs[i],"(Cq)")) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="none") +
      ggplot2::ggtitle(GOIs[i]) +
      ggplot2::labs(caption = gsub("^.*(Residual standard.*$)", "\\1", paste(capture.output(summary(Cq.list[[i]])), collapse = "\n"))) +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1))
  }
  ml <- gridExtra::marrangeGrob(grobs = Cq.plots, nrow = 2, ncol = 1)
  Cq.plots.output <- list("Cq.plots" = Cq.plots, "ml" = ml)
  return(Cq.plots.output)
}

# Variant for plotting from a Shiny app
rd.plot.Cq.Cq.s <- function(qPCR, Cq.list, GOIs, font.size) {
  pred.Cq <- list()
  for (i in 1:length(Cq.list)) {
    pred.Cq[[i]] <- predict(Cq.list[[i]])
  }
  max.length <- nrow(qPCR)
  d <- length(pred.Cq)
  for (i in 1:d) {
    pred.Cq[[d+1]] <- c(rep(NA, max.length))
    names(pred.Cq[[d+1]]) <- rownames(qPCR)
    for (j in names(pred.Cq[[i]])) {
      pred.Cq[[d+1]][j] <- pred.Cq[[i]][j]
    }
    pred.Cq[[i]] <- pred.Cq[[d+1]]
    pred.Cq[[d+1]] <- NULL
  }
  pred.Cq.df <- data.frame(pred.Cq, stringsAsFactors = FALSE)
  pred.Cq.names <- paste0(GOIs,"_pred")
  colnames(pred.Cq.df) <- pred.Cq.names

  Cq.plots <- list()

  for (i in 1:length(GOIs)) {
    Cq.plots[[i]] <- ggplot2::ggplot(data=cbind(qPCR, pred.Cq.df), ggplot2::aes_string(x = "NF", y = GOIs[i], colour="Replicates", fill="Replicates")) +
      ggplot2::geom_point(size=2) +
      ggplot2::geom_line(ggplot2::aes_(y = as.name(pred.Cq.names[i]))) +
      ggplot2::xlim(floor(min(qPCR$NF)), ceiling(max(qPCR$NF))) +
      ggplot2::ylim(floor(min(pred.Cq.df[,i])), ceiling(max(pred.Cq.df[,i]))) +
      ggplot2::coord_fixed(ratio=1) +
      ggplot2::xlab(c("NF (Cq)")) +
      ggplot2::ylab(paste(GOIs[i],"(Cq)")) +
      ggplot2::ggtitle(GOIs[i]) +
      ggplot2::labs(caption = gsub("^.*(Residual standard.*$)", "\\1", paste(capture.output(summary(Cq.list[[i]])), collapse = "\n"))) +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=4), fill=ggplot2::guide_legend(ncol=4)) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position="bottom") +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  }
  return(Cq.plots)
}


## Relative quantities table
rd.rel.quant <- function(qPCR, Cq.list, eff.df, GOIs) {
  rel.quantities <- list()
  for (i in GOIs) {
    m2 <- m1 <- m <- list()
    for (n in 3:length(Cq.list[[i]]$coefficients)) {
      m[[n]] <- exp(1)^-((Cq.list[[i]]$coefficients[1]+Cq.list[[i]]$coefficients[n])*log(eff.df[i,1]))
      # m1[[n]] <- exp(1)^-((Cq.list[[i]]$coefficients[1]+confint(Cq.list[[i]])[n,2])*log(eff.df[i,1]))
      # m2[[n]] <- exp(1)^-((Cq.list[[i]]$coefficients[1]+confint(Cq.list[[i]])[n,1])*log(eff.df[i,1]))
    }
    v <- gsub("Replicates","",names(Cq.list[[i]]$coefficients)[c(-1,-2)])
    m <-  unlist(m)
    # m1 <-  unlist(m1)
    # m2 <-  unlist(m2)
    rel.quantities[[i]] <- c(exp(1)^-(Cq.list[[i]]$coefficients[1]*log(eff.df[i,1])), m)
    # rel.quantities[[paste0(i, ".left.CI")]] <- c(exp(1)^-(confint(Cq.list[[i]])[1,2]*log(eff.df[i,1])), m1)
    # rel.quantities[[paste0(i, ".right.CI")]] <- c(exp(1)^-(confint(Cq.list[[i]])[1,1]*log(eff.df[i,1])), m2)
    names(rel.quantities[[i]]) <- c(levels(qPCR$Replicates)[1],v)
    # names(rel.quantities[[paste0(i, ".left.CI")]]) <- c(levels(qPCR$Replicates)[1],v)
    # names(rel.quantities[[paste0(i, ".right.CI")]]) <- c(levels(qPCR$Replicates)[1],v)
  }
  lengths <- lapply(rel.quantities, length)
  max.length <- lengths[[which.max(lengths)]]
  d <- length(rel.quantities)
  for (i in 1:d) {
    rel.quantities[[d+1]] <- c(rep(NA, max.length))
    names(rel.quantities[[d+1]]) <- as.character(unique(qPCR$Replicates))
    for (j in 1:max.length) {
      rel.quantities[[d+1]][j] <- rel.quantities[[i]][as.character(unique(qPCR$Replicates))[j]]
    }
    rel.quantities[[i]] <- rel.quantities[[d+1]]
    rel.quantities[[d+1]] <- NULL
  }

  rel.q.detailed <- as.data.frame(rel.quantities, stringsAsFactors = FALSE)
  if(ncol(rel.q.detailed) > 1) {
    rel.q <- rel.q.detailed[,GOIs]
  } else  {
    rel.q <- rel.q.detailed
  }
  rel.q.df <- as.data.frame(cbind(as.vector(as.matrix(rel.q)), rep(rownames(rel.q),ncol(rel.q)), rep(colnames(rel.q), each=nrow(rel.q))), stringsAsFactors = FALSE)
  colnames(rel.q.df) <- c("Rel.quant", "Replicates", "Genes")
  rel.q.df$Samples <- gsub("_.*","",rel.q.df$Replicates)
  rel.q.df$GeneSampleCombo <- paste0(rel.q.df$Genes, " in ", rel.q.df$Samples)
  rel.q.df$Rel.quant <- as.numeric(rel.q.df$Rel.quant)
  rel.q.df$Pairs <- qPCR[match(rel.q.df$Replicates, qPCR$Replicates), "Pairs"]
  rel.q.mean <- unique(rel.q.df[,c("Genes", "Samples", "Pairs")])
  for (i in 1:nrow(rel.q.mean)) {
    gene <- rel.q.mean[i, "Genes"]
    s.name <- rel.q.mean[i, "Samples"]
    rel.q.mean$Replicates[i] <- length(na.omit(rel.q.df[which(rel.q.df$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"]))
    rel.q.mean$Expression[i] <- exp(mean(log(as.numeric(rel.q.df[which(rel.q.df$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"])), na.rm = TRUE))
    # rel.q.mean$SD[i] <- sd(as.numeric(rel.q.df[which(rel.q.df$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"]), na.rm = TRUE)
  }
  rel.q.df$Samples <- factor(rel.q.df$Samples, levels = unique(rel.q.df$Samples))
  rel.q.mean$Samples <- factor(rel.q.mean$Samples, levels = unique(rel.q.mean$Samples))
  rel.q.results <- list("rel.q.detailed" = rel.q.detailed, "rel.q.df" = rel.q.df, "rel.q.mean" = rel.q.mean)
  return(rel.q.results)
}

rd.rel.quant.2 <- function(qPCR, GOIs) {
  rel.q <- qPCR[,GOIs]
  rel.q.df <- as.data.frame(cbind(as.vector(as.matrix(rel.q)), rep(rownames(rel.q),ncol(rel.q)), rep(colnames(rel.q), each=nrow(rel.q))), stringsAsFactors = FALSE)
  colnames(rel.q.df) <- c("Rel.quant", "Replicates", "Genes")
  rel.q.df$Samples <- gsub("_.*","",rel.q.df$Replicates)
  rel.q.df$GeneSampleCombo <- paste0(rel.q.df$Genes, " in ", rel.q.df$Samples)
  rel.q.df$Rel.quant <- as.numeric(rel.q.df$Rel.quant)
  rel.q.df$Pairs <- qPCR$Pairs
  rel.q.mean <- unique(rel.q.df[,c("Genes", "Samples", "Pairs")])
  for (i in 1:nrow(rel.q.mean)) {
    gene <- rel.q.mean[i, "Genes"]
    s.name <- rel.q.mean[i, "Samples"]
    rel.q.mean$Replicates[i] <- length(na.omit(rel.q.df[which(rel.q.df$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"]))
    rel.q.mean$Expression[i] <- exp(mean(log(as.numeric(rel.q.df[which(rel.q.df$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"])), na.rm = TRUE))
    # rel.q.mean$SD[i] <- sd(as.numeric(rel.q.df[which(rel.q.df$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"]), na.rm = TRUE)
  }
  rel.q.df$Samples <- factor(rel.q.df$Samples, levels = unique(rel.q.df$Samples))
  rel.q.mean$Samples <- factor(rel.q.mean$Samples, levels = unique(rel.q.mean$Samples))
  rel.q.results <- list("rel.q.detailed" = rel.q, "rel.q.df" = rel.q.df, "rel.q.mean" = rel.q.mean)
  return(rel.q.results)
}

# Logarithmic tables
rd.log <- function(rel.q.detailed, rel.q.df) {
  rel.q.detailed.log <- log2(rel.q.detailed)
  rel.q.log <- rel.q.df
  rel.q.log$Rel.quant <- log2(rel.q.log$Rel.quant)
  rel.q.mean.log <- unique(rel.q.log[,c("Genes", "Samples", "Pairs")])
  for (i in 1:nrow(rel.q.mean.log)) {
    gene <- rel.q.mean.log[i, "Genes"]
    s.name <- rel.q.mean.log[i, "Samples"]
    rel.q.mean.log$Expression[i] <- mean(as.numeric(rel.q.log[which(rel.q.log$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"]), na.rm = TRUE)
    rel.q.mean.log$SD[i] <- sd(as.numeric(rel.q.log[which(rel.q.log$GeneSampleCombo==paste0(gene, " in ", s.name)),"Rel.quant"]), na.rm = TRUE)
  }
  rel.q.log$Samples <- factor(rel.q.log$Samples, levels = unique(rel.q.log$Samples))
  rel.q.mean.log$Samples <- factor(rel.q.mean.log$Samples, levels = unique(rel.q.mean.log$Samples))
  rel.q.results.log <- list("rel.q.detailed.log" = rel.q.detailed.log, "rel.q.log" = rel.q.log, "rel.q.mean.log" = rel.q.mean.log)
  return(rel.q.results.log)
}

# Normalize relative quantities according to expression levels in a selected sample (if any)
rd.normalize <- function(rel.q.detailed, rel.q.detailed.log, rel.q.df, rel.q.log, rel.q.mean, rel.q.mean.log, ref.sample, GOIs) {
  if (ref.sample == "default") {
    ref.sample <- rel.q.mean[1, "Samples"]
  }
  if (ref.sample %in% rel.q.mean$Samples) {
    norm.factor.log <- norm.factor <- numeric()
    for (i in unique(rel.q.mean$Genes)) {
      norm.factor[i] <- rel.q.mean[which(rel.q.mean$Genes == i & rel.q.mean$Samples == ref.sample),"Expression"]
      norm.factor.log[i] <- rel.q.mean.log[which(rel.q.mean.log$Genes == i & rel.q.mean.log$Samples == ref.sample),"Expression"]
    }
    for (i in GOIs) {
      rel.q.detailed[, i] <- rel.q.detailed[, i]/norm.factor[i]
      # rel.q.detailed[, paste0(i, ".left.CI")] <- rel.q.detailed[, paste0(i, ".left.CI")]/norm.factor[i]
      # rel.q.detailed[, paste0(i, ".right.CI")] <- rel.q.detailed[, paste0(i, ".right.CI")]/norm.factor[i]
      rel.q.detailed.log[, i] <- rel.q.detailed.log[, i] - norm.factor.log[i]
      # rel.q.detailed.log[, paste0(i, ".left.CI")] <- rel.q.detailed.log[, paste0(i, ".left.CI")] - norm.factor.log[i]
      # rel.q.detailed.log[, paste0(i, ".right.CI")] <- rel.q.detailed.log[, paste0(i, ".right.CI")] - norm.factor.log[i]
    }
    for (i in unique(rel.q.df$Genes)) {
      rel.q.df[which(rel.q.df$Genes == i), "Rel.quant"] <- rel.q.df[which(rel.q.df$Genes == i), "Rel.quant"]/norm.factor[i]
      rel.q.log[which(rel.q.log$Genes == i), "Rel.quant"] <- rel.q.log[which(rel.q.log$Genes == i), "Rel.quant"] - norm.factor.log[i]
    }
    for (i in unique(rel.q.mean$Genes)) {
      rel.q.mean[which(rel.q.mean$Genes == i), c("Expression")] <- rel.q.mean[which(rel.q.mean$Genes == i), c("Expression")]/norm.factor[i]
      rel.q.mean.log[which(rel.q.mean.log$Genes == i), "Expression"] <- rel.q.mean.log[which(rel.q.mean.log$Genes == i), "Expression"] - norm.factor.log[i]
    }
  }
  rel.q.norm.results <- list("rel.q.detailed" = rel.q.detailed, "rel.q.df" = rel.q.df, "rel.q.mean" = rel.q.mean, "rel.q.detailed.log" = rel.q.detailed.log, "rel.q.log" = rel.q.log, "rel.q.mean.log" = rel.q.mean.log, "ref.sample" = ref.sample)
  return(rel.q.norm.results)
}


## Calculate confidence intervals
rd.confint <- function(rel.q.mean, rel.q.mean.log, p) {
  for (i in 1:nrow(rel.q.mean)) {
    rel.q.mean.log$left.CI[i] <- rel.q.mean.log$Expression[i] - qt(1-p/2,rel.q.mean$Replicates[i]-1)*(rel.q.mean.log$SD[i]/rel.q.mean$Replicates[i])
    rel.q.mean.log$right.CI[i] <- rel.q.mean.log$Expression[i] + qt(1-p/2,rel.q.mean$Replicates[i]-1)*(rel.q.mean.log$SD[i]/rel.q.mean$Replicates[i])
  }
  rel.q.mean$left.CI <- 2^rel.q.mean.log$left.CI
  rel.q.mean$right.CI <- 2^rel.q.mean.log$right.CI
  rel.q.confint <- list("rel.q.mean" = rel.q.mean, "rel.q.mean.log" = rel.q.mean.log)
  return(rel.q.confint)
}


## Statistical tests
rd.statistics <- function(rel.q.df, rel.q.log, rel.q.mean, rel.q.mean.log, statistics, test.type, posthoc, ref.sample, p, sp.f) {
  frw <- 0
  few.repl.warn <- c()
  if (statistics == TRUE && min(table(rel.q.df$Samples))/length(unique(rel.q.df$Genes)) < 3) {
    statistics <- FALSE
    few.repl.warn <- "At least one of the groups contains too few replicates and it is not possible to perform statistical tests.\n"
    frw <- 1
  }
  res.posthoc <- t.pairs <- res.all.to.one <- res.pairs <- t.t <- res.KW <- res.aov <- list()
  if (statistics == TRUE) {
    if (posthoc == "all to one" && ref.sample %in% rel.q.mean$Samples) {
      rel.q.log$Samples <- relevel(rel.q.log$Samples, ref = as.character(ref.sample))
      # rel.q.log <- rel.q.log[order(rel.q.log$Genes, rel.q.log$Samples),]
      for (i in unique(rel.q.log$Genes)) {
        if (test.type == "parametric") {
          res.aov[[i]] <- aov(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
          if (length(unique(as.character(rel.q.log$Samples))) > 2) {
            if (car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05 && summary(res.aov[[i]])[[1]]$"Pr(>F)"[1] <= p) {
              set.seed(40075017)
              res.all.to.one[[i]] <- PMCMRplus::dunnettTest(res.aov[[i]])
            } else {
              if (car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] <= 0.05 && oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value <= p) {
                set.seed(40075017)
                res.all.to.one[[i]] <- PMCMRplus::tamhaneDunnettTest(res.aov[[i]])
              }
            }
            res.posthoc[[i]] <- as.data.frame(res.all.to.one[[i]]$p.value, stringsAsFactors = FALSE)
          } else {
            if (car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05) {
              t.t[[i]] <- t.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = TRUE)
            } else {
              t.t[[i]] <- t.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)
            }
            res.posthoc[[i]] <- as.data.frame(t.t[[i]]$p.value, stringsAsFactors = FALSE)
            colnames(res.posthoc[[i]]) <- ref.sample
            rownames(res.posthoc[[i]]) <- gsub("mean in group ", "", names(t.t[[i]]$estimate[2]))
          }
        }
        if (test.type == "non-parametric") {
          if (length(unique(as.character(rel.q.log$Samples))) > 2) {
            set.seed(40075017)
            res.KW[[i]] <-  PMCMRplus::kruskalTest(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            if (res.KW[[i]]$p.value <= p) {
              res.all.to.one[[i]] <- PMCMRplus::kwManyOneDunnTest(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            }
            res.posthoc[[i]] <- as.data.frame(res.all.to.one[[i]]$p.value, stringsAsFactors = FALSE)
          } else {
            t.t[[i]] <- wilcox.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            res.posthoc[[i]] <- as.data.frame(t.t[[i]]$p.value, stringsAsFactors = FALSE)
            colnames(res.posthoc[[i]]) <- ref.sample
            rownames(res.posthoc[[i]]) <- levels(rel.q.log$Samples)[2]
          }
        }
        res.posthoc[[i]]$Genes <- i
        res.posthoc[[i]]$Comparisons <- rownames(res.posthoc[[i]])
        # for (j in 1:nrow(res.posthoc[[i]])) {
        #   a <- which(rel.q.log$Genes == i & rel.q.log$Samples == res.posthoc[[i]][j, "Comparisons"])
        #   res.posthoc[[i]][j, "y"] <- max(rel.q.log[a, "Rel.quant"])
      }
      res.posthoc <- do.call("rbind", res.posthoc)
      res.posthoc$GeneSampleCombo <- paste0(res.posthoc$Genes, " in ", res.posthoc$Comparisons)
      for (k in unique(rel.q.log$GeneSampleCombo)) {
        if (length(which(res.posthoc$GeneSampleCombo == k)) > 0) {
          a <- which(rel.q.log$GeneSampleCombo == k)
          b <- rel.q.log[a, "Rel.quant"]
          rel.q.log[a[which(b == max(b))], "p.value"] <- res.posthoc[which(res.posthoc$GeneSampleCombo == k),as.character(ref.sample)]
        }
      }
      for (l in 1:nrow(rel.q.mean.log)) {
        s.name <- paste0(rel.q.mean.log[l, "Genes"], " in ", rel.q.mean.log[l, "Samples"])
        if (length(which(res.posthoc$GeneSampleCombo == s.name)) > 0) {
          rel.q.mean.log$p.value[l] <- res.posthoc[which(res.posthoc$GeneSampleCombo == s.name),as.character(ref.sample)]
        }
        else rel.q.mean.log$p.value[l] <- ""
      }
      rel.q.log$p.value <- get.rid.zeros(rel.q.log$p.value)
      rel.q.log$p.val.exp <- parse(text = gsub("e", "%*%10^", signif(rel.q.log$p.value, digits = 2)))
      rel.q.log$p.val.exp <- gsub("NA", "", rel.q.log$p.val.exp)
      rel.q.log$Samples <- factor(rel.q.log$Samples, levels = unique(rel.q.log$Samples))
      rel.q.log$asterisks <- as.character(symnum(rel.q.log$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", " ")))
      rel.q.mean.log$p.value <- as.numeric(rel.q.mean.log$p.value)
      rel.q.mean.log$p.value <- get.rid.zeros(rel.q.mean.log$p.value)
      rel.q.mean.log$p.val.exp <- parse(text = gsub("e", "%*%10^", signif(rel.q.mean.log$p.value, digits = 2)))
      rel.q.mean.log$p.val.exp <- gsub("NA", "", rel.q.mean.log$p.val.exp)
      rel.q.mean.log$asterisks <- as.character(symnum(rel.q.mean.log$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", " ")))
      rel.q.df <- cbind(rel.q.df, rel.q.log[,c("p.value", "p.val.exp", "asterisks")])
      rel.q.mean <- cbind(rel.q.mean, rel.q.mean.log[,c("p.value", "p.val.exp", "asterisks")])
      # colnames(rel.q.mean.log) <- c("Genes", "Samples", "Expression", "SD", "p.value", "p.val.exp", "asterisks")
    }
    for (k in 1:max(rel.q.df$Pairs)) {
      t.pairs[[k]] <- unique(subset(rel.q.df, Pairs == k)$Samples)
    }
    if (posthoc == "all pairs" | posthoc == "selected pairs") {
      for (i in unique(rel.q.log$Genes)) {
        if (posthoc == "all pairs") {
          # res.pairs[[i]] <- ""
          if (test.type == "parametric") {
            res.aov[[i]] <- aov(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            if (car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05 && summary(res.aov[[i]])[[1]]$"Pr(>F)"[1] <= p) {
              set.seed(40075017)
              res.pairs[[i]] <- PMCMRplus::tukeyTest(res.aov[[i]])
            } else {
              if (car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] <= 0.05 && oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value <= p) {
                if (length(res.aov[[1]]$effects)/length(res.aov[[i]]$coefficients) < 50) {
                  res.pairs[[i]] <- PMCMRplus::dunnettT3Test(res.aov[[i]])
                } else {
                  res.pairs[[i]] <- PMCMRplus::gamesHowellTest(res.aov[[i]])
                }
              }
            }
          }
          if (test.type == "non-parametric") {
            res.KW[[i]] <-  PMCMRplus::kruskalTest(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            if (res.KW[[i]]$p.value <= p) {
              res.pairs[[i]] <- PMCMRplus::kwAllPairsDunnTest(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            }
          }
        } else {
          if (test.type == "parametric") {
            res.aov[[i]] <- aov(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            if (car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05) {
              res.pairs[[i]] <- pairwise.t.test(subset(rel.q.log, Genes == i)$Rel.quant, subset(rel.q.log, Genes == i)$Samples, p.adjust.method = "none", paired = FALSE)
            } else {res.pairs[[i]] <- pairwise.t.test(subset(rel.q.log, Genes == i)$Rel.quant, subset(rel.q.log, Genes == i)$Samples, p.adjust.method = "none", pool.sd = FALSE, paired = FALSE)
            }
          }
          if (test.type == "non-parametric") {
            # res.KW[[i]] <-  PMCMRplus::kruskalTest(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i))
            # if (res.KW[[i]]$p.value <= p) {
              res.pairs[[i]] <- pairwise.wilcox.test(subset(rel.q.log, Genes == i)$Rel.quant, subset(rel.q.log, Genes == i)$Samples, p.adjust.method = "none", paired = FALSE)
            # }
          }
          for (k in colnames(res.pairs[[i]]$p.value)) {
            for (l in rownames(res.pairs[[i]]$p.value)) {
              w <- 0
              for (m in t.pairs) {
                ifelse (k %in% m & l %in% m == TRUE, v <- 1, v <- 0)
                w <- w + v
              }
              if (w == 0) {
                res.pairs[[i]]$p.value[l,k] <- NA
              }
            }
          }
          coln <- colnames(res.pairs[[i]]$p.value)
          rown <- rownames(res.pairs[[i]]$p.value)
          dimn <- dim(res.pairs[[i]]$p.value)
          res.pairs[[i]]$p.value <- matrix(p.adjust(res.pairs[[i]]$p.value, method = "bonferroni"), nrow = dimn[1], ncol = dimn[2])
          colnames(res.pairs[[i]]$p.value) <- coln
          rownames(res.pairs[[i]]$p.value) <- rown
        }
        if (rlang::is_empty(res.pairs[[i]]) == TRUE) {
          res.posthoc[[i]] <- NULL
        } else {
          res.posthoc[[i]] <- as.data.frame(res.pairs[[i]]$p.value, stringsAsFactors = FALSE)
          res.posthoc[[i]] <- data.frame(Genes = i, Group1=colnames(res.posthoc[[i]])[col(res.posthoc[[i]])[lower.tri(res.posthoc[[i]], diag = TRUE)]],
                                         Group2=rownames(res.posthoc[[i]])[row(res.posthoc[[i]])[lower.tri(res.posthoc[[i]], diag = TRUE)]],
                                         p.value=res.posthoc[[i]][lower.tri(res.posthoc[[i]], diag = TRUE)], stringsAsFactors = FALSE)
          # res.posthoc[[i]]$Samples <- res.posthoc[[i]]$Rel.quant <- res.posthoc[[i]]$y <- ""
          res.posthoc[[i]]$y <- ""
          for (j in 1:nrow(res.posthoc[[i]])) {
            a <- which(rel.q.log$Genes == i & rel.q.log$Samples == res.posthoc[[i]][j, "Group1"] | rel.q.log$Genes == i & rel.q.log$Samples == res.posthoc[[i]][j, "Group2"])
            b <- which(rel.q.mean.log$Genes == i & rel.q.mean.log$Samples == res.posthoc[[i]][j, "Group1"] | rel.q.mean.log$Genes == i & rel.q.mean.log$Samples == res.posthoc[[i]][j, "Group2"])
            res.posthoc[[i]][j, "y"] <- max(rel.q.log[seq(a[1],a[length(a)]), "Rel.quant"])
            res.posthoc[[i]][j, "y.s"] <- max(rel.q.mean.log[seq(b[1],b[length(b)]), "Expression"] + rel.q.mean.log[seq(b[1],b[length(b)]), "SD"])
            res.posthoc[[i]][j, "y.s.lin"] <- 2^max(rel.q.mean.log[seq(b[1],b[length(b)]), "right.CI"])
          }
        }
        # colnames(rel.q.mean.log) <- c("Genes", "Samples", "Expression", "SD")
      }
      # res.all.to.one.df[[i]]$Genes <- i
      # res.all.to.one.df[[i]]$Comparisons <- rownames(res.all.to.one.df[[i]])
      # res.Tukey[[i]] <- as.data.frame(TukeyHSD(res.aov[[i]], "Samples")[[1]])
      # res.Tukey[[i]]$Genes <- i
      # res.Tukey[[i]]$Comparisons <- rownames(res.Tukey[[i]])
      res.posthoc <- do.call("rbind", res.posthoc)
      res.posthoc$y <- as.numeric(res.posthoc$y)
      res.posthoc <- subset(res.posthoc, p.value <= p)
      res.posthoc$p.value <- get.rid.zeros(res.posthoc$p.value)
      res.posthoc <- res.posthoc[order(res.posthoc$Genes, res.posthoc$y),]
      res.posthoc$y.lin <- 2^res.posthoc$y
      # res.posthoc$y.s.lin <- 2^res.posthoc$y.s
      # The next loop is just a complicated logic to ensure that statistical significance bars on figures do not overlap and are adequately positioned regarding data and plot margins
      un.diff.c.lin <- un.diff.c <- un.diff.b.lin <- un.diff.b <- un.diff.lin <- un.diff <- span.lin <- span <- mx.lin <- mx <- mn.lin <- mn <- vector()
      for (i in unique(res.posthoc$Genes)) {
        mn[i] <- min(rel.q.log[which(rel.q.log$Genes == i), "Rel.quant"])
        mx[i] <- max(rel.q.log[which(rel.q.log$Genes == i), "Rel.quant"])
        span[i] <- mx[i] - mn[i]
        mn.lin[i] <- min(rel.q.df[which(rel.q.df$Genes == i), "Rel.quant"])
        mx.lin[i] <- max(rel.q.df[which(rel.q.df$Genes == i), "Rel.quant"])
        span.lin[i] <- mx.lin[i] - mn.lin[i]
        res.posthoc[which(res.posthoc$Genes == i), "y"] <- res.posthoc[which(res.posthoc$Genes == i), "y"] + 0.1*span[i]
        res.posthoc[which(res.posthoc$Genes == i), "y.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.lin"] + 0.1*span.lin[i]
        res.posthoc[which(res.posthoc$Genes == i), "y.s"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s"] + 0.1*span[i]
        res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"] + 0.1*span.lin[i]
        if (length(which(res.posthoc$Genes == i)) > 1 && posthoc != "selected pairs") {
          res.posthoc[which(res.posthoc$Genes == i), "y1"] <- seq(min(res.posthoc[which(res.posthoc$Genes == i), "y"]), ceiling(mx[i] + 0.1*span[i]), length.out = nrow(res.posthoc[which(res.posthoc$Genes == i),]))
          res.posthoc[which(res.posthoc$Genes == i), "y1.lin"] <- seq(min(res.posthoc[which(res.posthoc$Genes == i), "y.lin"]), ceiling(mx.lin[i] + 0.1*span.lin[i]), length.out = nrow(res.posthoc[which(res.posthoc$Genes == i),]))
          un.diff[i] <- sp.f*(res.posthoc[which(res.posthoc$Genes == i), "y1"][2] - res.posthoc[which(res.posthoc$Genes == i), "y1"][1])
          un.diff.lin[i] <- sp.f*(res.posthoc[which(res.posthoc$Genes == i), "y1.lin"][2] - res.posthoc[which(res.posthoc$Genes == i), "y1.lin"][1])
          res.posthoc[which(res.posthoc$Genes == i), "y2"] <- res.posthoc[which(res.posthoc$Genes == i), "y"]
          res.posthoc[which(res.posthoc$Genes == i), "y2.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.lin"]
          for (j in 2:nrow(res.posthoc[which(res.posthoc$Genes == i),])) {
            if (res.posthoc[which(res.posthoc$Genes == i), "y2"][j] <= res.posthoc[which(res.posthoc$Genes == i), "y2"][j-1] + un.diff[i]) {
              res.posthoc[which(res.posthoc$Genes == i), "y2"][j] <- res.posthoc[which(res.posthoc$Genes == i), "y2"][j-1] + un.diff[i]
            }
            if (res.posthoc[which(res.posthoc$Genes == i), "y2.lin"][j] <= res.posthoc[which(res.posthoc$Genes == i), "y2.lin"][j-1] + un.diff.lin[i]) {
              res.posthoc[which(res.posthoc$Genes == i), "y2.lin"][j] <- res.posthoc[which(res.posthoc$Genes == i), "y2.lin"][j-1] + un.diff.lin[i]
            }
          }
          res.posthoc[which(res.posthoc$Genes == i), "y3"] <- seq(min(res.posthoc[which(res.posthoc$Genes == i), "y.s"]), ceiling(mx[i] + 0.1*span[i]), length.out = nrow(res.posthoc[which(res.posthoc$Genes == i),]))
          res.posthoc[which(res.posthoc$Genes == i), "y3.lin"] <- seq(min(res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"]), ceiling(mx.lin[i] + 0.1*span.lin[i]), length.out = nrow(res.posthoc[which(res.posthoc$Genes == i),]))
          un.diff.b[i] <- sp.f*(res.posthoc[which(res.posthoc$Genes == i), "y3"][2] - res.posthoc[which(res.posthoc$Genes == i), "y3"][1])
          un.diff.b.lin[i] <- sp.f*(res.posthoc[which(res.posthoc$Genes == i), "y3.lin"][2] - res.posthoc[which(res.posthoc$Genes == i), "y3.lin"][1])
          res.posthoc[which(res.posthoc$Genes == i), "y4"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s"]
          res.posthoc[which(res.posthoc$Genes == i), "y4.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"]
          for (k in 2:nrow(res.posthoc[which(res.posthoc$Genes == i),])) {
            res.posthoc[which(res.posthoc$Genes == i), "y4"][k] <- max(res.posthoc[which(res.posthoc$Genes == i), "y4"][k], res.posthoc[which(res.posthoc$Genes == i), "y4"][k-1] + un.diff.b[i])
            res.posthoc[which(res.posthoc$Genes == i), "y4.lin"][k] <- max(res.posthoc[which(res.posthoc$Genes == i), "y4.lin"][k], res.posthoc[which(res.posthoc$Genes == i), "y4.lin"][k-1] + un.diff.b.lin[i])
          }
          res.posthoc[which(res.posthoc$Genes == i), "y5"] <- seq(max(0, min(res.posthoc[which(res.posthoc$Genes == i), "y.s"])), ceiling(max(0.5, mx[i]) + 0.1*span[i]), length.out = nrow(res.posthoc[which(res.posthoc$Genes == i),]))
          res.posthoc[which(res.posthoc$Genes == i), "y5.lin"] <- seq(max(0, min(res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"])), ceiling(max(0.5, mx.lin[i]) + 0.1*span.lin[i]), length.out = nrow(res.posthoc[which(res.posthoc$Genes == i),]))
          un.diff.c[i] <- sp.f*(res.posthoc[which(res.posthoc$Genes == i), "y5"][2] - res.posthoc[which(res.posthoc$Genes == i), "y5"][1])
          un.diff.c.lin[i] <- sp.f*(res.posthoc[which(res.posthoc$Genes == i), "y5.lin"][2] - res.posthoc[which(res.posthoc$Genes == i), "y5.lin"][1])
          res.posthoc[which(res.posthoc$Genes == i), "y6"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s"]
          res.posthoc[which(res.posthoc$Genes == i), "y6.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"]
          res.posthoc[which(res.posthoc$Genes == i), "y6"][1] <- max(0.5, res.posthoc[which(res.posthoc$Genes == i), "y6"][1])
          # res.posthoc[which(res.posthoc$Genes == i), "y6.lin"][1] <- max(0.5, res.posthoc[which(res.posthoc$Genes == i), "y6.lin"][1])
          for (l in 2:nrow(res.posthoc[which(res.posthoc$Genes == i),])) {
            res.posthoc[which(res.posthoc$Genes == i), "y6"][l] <- max(res.posthoc[which(res.posthoc$Genes == i), "y6"][l], res.posthoc[which(res.posthoc$Genes == i), "y6"][l-1] + un.diff.c[i])
            res.posthoc[which(res.posthoc$Genes == i), "y6.lin"][l] <- max(res.posthoc[which(res.posthoc$Genes == i), "y6.lin"][l], res.posthoc[which(res.posthoc$Genes == i), "y6.lin"][l-1] + un.diff.c.lin[i])
          }
        } else {
          res.posthoc[which(res.posthoc$Genes == i), "y2"] <- res.posthoc[which(res.posthoc$Genes == i), "y1"] <- res.posthoc[which(res.posthoc$Genes == i), "y"]
          res.posthoc[which(res.posthoc$Genes == i), "y2.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y1.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.lin"]
          res.posthoc[which(res.posthoc$Genes == i), "y5"] <- res.posthoc[which(res.posthoc$Genes == i), "y4"] <- res.posthoc[which(res.posthoc$Genes == i), "y3"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s"]
          res.posthoc[which(res.posthoc$Genes == i), "y6.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y5.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y4.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y3.lin"] <- res.posthoc[which(res.posthoc$Genes == i), "y.s.lin"]
          res.posthoc[which(res.posthoc$Genes == i), "y6"] <- pmax(0.5, res.posthoc[which(res.posthoc$Genes == i), "y5"])
          # res.posthoc[which(res.posthoc$Genes == i), "y6.lin"] <- pmax(0.5, res.posthoc[which(res.posthoc$Genes == i), "y5.lin"])
        }
      }
      res.posthoc$p.val.exp <- parse(text = gsub("e", "%*%10^", signif(res.posthoc$p.value, digits = 2)))
      res.posthoc$p.val.exp <- gsub("NA", "", res.posthoc$p.val.exp)
      res.posthoc$asterisks <- as.character(symnum(res.posthoc$p.value, corr = FALSE, na = FALSE, cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 1), symbols = c("****", "***", "**", "*", ".", " ")))
    }
  }
  #   statistics.results <- list("rel.q.df" = rel.q.df, "rel.q.mean" = rel.q.mean, "rel.q.log" = rel.q.log, "rel.q.mean.log" = rel.q.mean.log, "res.posthoc" = res.posthoc, "res.all.to.one" = res.all.to.one, "res.pairs" = res.pairs, "t.t" = t.t, "res.KW" = res.KW, "res.aov" = res.aov, "frw" = frw, "few.repl.warn" = few.repl.warn)
  #   return(statistics.results)
  # }

  stat.test <- list()
  if (statistics == TRUE) {
    if (posthoc == "all to one" && ref.sample %in% rel.q.mean$Samples) {
      for (i in unique(rel.q.log$Genes)) {
        if (length(unique(as.character(rel.q.log$Samples))) > 2) {
          if (test.type == "parametric" && car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05 && summary(res.aov[[i]])[[1]]$"Pr(>F)"[1] <= p) {
            # stat.test[[i]] <- signif(summary(res.aov[[i]])[[1]]$"Pr(>F)"[1], digits = 2)
            stat.test[[i]] <- paste0("Statistical test(s): one-way ANOVA (p = ", signif(summary(res.aov[[i]])[[1]]$"Pr(>F)"[1], digits = 2),"), ", gsub("\t\n *", " ", res.all.to.one[[i]]$method))
          }
          if (test.type == "parametric" && car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] <= 0.05 && oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value <= p) {
            stat.test[[i]] <- paste0("Statistical test(s): Welch's ANOVA (p = ", signif(oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value, digits = 2),"), ", gsub("\t\n *", " ", res.all.to.one[[i]]$method))
            # stat.test[[i]] <- signif(summary(res.aov[[i]])[[1]]$"Pr(>F)"[1], digits = 2)
          }
          if (test.type == "non-parametric") {
            stat.test[[i]] <- paste0("Statistical test(s): Kruskal-Wallis (p = ", signif(res.KW[[i]]$p.value, digits = 2),"), ", gsub("\t\n *", " ", res.all.to.one[[i]]$method))
          }
        } else {
          stat.test[[i]] <- paste0("Statistical test: ", t.t[[i]]$method)
        }
      }
    }
    if (posthoc == "all pairs") {
      for (i in unique(rel.q.log$Genes)) {
        if (test.type == "parametric" && car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05 && summary(res.aov[[i]])[[1]]$"Pr(>F)"[1] <= p) {
          stat.test[[i]] <- paste0("Statistical test(s): one-way ANOVA (p = ", signif(summary(res.aov[[i]])[[1]]$"Pr(>F)"[1], digits = 2),"), ", gsub("\n\t\t", " ", res.pairs[[i]]$method))
        }
        if (test.type == "parametric" && car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] <= 0.05 && oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value <= p) {
          stat.test[[i]] <- paste0("Statistical test(s): Welch's ANOVA (p = ", signif(oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value, digits = 2),"), ", gsub("\n\t\t", " ", res.pairs[[i]]$method))
        }
        if (test.type == "non-parametric") {
          stat.test[[i]] <- paste0("Statistical test(s): Kruskal-Wallis (p = ", signif(res.KW[[i]]$p.value, digits = 2),"), ", gsub("\n\t\t", " ", res.pairs[[i]]$method))
        }
      }
    }
    if (posthoc == "selected pairs") {
      for (i in unique(rel.q.log$Genes)) {
        if (test.type == "parametric" && car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] > 0.05 && summary(res.aov[[i]])[[1]]$"Pr(>F)"[1] <= p) {
          stat.test[[i]] <- paste0("Statistical test(s): one-way ANOVA (p = ", signif(summary(res.aov[[i]])[[1]]$"Pr(>F)"[1], digits = 2),"), pairwise ", gsub("\n\t\t", " ", res.pairs[[i]]$method))
        }
        if (test.type == "parametric" && car::leveneTest(res.aov[[i]])$"Pr(>F)"[1] <= 0.05 && oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value <= p) {
          stat.test[[i]] <- paste0("Statistical test(s): Welch's ANOVA (p = ", signif(oneway.test(Rel.quant ~ Samples, data = subset(rel.q.log, Genes == i), var.equal = FALSE)$p.value, digits = 2),"), pairwise ", gsub("\n\t\t", " ", res.pairs[[i]]$method))
        }
        if (test.type == "non-parametric") {
          stat.test[[i]] <- paste0("Statistical test(s): pairwise ", gsub("\n\t\t", " ", res.pairs[[i]]$method))
        }
      }
    }
  }
  statistics.results <- list("rel.q.df" = rel.q.df, "rel.q.mean" = rel.q.mean, "rel.q.log" = rel.q.log, "rel.q.mean.log" = rel.q.mean.log, "res.posthoc" = res.posthoc, "ref.sample" = ref.sample, "statistics" = statistics, "stat.test" = stat.test, "frw" = frw, "few.repl.warn" = few.repl.warn)
  return(statistics.results)
}


## Plot relative expression
rd.plot.p1 <- function(rel.q.df, rel.q.mean, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
  p1 <- list()
  if (statistics == FALSE | posthoc == "all to one") {
    for (i in GOIs) {
      p1[[i]] <-ggplot2::ggplot(subset(rel.q.df, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples, fill=Samples)) +
        ggbeeswarm::geom_beeswarm(cex = 2, size = 2, groupOnX = TRUE) +
        ggplot2::ylim(0, 1.2*max(subset(rel.q.df, Genes == i)$Rel.quant, na.rm = TRUE)) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "values") {
      for (i in unique((rel.q.df %>% filter(p.value <= p))$Genes)) {
        p1[[i]] <- p1[[i]] + ggplot2::geom_text(data = rel.q.df %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.df, Genes == i)$Rel.quant) - min(subset(rel.q.df, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
      }
      for (i in unique((rel.q.df %>% filter(p.value <= .Machine$double.eps))$Genes)) {
        p1[[i]] <- p1[[i]] + ggplot2::geom_text(data = rel.q.df %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.df, Genes == i)$Rel.quant) - min(subset(rel.q.df, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
      }
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "asterisks") {
      for (i in unique((rel.q.df %>% filter(p.value <= p))$Genes)) {
        p1[[i]] <- p1[[i]] + ggplot2::geom_text(data = rel.q.df %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = 0.1*(max(subset(rel.q.df, Genes == i)$Rel.quant) - min(subset(rel.q.df, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3 + 1)
      }
    }
    for (i in GOIs) {
      p1[[i]] <- p1[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  } else {
    for (i in GOIs) {
      p1[[i]] <-ggplot2::ggplot(subset(rel.q.df, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples, fill=Samples)) +
        ggbeeswarm::geom_beeswarm(cex = 2, size = 2, groupOnX = TRUE) +
        ggplot2::ylim(0, max(ceiling(max(subset(rel.q.df, Genes == i)$Rel.quant, na.rm = TRUE)*1.025), max(subset(res.posthoc, Genes == i)$y2.lin)*1.025)) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p1[[i]] <- p1[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y2.lin), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
      }
    }
    if (statistics == TRUE && sign.repr == "asterisks"  && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p1[[i]] <- p1[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y2.lin), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
      }
    }
    for (i in GOIs) {
      p1[[i]] <- p1[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  }
  ml <- gridExtra::marrangeGrob(grobs = p1, nrow = 3, ncol = 1)
  p1.results <- list("p1" = p1, "ml" = ml)
  return(p1.results)
}

rd.plot.p2 <- function(rel.q.mean, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
  p2 <- list()
  if (statistics == FALSE | posthoc == "all to one") {
    for (i in GOIs) {
      p2[[i]] <- ggplot2::ggplot(subset(rel.q.mean, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
        ggplot2::geom_point(size=2) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin=left.CI, ymax=right.CI), width=.2) +
        ggplot2::ylim(0, 1.2*max(subset(rel.q.mean, Genes == i)$right.CI)) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "values") {
      for (i in unique((rel.q.mean %>% filter(p.value <= p))$Genes)) {
        p2[[i]] <- p2[[i]] + ggplot2::geom_text(data = rel.q.mean %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = subset(rel.q.mean, Genes == i & p.value <= p & p.value > .Machine$double.eps)$right.CI - subset(rel.q.mean, Genes == i & p.value <= p & p.value > .Machine$double.eps)$Expression + 0.1*max(subset(rel.q.mean, Genes == i)$right.CI), colour = "gray20", size = font.size/3)
      }
      for (i in unique((rel.q.mean %>% filter(p.value <= .Machine$double.eps))$Genes)) {
        p2[[i]] <- p2[[i]] + ggplot2::geom_text(data = rel.q.mean %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = subset(rel.q.mean, Genes == i & p.value <= .Machine$double.eps)$right.CI - subset(rel.q.mean, Genes == i & p.value <= .Machine$double.eps)$Expression + 0.1*max(subset(rel.q.mean, Genes == i)$right.CI), colour = "gray20", size = font.size/3)
      }
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "asterisks") {
      for (i in unique((rel.q.mean %>% filter(p.value <= p))$Genes)) {
        p2[[i]] <- p2[[i]] + ggplot2::geom_text(data = rel.q.mean %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = subset(rel.q.mean, Genes == i & p.value <= p)$right.CI - subset(rel.q.mean, Genes == i & p.value <= p)$Expression + 0.1*max(subset(rel.q.mean, Genes == i)$right.CI), colour = "gray20", size = font.size/3 + 1)
      }
    }
    for (i in GOIs) {
      p2[[i]] <- p2[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  } else {
    for (i in GOIs) {
      p2[[i]] <- ggplot2::ggplot(subset(rel.q.mean, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
        ggplot2::geom_point(size=2) +
        ggplot2::geom_errorbar(aes(ymin=left.CI, ymax=right.CI), width=.2) +
        ggplot2::ylim(0, max(ceiling(max(subset(rel.q.mean, Genes == i)$right.CI)*1.025), max(subset(res.posthoc, Genes == i)$y4.lin)*1.025)) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p2[[i]] <- p2[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y4.lin), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
      }
    }
    if (statistics == TRUE && sign.repr == "asterisks" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p2[[i]] <- p2[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y4.lin), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
      }
    }
    for (i in GOIs) {
      p2[[i]] <- p2[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  }
  ml <- gridExtra::marrangeGrob(grobs = p2, nrow = 3, ncol = 1)
  p2.results <- list("p2" = p2, "ml" = ml)
  return(p2.results)
}

rd.plot.p3 <- function(rel.q.mean, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
  p3 <- list()
  if (statistics == FALSE | posthoc == "all to one") {
    for (i in GOIs) {
      p3[[i]] <-ggplot2::ggplot(subset(rel.q.mean, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::geom_errorbar(aes(ymin=left.CI, ymax=right.CI), width=.2, position=position_dodge(.9), color="gray20") +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "values") {
      for (i in unique((rel.q.mean %>% filter(p.value <= p))$Genes)) {
        p3[[i]] <- p3[[i]] + ggplot2::geom_text(data = rel.q.mean %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = subset(rel.q.mean, Genes == i & p.value <= p & p.value > .Machine$double.eps)$right.CI - subset(rel.q.mean, Genes == i & p.value <= p & p.value > .Machine$double.eps)$Expression + 0.05*max(subset(rel.q.mean, Genes == i)$right.CI), colour = "gray20", size = font.size/3)
      }
      for (i in unique((rel.q.mean %>% filter(p.value <= .Machine$double.eps))$Genes)) {
        p3[[i]] <- p3[[i]] + ggplot2::geom_text(data = rel.q.mean %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = subset(rel.q.mean, Genes == i & p.value <= .Machine$double.eps)$right.CI - subset(rel.q.mean, Genes == i & p.value <= .Machine$double.eps)$Expression + 0.05*max(subset(rel.q.mean, Genes == i)$right.CI), colour = "gray20", size = font.size/3)
      }
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "asterisks") {
      for (i in unique((rel.q.mean %>% filter(p.value <= p))$Genes)) {
        p3[[i]] <- p3[[i]] + ggplot2::geom_text(data = rel.q.mean %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = subset(rel.q.mean, Genes == i & p.value <= p)$right.CI - subset(rel.q.mean, Genes == i & p.value <= p)$Expression + 0.05*max(subset(rel.q.mean, Genes == i)$right.CI), colour = "gray20", size = font.size/3 + 1)
      }
    }
    for (i in GOIs) {
      p3[[i]] <- p3[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  } else {
    for (i in GOIs) {
      p3[[i]] <- ggplot2::ggplot(subset(rel.q.mean, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
        ggplot2::geom_col(position = "dodge") +
        ggplot2::geom_errorbar(aes(ymin=left.CI, ymax=right.CI), width=.2, position=position_dodge(.9), color="gray20") +
        ggplot2::ylim(0, max(subset(res.posthoc, Genes == i)$y6.lin)*1.025) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p3[[i]] <- p3[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y6.lin), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
      }
    }
    if (statistics == TRUE && sign.repr == "asterisks" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p3[[i]] <- p3[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y6.lin), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
      }
    }
    for (i in GOIs) {
      p3[[i]] <- p3[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  }
  ml <- gridExtra::marrangeGrob(grobs = p3, nrow = 3, ncol = 1)
  p3.results <- list("p3" = p3, "ml" = ml)
  return(p3.results)
}

rd.plot.p2n <- function(rel.q.df, rel.q.mean, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
  p2n <- list()
  if (statistics == FALSE | posthoc == "all to one") {
    for (i in GOIs) {
      p2n[[i]] <- ggplot2::ggplot(subset(rel.q.df, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples)) +
        ggplot2::geom_boxplot() +
        ggplot2::ylim(0, 1.2*max(subset(rel.q.df, Genes == i)$Rel.quant, na.rm = TRUE)) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "values") {
      for (i in unique((rel.q.df %>% filter(p.value <= p))$Genes)) {
        p2n[[i]] <- p2n[[i]] + ggplot2::geom_text(data = rel.q.df %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.df, Genes == i)$Rel.quant) - min(subset(rel.q.df, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
      }
      for (i in unique((rel.q.df %>% filter(p.value <= .Machine$double.eps))$Genes)) {
        p2n[[i]] <- p2n[[i]] + ggplot2::geom_text(data = rel.q.df %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.df, Genes == i)$Rel.quant) - min(subset(rel.q.df, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
      }
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "asterisks") {
      for (i in unique((rel.q.df %>% filter(p.value <= p))$Genes)) {
        p2n[[i]] <- p2n[[i]] + ggplot2::geom_text(data = rel.q.df %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = 0.1*(max(subset(rel.q.df, Genes == i)$Rel.quant) - min(subset(rel.q.df, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3 + 1)
      }
    }
    for (i in GOIs) {
      p2n[[i]] <- p2n[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  } else {
    for (i in GOIs) {
      p2n[[i]] <- ggplot2::ggplot(subset(rel.q.df, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples)) +
        ggplot2::geom_boxplot() +
        ggplot2::ylim(0, max(ceiling(max(subset(rel.q.df, Genes == i)$Rel.quant, na.rm = TRUE)*1.025), max(subset(res.posthoc, Genes == i)$y2.lin)*1.025)) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab("Relative expression") +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p2n[[i]] <- p2n[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y2.lin), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
      }
    }
    if (statistics == TRUE && sign.repr == "asterisks"  && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p2n[[i]] <- p2n[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y2.lin), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
      }
    }
    for (i in GOIs) {
      p2n[[i]] <- p2n[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  }
  ml <- gridExtra::marrangeGrob(grobs = p2n, nrow = 3, ncol = 1)
  p2n.results <- list("p2n" = p2n, "ml" = ml)
  return(p2n.results)
}


rd.plot.p4 <- function(rel.q.log, rel.q.mean, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
  p4 <- list()
  if (statistics == FALSE | posthoc == "all to one") {
    for (i in GOIs) {
      p4[[i]] <-ggplot2::ggplot(subset(rel.q.log, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples, fill=Samples)) +
        ggbeeswarm::geom_beeswarm(cex = 2, size = 2, groupOnX = TRUE) +
        ggplot2::ylim(floor(min(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)), ceiling(max(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE) + 0.025*(max(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE) - min(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "values") {
      for (i in unique((rel.q.log %>% filter(p.value <= p))$Genes)) {
        p4[[i]] <- p4[[i]] + ggplot2::geom_text(data = rel.q.log %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
      }
      for (i in unique((rel.q.log %>% filter(p.value <= .Machine$double.eps))$Genes)) {
        p4[[i]] <- p4[[i]] + ggplot2::geom_text(data = rel.q.log %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
      }
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "asterisks") {
      for (i in unique((rel.q.log %>% filter(p.value <= p))$Genes)) {
        p4[[i]] <- p4[[i]] + ggplot2::geom_text(data = rel.q.log %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = 0.1*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3 + 1)
      }
    }
    for (i in GOIs) {
      p4[[i]] <- p4[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  } else {
    for (i in GOIs) {
      p4[[i]] <-ggplot2::ggplot(subset(rel.q.log, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples, fill=Samples)) +
        ggbeeswarm::geom_beeswarm(cex = 2, size = 2, groupOnX = TRUE) +
        ggplot2::ylim(floor(min(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)), max(ceiling(1.025*max(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)), 1.025*max(subset(res.posthoc, Genes == i)$y2))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p4[[i]] <- p4[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y2), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
      }
    }
    if (statistics == TRUE && sign.repr == "asterisks"  && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p4[[i]] <- p4[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y2), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
      }
    }
    for (i in GOIs) {
      p4[[i]] <- p4[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  }
  ml <- gridExtra::marrangeGrob(grobs = p4, nrow = 3, ncol = 1)
  p4.results <- list("p4" = p4, "ml" = ml)
  return(p4.results)
}

rd.plot.p5 <- function(rel.q.mean.log, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
  p5 <- list()
  if (statistics == FALSE | posthoc == "all to one") {
    for (i in GOIs) {
      p5[[i]] <- ggplot2::ggplot(subset(rel.q.mean.log, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
        ggplot2::geom_point(size=2) +
        ggplot2::geom_errorbar(aes(ymin=Expression-SD, ymax=Expression+SD), width=.2) +
        ggplot2::ylim(floor(min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD)), ceiling(max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD)) + 0.025*(max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD) - min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean.log$Samples && sign.repr == "values") {
      for (i in unique((rel.q.mean.log %>% filter(p.value <= p))$Genes)) {
        p5[[i]] <- p5[[i]] + ggplot2::geom_text(data = rel.q.mean.log %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = subset(rel.q.mean.log, Genes == i & p.value <= p & p.value > .Machine$double.eps)$SD + 0.1*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression)), colour = "gray20", size = font.size/3)
      }
      for (i in unique((rel.q.mean.log %>% filter(p.value <= .Machine$double.eps))$Genes)) {
        p5[[i]] <- p5[[i]] + ggplot2::geom_text(data = rel.q.mean.log %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = subset(rel.q.mean.log, Genes == i & p.value <= .Machine$double.eps)$SD + 0.1*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression)), colour = "gray20", size = font.size/3)
      }
    }
    if (statistics == TRUE && ref.sample %in% rel.q.mean.log$Samples && sign.repr == "asterisks") {
      for (i in unique((rel.q.mean.log %>% filter(p.value <= p))$Genes)) {
        p5[[i]] <- p5[[i]] + ggplot2::geom_text(data = rel.q.mean.log %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = subset(rel.q.mean.log, Genes == i & p.value <= p)$SD + 0.1*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression)), colour = "gray20", size = font.size/3 + 1)
      }
    }
    for (i in GOIs) {
      p5[[i]] <- p5[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  } else {
    for (i in GOIs) {
      p5[[i]] <- ggplot2::ggplot(subset(rel.q.mean.log, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
        ggplot2::geom_point(size=2) +
        ggplot2::geom_errorbar(aes(ymin=Expression-SD, ymax=Expression+SD), width=.2) +
        ggplot2::ylim(floor(min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD)), max(ceiling(1.025*max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD)), 1.025*max(subset(res.posthoc, Genes == i)$y4))) +
        ggplot2::theme_bw() +
        ggplot2::ggtitle(i) +
        ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
        ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
        ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
        ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
        ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
        ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
    }
    if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p5[[i]] <- p5[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y4), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
      }
    }
    if (statistics == TRUE && sign.repr == "asterisks" && length(res.posthoc$Genes) != 0) {
      for (i in unique(res.posthoc$Genes)) {
        p5[[i]] <- p5[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y4), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
      }
    }
    for (i in GOIs) {
      p5[[i]] <- p5[[i]] + ggplot2::labs(caption = stat.test[[i]])
    }
  }
  ml <- gridExtra::marrangeGrob(grobs = p5, nrow = 3, ncol = 1)
  p5.results <- list("p5" = p5, "ml" = ml)
  return(p5.results)
}

rd.plot.p6 <- function(rel.q.mean.log, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
p6 <- list()
if (statistics == FALSE | posthoc == "all to one") {
  for (i in GOIs) {
  p6[[i]] <- ggplot2::ggplot(subset(rel.q.mean.log, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::geom_errorbar(aes(ymin=Expression-SD, ymax=Expression+SD), width=.2, position=position_dodge(.9), color="gray20") +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(i) +
    ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
    ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
    ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
    ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
    ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
    ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
    ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  }
  if (statistics == TRUE && ref.sample %in% rel.q.mean.log$Samples && sign.repr == "values") {
    for (i in unique((rel.q.mean.log %>% filter(p.value <= p))$Genes)) {
      p6[[i]] <- p6[[i]] + ggplot2::geom_text(data = rel.q.mean.log %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = ifelse(subset(rel.q.mean.log, Genes == i & p.value <= p)$Expression < 0, -subset(rel.q.mean.log, Genes == i & p.value <= p)$SD - 0.05*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression)), subset(rel.q.mean.log, Genes == i & p.value <= p)$SD + 0.1*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression))), colour = "gray20", size = font.size/3)
    }
    for (i in unique((rel.q.mean.log %>% filter(p.value <= .Machine$double.eps))$Genes)) {
      p6[[i]] <- p6[[i]] + ggplot2::geom_text(data = rel.q.mean.log %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = ifelse(subset(rel.q.mean.log, Genes == i & p.value <= p)$Expression < 0, -subset(rel.q.mean.log, Genes == i & p.value <= p)$SD - 0.05*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression)), subset(rel.q.mean.log, Genes == i & p.value <= p)$SD + 0.1*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression))), colour = "gray20", size = font.size/3)
    }
  }
  if (statistics == TRUE && ref.sample %in% rel.q.mean.log$Samples && sign.repr == "asterisks") {
    for (i in unique((rel.q.mean.log %>% filter(p.value <= p))$Genes)) {
      p6[[i]] <- p6[[i]] + ggplot2::geom_text(data = rel.q.mean.log %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = ifelse(subset(rel.q.mean.log, Genes == i & p.value <= p)$Expression < 0, -subset(rel.q.mean.log, Genes == i & p.value <= p)$SD - 0.05*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression)), subset(rel.q.mean.log, Genes == i & p.value <= p)$SD + 0.1*(max(subset(rel.q.mean.log, Genes == i)$Expression) - min(subset(rel.q.mean.log, Genes == i)$Expression))), colour = "gray20", size = font.size/3 + 1)
    }
  }
  for (i in GOIs) {
    p6[[i]] <- p6[[i]] + ggplot2::labs(caption = stat.test[[i]])
  }
} else {
  for (i in GOIs) {
    p6[[i]] <- ggplot2::ggplot(subset(rel.q.mean.log, Genes == i), ggplot2::aes(x=Samples, y=Expression, colour=Samples, fill=Samples)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::geom_errorbar(aes(ymin=Expression-SD, ymax=Expression+SD), width=.2, position=position_dodge(.9), color="gray20") +
      # ggplot2::ylim(floor(min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD)), max(ceiling(max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD) + 0.025*(max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD) - min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD))), max(subset(res.posthoc, Genes == i)$y6) + 0.025*(max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD) - min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD)))) +
      ggplot2::ylim(floor(min(subset(rel.q.mean.log, Genes == i)$Expression - subset(rel.q.mean.log, Genes == i)$SD)), max(ceiling(1.025*max(subset(rel.q.mean.log, Genes == i)$Expression + subset(rel.q.mean.log, Genes == i)$SD)), 1.025*max(subset(res.posthoc, Genes == i)$y6))) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(i) +
      ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  }
  if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
    for (i in unique(res.posthoc$Genes)) {
      p6[[i]] <- p6[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y6), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
    }
  }
  if (statistics == TRUE && sign.repr == "asterisks" && length(res.posthoc$Genes) != 0) {
    for (i in unique(res.posthoc$Genes)) {
      p6[[i]] <- p6[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y6), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
    }
  }
  for (i in GOIs) {
    p6[[i]] <- p6[[i]] + ggplot2::labs(caption = stat.test[[i]])
  }
}
ml <- gridExtra::marrangeGrob(grobs = p6, nrow = 3, ncol = 1)
p6.results <- list("p6" = p6, "ml" = ml)
return(p6.results)
}

rd.plot.p5n <- function(rel.q.log, rel.q.mean, res.posthoc, ref.sample, GOIs, statistics, posthoc, sign.repr, p, stat.test, font.size) {
p5n <- list()
if (statistics == FALSE | posthoc == "all to one") {
  for (i in GOIs) {
    p5n[[i]] <- ggplot2::ggplot(subset(rel.q.log, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples)) +
      ggplot2::geom_boxplot() +
      ggplot2::ylim(floor(min(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)), ceiling(max(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE) + 0.025*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)))) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(i) +
      ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  }
  if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "values") {
    for (i in unique((rel.q.log %>% filter(p.value <= p))$Genes)) {
      p5n[[i]] <- p5n[[i]] + ggplot2::geom_text(data = rel.q.log %>% filter(Genes == i & p.value <= p & p.value > .Machine$double.eps), ggplot2::aes(label = paste0("p == ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
    }
    for (i in unique((rel.q.log %>% filter(p.value <= .Machine$double.eps))$Genes)) {
      p5n[[i]] <- p5n[[i]] + ggplot2::geom_text(data = rel.q.log %>% filter(Genes == i & p.value <= .Machine$double.eps), ggplot2::aes(label = paste0("p <= ", p.val.exp)), parse = TRUE, nudge_y = 0.1*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3)
    }
  }
  if (statistics == TRUE && ref.sample %in% rel.q.mean$Samples && sign.repr == "asterisks") {
    for (i in unique((rel.q.log %>% filter(p.value <= p))$Genes)) {
      p5n[[i]] <- p5n[[i]] + ggplot2::geom_text(data = rel.q.log %>% filter(Genes == i & p.value <= p), ggplot2::aes(label = asterisks), nudge_y = 0.1*(max(subset(rel.q.log, Genes == i)$Rel.quant) - min(subset(rel.q.log, Genes == i)$Rel.quant)), colour = "gray20", size = font.size/3 + 1)
    }
  }
  for (i in GOIs) {
    p5n[[i]] <- p5n[[i]] + ggplot2::labs(caption = stat.test[[i]])
  }
} else {
  for (i in GOIs) {
    p5n[[i]] <- ggplot2::ggplot(subset(rel.q.log, Genes == i), ggplot2::aes(x=Samples, y=Rel.quant, colour=Samples)) +
      ggplot2::geom_boxplot() +
      ggplot2::ylim(floor(min(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)), max(ceiling(1.025*max(subset(rel.q.log, Genes == i)$Rel.quant, na.rm = TRUE)), 1.025*max(subset(res.posthoc, Genes == i)$y2))) +
      ggplot2::theme_bw() +
      ggplot2::ggtitle(i) +
      ggplot2::ylab(expression("log"[2]~"(relative expression)")) +
      ggplot2::theme(axis.text=ggplot2::element_text(size = font.size)) +
      ggplot2::theme(axis.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(plot.title=ggplot2::element_text(size = font.size + 4, face = "bold")) +
      ggplot2::theme(plot.caption=ggplot2::element_text(size = font.size - 1)) +
      ggplot2::theme(legend.title=ggplot2::element_text(size = font.size + 2)) +
      ggplot2::theme(legend.text=ggplot2::element_text(size = font.size))
  }
  if (statistics == TRUE && sign.repr == "values" && length(res.posthoc$Genes) != 0) {
    for (i in unique(res.posthoc$Genes)) {
      p5n[[i]] <- p5n[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = paste0("p == ", p.val.exp), y_position = y2), size = 0.25, textsize = font.size/3, colour = "gray20", vjust = 0.2, tip_length = 0.02, show.legend = FALSE, parse = TRUE, manual = TRUE)
    }
  }
  if (statistics == TRUE && sign.repr == "asterisks"  && length(res.posthoc$Genes) != 0) {
    for (i in unique(res.posthoc$Genes)) {
      p5n[[i]] <- p5n[[i]] + ggsignif::geom_signif(data = res.posthoc %>% filter(Genes == i & p.value <= p), inherit.aes = FALSE, ggplot2::aes(xmin = Group1, xmax = Group2, annotations = asterisks, y_position = y2), size = 0.25, textsize = font.size/3 + 1, colour = "gray20", vjust = 0.4, tip_length = 0.02, show.legend = FALSE, manual = TRUE)
    }
  }
  for (i in GOIs) {
    p5n[[i]] <- p5n[[i]] + ggplot2::labs(caption = stat.test[[i]])
  }
}
ml <- gridExtra::marrangeGrob(grobs = p5n, nrow = 3, ncol = 1)
p5n.results <- list("p5n" = p5n, "ml" = ml)
return(p5n.results)
}


## Save tables
# confint.left <- paste0(" (lower bound of the ", (1-p)*100, "% confidence interval")
# confint.right <- paste0(" (upper bound of the ", (1-p)*100, "% confidence interval")
# colnames(rel.q.detailed) <- gsub(".left.CI", confint.left, colnames(rel.q.detailed))
# colnames(rel.q.detailed) <- gsub(".right.CI", confint.right, colnames(rel.q.detailed))
rd.save.tables <- function(input.table, rel.q.detailed, rel.q.detailed.log, rel.q.mean, rel.q.mean.log, p) {
write.csv(rel.q.detailed, paste0(gsub(".csv", "", input.table), "_relative_expression_replicates.csv"))
colnames(rel.q.detailed.log) <- gsub("$", " (log2)", colnames(rel.q.detailed.log))
write.csv(rel.q.detailed.log, paste0(gsub(".csv", "", input.table), "_relative_log_expression_replicates.csv"))
colnames(rel.q.mean) <- c("Genes", "Samples", "Pairs", "Replicates", "Relative expression", paste0("Lower bound of the ", (1-p)*100, "% confidence interval"), paste0("Upper bound of the ", (1-p)*100, "% confidence interval"))
write.csv(rel.q.mean[,1:7], paste0(gsub(".csv", "", input.table), "_relative_expression_averaged.csv"), row.names = FALSE)
colnames(rel.q.mean.log) <- c("Genes", "Samples", "Pairs", "log2(relative expression)", "Standard deviation", paste0("Lower bound of the ", (1-p)*100, "% confidence interval"), paste0("Upper bound of the ", (1-p)*100, "% confidence interval"))
write.csv(rel.q.mean.log[,1:7], paste0(gsub(".csv", "", input.table), "_relative_log_expression_averaged.csv"), row.names = FALSE)
save.tables <- list("rel.q.detailed.log" = rel.q.detailed.log, "rel.q.mean" = rel.q.mean[,1:7], "rel.q.mean.log" = rel.q.mean.log[,1:7])
return(save.tables)
}


## Print warning messages if any
rd.warn <- function(ref.sample, rel.q.mean, noref.warn, statistics, posthoc, nostatref.warn, frw, few.repl.warn) {
  warnings <- list()
  if ((ref.sample %in% rel.q.mean$Samples) == FALSE) {
    warnings[["noref.warn"]] <- noref.warn
    if (statistics == TRUE && posthoc == "all to one") {
      warnings[["nostatref.warn"]] <- nostatref.warn
    }
  }
  if (frw == 1) {
    warnings[["few.repl.warn"]] <- few.repl.warn
  }
  if (length(warnings) != 0) {
    warning(warnings)
    return(warnings)
  }
}