## Title: repDilPCR - a Shiny App to Analyze qPCR Data by the Dilution-replicate Method
## File name: app.R
## Version: 1.0.4
## Date: 2021-10-22
## Author: Deyan Yordanov Yosifov
## Maintainer: Deyan Yordanov Yosifov <deyan.yosifov@uniklinik-ulm.de>
## Copyright: University Hospital Ulm, Germany, 2021
## License: GNU Affero General Public License, https://www.gnu.org/licenses/agpl-3.0.html

source("repDilPCR_lib.R")
library(shinycssloaders)
library(shinyalert)
options(spinner.type = 7)

# User interface ----
ui <- fluidPage(
  useShinyalert(),
  tags$style("[type = 'number'] {width: 100px;}"),
  titlePanel("repDilPCR"),
  sidebarLayout(
    sidebarPanel(

      # div(img(src="repDilPCR_logo.png", height = 110), img(src="UKU_Logo_RGB_Optinal_Farbe.png", height = 110, style = "padding-left:25px")),
      div(img(src="repDilPCR_logo.png", height = 110)),

      hr(style = "border-top: lpx solid #c0c0c0;"),

      fileInput(inputId = "input.table", label = h3("Upload your raw data as a .csv file", accept = ".csv")),

      numericInput("RG",
                   h3("Number of reference genes"),
                   value = 3),

      checkboxInput("impute", "Impute missing Cq values of reference genes", value = TRUE),
      h6(""), helpText("Do not use if your experiment does not contain replicates."),
      br(),
      actionButton("analyze", "Analyze", width = "100%", style="font-size:16pt; color: #fff; background-color: #5eba7d; border-color: #5e6600"),
      br(),
      h3("Optional parameters"),
      checkboxInput("statistics", "Test for statistically significant differences between samples or experimental groups", value = TRUE),

      textInput("ref.sample", h4("Reference sample or experimental group"), value = "default"),
      h6(""), helpText("Reference sample, in which gene expression will be regarded as 1 (100%) on linear scale and 0 on log2-scale, respectively. If left to \"default\", this will be the first sample in the table, resp. the leftmost sample on the plots. Change to the name of another sample (without a trailing underscore and replicate number) to make it the reference sample. If left empty, results will be shown in their original form, without forcing any particular sample to be 1 (100%) or 0."),

      radioButtons("test.type", h4("Type of statistical test(s)"),
                   choices = list("parametric", "non-parametric"), selected = "parametric"),

      radioButtons("posthoc", h4("Comparisons to test for statistically significant differences"),
                   choices = list("all to one (all to reference)" = "all to one", "all pairs", "selected pairs"), selected = "all to one"),

      splitLayout(numericInput("p",
                   h4("Significance level (\u03B1)"),
                   value = 0.05, min = 0, max = 1, step = 0.01),
                  numericInput("font.size",
                   h4("Font size of text on plots"),
                   value = 9, min = 2, step = 1)),

      radioButtons("sign.repr", h4("Display format of statistical significance on plots"),
                   choices = list("numeric p-values" = "values", "significance levels (asterisks)" = "asterisks"), selected = "values"),

      numericInput("sp.f", h5("Distance between significance bars on plots"), value = 1.5, min = 0, step = 0.5),
      h6(""), helpText("Spacing factor influencing the distance between significance bars on plots. In most cases, repDilPCR will succeed to distribute significance bars so that they will not overlap. If your significance bars overlap (which can be the case if you compare a lot of experimental groups), you can try increasing this heuristic parameter. Conversely, if the distances between significance bars are too big and they are wasting space on plots, you can try decreasing the factor. The default value is 1.5."),


      radioButtons("plot.format", h4("Format of graphical output (only for downloadable files)"),
                   choices = list("PDF", "PNG", "both", "no graphics" = "none"), selected = "PDF"),
      h5("Settings for PNG plots:"),
      splitLayout(#cellWidths = c("25%", "25%", "25%"),
                  numericInput("png.width", h5("Width in mm"), value = 190, min = 30),
                  numericInput("png.height", h5("Height in mm"), value = 134, min = 30),
                  numericInput("png_dpi", h5("Resolution in dpi"), value = 96, min = 72))


    ),
    mainPanel(
      tabsetPanel(id = "panels",
        tabPanel("Input data", tableOutput("inp.data")),
        tabPanel(title = "Preprocessed input data", shinycssloaders::withSpinner(tableOutput("qPCR")), value = 9),
        tabPanel(title = "Regression plots",
                 tabsetPanel(
                   tabPanel("Standard curves", uiOutput("stand.curves")),
                   tabPanel("Cq-Cq plots", uiOutput("Cq.plots"))
                 ), value = 10),
        tabPanel("Results",
                 tabsetPanel(
                   tabPanel("Plots in linear scale",
                            tabsetPanel(id = "linear",
                                        tabPanel(title = "Dot plots (all points)", uiOutput("p1"), value = 1),
                                        tabPanel(title = "Dot plots (means and confidence intervals)", uiOutput("p2"), value = 2),
                                        tabPanel(title = "Bar graphs (means and confidence intervals)", uiOutput("p3"), value = 3),
                                        tabPanel(title = "Box plots", uiOutput("p2n"), value = 4)
                                        # tabPanel("Cq-Cq plots", uiOutput("Cq.plots"))
                            )),
                   tabPanel("Plots in logarithmic scale",
                            tabsetPanel(id = "logarithmic",
                                        tabPanel(title = "Dot plots (all points)", uiOutput("p4"), value = 5),
                                        tabPanel(title = "Dot plots (means and standard deviations)", uiOutput("p5"), value = 6),
                                        tabPanel(title = "Bar graphs (means and standard deviations)", uiOutput("p6"), value = 7),
                                        tabPanel(title = "Box plots", uiOutput("p5n"), value = 8)
                                        )
                            ))
                 ),
        tabPanel("Download results",
                 tabsetPanel(id = "downloads",
                   tabPanel("Plots",
                            tabsetPanel(
                              tabPanel("Plots in linear scale",
                                       uiOutput("download.p1.c"),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p1")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p1.zip")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p2")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p2.zip")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p3")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p3.zip")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p2n")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p2n.zip"))
                                       ),
                              tabPanel("Plots in logarithmic scale",
                                       uiOutput("download.p4.c"),
                                       shinycssloaders::withSpinner(uiOutput("download.p4")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p4.zip")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p5")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p5.zip")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p6")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p6.zip")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p5n")),
                                       h6(""),
                                       shinycssloaders::withSpinner(uiOutput("download.p5n.zip"))
                                       ))),
                   tabPanel("Tables",
                            uiOutput("download.rel.q.detailed.c"),
                            shinycssloaders::withSpinner(uiOutput("download.rel.q.detailed")),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.rel.q.detailed.log")),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.rel.q.mean")),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.rel.q.mean.log"))
                   ),
                   tabPanel(title = "Intermediate data",
                            uiOutput("download.qPCR.b"),
                            uiOutput("download.qPCR"),
                            uiOutput("download.qPCR.c"),
                            uiOutput("download.eff.df"),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.stand.curves")),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.stand.curves.zip")),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.Cq.plots")),
                            h6(""),
                            shinycssloaders::withSpinner(uiOutput("download.Cq.plots.zip")), value = 11
                   )
                 )
        ),
        tabPanel("About", h3(tagList("The repDilPCR program was written by Deyan Yordanov Yosifov at the Department of Internal Medicine III of the University Hospital in Ulm, Germany.
                                     The program is inspired by the dilution-replicate approach for design and analysis of real-time PCR assays (Kwokyin Hui & Zhong-Ping Feng (2013)
                                     Efficient experimental design and analysis of real-time PCR assays, Channels, 7:3, 160-170, DOI: ", a("10.4161/chan.24024", href = "https://doi.org/10.4161/chan.24024"), ")."), style="font-size:12pt"),
                          h3(tagList("\u2003")),
                          h3(tagList("Overview:"), style="font-size:12pt"),
                          h3(tagList("In a qPCR experiment, it is of key importance to determine the efficiency of the PCR reaction for each amplicon and primer pair for correct evaluation and interpretation of the data.
                                     Different approaches to determine efficiency have been developed, from the classical calibration curve-based method to sophisticated methods that rely on fitting
                                     linear or non-linear models on individual amplification curves. Occupying the middle ground between these two extremes is the dilution-replicate experimental
                                     design of Hui and Feng that has remained somehow overlooked, most probably due to the lack up to now of a dedicated software tool to apply the method. This is a multiple
                                     linear regression-based approach with a number of advantages. It requires fewer reactions than the traditional approach with calibration curves produced by a
                                     separate set of dilutions of a standard sample. In the dilution-replicate design, standard curves are determined from so-called dilution-replicates of
                                     experimental reactions that serve both to control technical variance and to determine efficiency. Like this, all samples contribute to the efficiency estimate,
                                     thus precision increases with the number of samples on a plate. Furthermore, the traditional approach requires that the linear dynamic range of the independent
                                     standard curve covers all sample Cq values which sometimes leads to the necessity to repeat experiments using different dilutions. In contrast, with the
                                     dilution-replicate design it is guaranteed that the sample Cq values will be within range."), style="font-size:12pt"),
                          h3(tagList("repDilPCR utilizes the described dilution-replicate analytical method and extends it by adding the possibility to use multiple reference genes. It also offers capabilities for
                                     performing statistical tests and plotting publication-ready graphs. The program has been designed with the philosophy to automate and speed up analysis of
                                     qPCR data (typically less than one minute from raw Cq values to publication-ready plots) and to help users with little knowledge of statistics to select and
                                     perform the appropriate statistical tests, at least in the case of one-factor experimental designs. At the same time, the program allows experienced users to
                                     export intermediate data and perform more sophisticated analyses with external statistical software, e.g. if two-way ANOVA is necessary."), style="font-size:12pt"),
                          h3(tagList("Detailed user manual can be found ", a("here", href = "https://gitfront.io/r/deyanyosifov/c8e8e53b2f70690abb47d0847dabe55c31e73afd/repDilPCR/"), " (still in the process of writing). Potential users are advised to read the ", a("article", href = "https://doi.org/10.4161/chan.24024")," by Hui and Feng before setting up their experiment and using the program."), style="font-size:12pt"),
                          h3(tagList("\u2003")),
                          uiOutput("download.test.data"),
                          h6(""),
                          uiOutput("download.test.data.precalc"),
                          h6(textOutput("count"), style="font-size:10pt; color: #fff"))
      )
    )
  )
)

# Server logic ----
server <- function(input, output, session) {

  # Prepare data ----
  inp.data <- reactive({
    file <- input$input.table
    ext <- tools::file_ext(file$datapath)
    req(file)
    validate(need(ext == "csv", "Please upload a csv file"))
    inp.data <- read.csv(file$datapath, sep = ",", dec = ".", stringsAsFactors = TRUE)
    if (!file.exists("counter.txt")) {
      write.table(0, "counter.txt", row.names = FALSE, col.names = FALSE)
    }
    counter <- read.table("counter.txt", stringsAsFactors = FALSE)
    counter[1,1] <- counter[1,1] + 1
    write.table(counter, "counter.txt", row.names = FALSE, col.names = FALSE)
    if (colnames(inp.data)[3] == "Dilution") {
      rd.preprocess.output <- rd.preprocess(inp.data, input$RG)
    } else {
      rd.preprocess.output <- rd.preprocess.2(inp.data)
    }
  })

  output$inp.data <- renderTable({
    inp.data()$qPCR
  })

  # Impute missing Cq values of reference genes ----
  qPCR <- reactive({
    if (colnames(inp.data()$qPCR)[3] == "Dilution" && input$impute == TRUE) {rd.impute(inp.data()$qPCR, inp.data()$ref.genes)} else {inp.data()$qPCR}
  })

  output$qPCR <- renderTable({
    qPCR()
  })

  # Make a virtual reference gene by averaging the Cq values of the real reference genes at each sample/dilution/replicate combination ----
  qPCR.NF <- reactive({
    if (input$analyze > 0) {
      if (colnames(inp.data()$qPCR)[3] == "Dilution") {
        rd.ref(qPCR(), input$RG)
      } else {
        qPCR()
        }
    }
  })

  update.qPCR <- observeEvent(qPCR.NF(), {
    output$qPCR <- renderTable({qPCR.NF()})
  })

  observe({
    if (colnames(inp.data()$qPCR)[3] != "Dilution") {
      hideTab(inputId = "panels", target = "9")
      hideTab(inputId = "panels", target = "10")
      hideTab(inputId = "downloads", target = "11")
    } else {
      showTab(inputId = "panels", target = "9")
      showTab(inputId = "panels", target = "10")
      showTab(inputId = "downloads", target = "11")
    }
  })


  # Multiple linear regression for standard curves ----
  model.list <- reactive({
    req(qPCR.NF())
    if (colnames(inp.data()$qPCR)[3] == "Dilution"){
      rd.mlr.model(qPCR.NF(), inp.data()$all.genes)
    }
  })

  # Calculate efficiencies ----
  eff.df <- reactive({
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      rd.eff(model.list(), inp.data()$all.genes)
    }
  })

  # Plot multiple regressions with separate regression curves but common slope (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), model.list(), eff.df())},{
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      output$stand.curves <- renderUI({
        stand_curves_list <- lapply(1:length(inp.data()$all.genes), function(i) {
          plotname <- paste("plot", i, sep="")
          shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
        })
        do.call(tagList, stand_curves_list)
      })

      for (i in 1:length(inp.data()$all.genes)) {
        local({
          stand_i <- i
          plotname <- paste("plot", stand_i, sep="")
          output[[plotname]] <- renderPlot({
            rd.plot.mlr.s(qPCR = qPCR.NF(), model.list = model.list(), eff.df = eff.df(), all.genes = inp.data()$all.genes, font.size = input$font.size + 4)[[stand_i]]
          })
        })
      }}
  })

  # Plot multiple regressions with separate regression curves but common slope (as files for download) ----
  stand.curves.output <- reactive({
    req(qPCR.NF())
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      rd.plot.mlr(qPCR = qPCR.NF(), model.list = model.list(), eff.df = eff.df(), all.genes = inp.data()$all.genes, font.size = input$font.size)
    }
  })

  # Multiple linear regression for Cq-Cq plots ----
  Cq.list <- reactive({
    req(qPCR.NF())
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      rd.Cq.Cq(qPCR.NF(), inp.data()$GOIs)
    }
  })

  ## Plot Cq-Cq plots for display in a web browser ----
  observeEvent({c(qPCR.NF(), Cq.list())},{
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      output$Cq.plots <- renderUI({
        Cq.plots.list <- lapply(1:length(inp.data()$GOIs), function(i) {
          Cq.plotname <- paste("Cq.plot", i, sep="")
          shinycssloaders::withSpinner(plotOutput(Cq.plotname, height = 1200))
        })
        do.call(tagList, Cq.plots.list)
      })

      for (i in 1:length(inp.data()$GOIs)) {
        local({
          Cq_i <- i
          Cq.plotname <- paste("Cq.plot", Cq_i, sep="")
          output[[Cq.plotname]] <- renderPlot({
            rd.plot.Cq.Cq.s(qPCR = qPCR.NF(), Cq.list = Cq.list(), GOIs = inp.data()$GOIs, font.size = input$font.size + 4)[[Cq_i]]
          })
        })
      }}
  })

  ## Plot Cq-Cq plots as files for download ----
  Cq.plots.output <- reactive({
    req(qPCR.NF())
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      rd.plot.Cq.Cq(qPCR = qPCR.NF(), Cq.list = Cq.list(), GOIs = inp.data()$GOIs, font.size = input$font.size)
    }
  })

  ## Relative quantities ----
  rel.q.results <- reactive({
    req(qPCR.NF())
    if (colnames(inp.data()$qPCR)[3] == "Dilution") {
      rd.rel.quant(qPCR = qPCR.NF(), Cq.list = Cq.list(), eff.df = eff.df(), GOIs = inp.data()$GOIs)
    } else {
    rd.rel.quant.2(qPCR = qPCR.NF(), GOIs = inp.data()$GOIs)
    }
  })

  ## Logarithmic tables ----
  rel.q.results.log <- reactive({
    rd.log(rel.q.detailed = rel.q.results()$rel.q.detailed, rel.q.df = rel.q.results()$rel.q.df)
  })

  # Normalize relative quantities according to expression levels in a selected sample (if any) ----
  rel.q.norm.results.c <- reactive({c(rel.q.results(), rel.q.results.log())})
  rel.q.norm.results <- reactive({
    rd.normalize(rel.q.detailed = rel.q.norm.results.c()$rel.q.detailed, rel.q.detailed.log = rel.q.norm.results.c()$rel.q.detailed.log, rel.q.df = rel.q.norm.results.c()$rel.q.df, rel.q.log = rel.q.norm.results.c()$rel.q.log, rel.q.mean = rel.q.norm.results.c()$rel.q.mean, rel.q.mean.log = rel.q.norm.results.c()$rel.q.mean.log, ref.sample = isolate(input$ref.sample), GOIs = inp.data()$GOIs)
  })

  ## Calculate confidence intervals
  rel.q.confint <- reactive({
    rd.confint(rel.q.mean = rel.q.norm.results()$rel.q.mean, rel.q.mean.log = rel.q.norm.results()$rel.q.mean.log, p = input$p)
  })

  ## Statistical tests
  statistics.results.c <- reactive({rel.q.norm.results()})
  statistics.results <- reactive({
    req(statistics.results.c())
    rd.statistics(rel.q.df = statistics.results.c()$rel.q.df, rel.q.log = statistics.results.c()$rel.q.log, rel.q.mean = rel.q.confint()$rel.q.mean, rel.q.mean.log = rel.q.confint()$rel.q.mean.log, statistics = input$statistics, test.type = input$test.type, posthoc = input$posthoc, ref.sample = statistics.results.c()$ref.sample, p = input$p, sp.f = input$sp.f)
  })

  # Plot relative expression dotplots in linear scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p1 <- renderUI({
      p1.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p1.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p1.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p1_i <- i
        plotname <- paste("p1.plot", p1_i, sep="")
        output[[plotname]] <- renderPlot({
          rd.plot.p1(rel.q.df = statistics.results()$rel.q.df, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p1[[p1_i]]
        })
      })
    }
  })

  # Plot relative expression dotplots (means + CI) in linear scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p2 <- renderUI({
      p2.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p2.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p2.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p2_i <- i
        plotname <- paste("p2.plot", p2_i, sep="")
        output[[plotname]] <- renderPlot({
          req(statistics.results(), input$test.type == "parametric")
          rd.plot.p2(rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p2[[p2_i]]
        })
      })
    }
  })

  # Plot relative expression bar graphs in linear scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p3 <- renderUI({
      p3.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p3.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p3.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p3_i <- i
        plotname <- paste("p3.plot", p3_i, sep="")
        output[[plotname]] <- renderPlot({
          req(statistics.results(), input$test.type == "parametric")
          rd.plot.p3(rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p3[[p3_i]]
        })
      })
    }
  })

  # Plot relative expression box plots in linear scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p2n <- renderUI({
      p2n.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p2n.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p2n.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p2n_i <- i
        plotname <- paste("p2n.plot", p2n_i, sep="")
        output[[plotname]] <- renderPlot({
          req(statistics.results(), input$test.type == "non-parametric")
          rd.plot.p2n(rel.q.df = statistics.results()$rel.q.df, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p2n[[p2n_i]]
        })
      })
    }
  })


  observe({
    if (input$test.type == "non-parametric") {
      hideTab(inputId = "linear", target = "2")
      hideTab(inputId = "linear", target = "3")
      showTab(inputId = "linear", target = "4")
    } else {
      showTab(inputId = "linear", target = "2")
      showTab(inputId = "linear", target = "3")
      hideTab(inputId = "linear", target = "4")
    }
  })


  # Plot relative expression dotplots in logarithmic scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p4 <- renderUI({
      p4.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p4.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p4.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p4_i <- i
        plotname <- paste("p4.plot", p4_i, sep="")
        output[[plotname]] <- renderPlot({
          rd.plot.p4(rel.q.log = statistics.results()$rel.q.log, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p4[[p4_i]]
        })
      })
    }
  })

  # Plot relative expression dotplots (means + SD) in logarithmic scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p5 <- renderUI({
      p5.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p5.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p5.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p5_i <- i
        plotname <- paste("p5.plot", p5_i, sep="")
        output[[plotname]] <- renderPlot({
          req(statistics.results(), input$test.type == "parametric")
          rd.plot.p5(rel.q.mean.log = statistics.results()$rel.q.mean.log, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p5[[p5_i]]
        })
      })
    }
  })

  # Plot relative expression bar graphs in logarithmic scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p6 <- renderUI({
      p6.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p6.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p6.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p6_i <- i
        plotname <- paste("p6.plot", p6_i, sep="")
        output[[plotname]] <- renderPlot({
          req(statistics.results(), input$test.type == "parametric")
          rd.plot.p6(rel.q.mean.log = statistics.results()$rel.q.mean.log, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p6[[p6_i]]
        })
      })
    }
  })

  # Plot relative expression box plots in logarithmic scale (for display in a web browser) ----
  observeEvent({c(qPCR.NF(), statistics.results())},{
    output$p5n <- renderUI({
      p5n.results_list <- lapply(inp.data()$GOIs, function(i) {
        plotname <- paste("p5n.plot", i, sep="")
        shinycssloaders::withSpinner(plotOutput(plotname, height = 600))
      })
      do.call(tagList, p5n.results_list)
    })

    for (i in inp.data()$GOIs) {
      local({
        p5n_i <- i
        plotname <- paste("p5n.plot", p5n_i, sep="")
        output[[plotname]] <- renderPlot({
          req(statistics.results(), input$test.type == "non-parametric")
          rd.plot.p5n(rel.q.log = statistics.results()$rel.q.log, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size + 4)$p5n[[p5n_i]]
        })
      })
    }
  })

  observe({
    if (input$test.type == "non-parametric") {
      hideTab(inputId = "logarithmic", target = "6")
      hideTab(inputId = "logarithmic", target = "7")
      showTab(inputId = "logarithmic", target = "8")
    } else {
      showTab(inputId = "logarithmic", target = "6")
      showTab(inputId = "logarithmic", target = "7")
      hideTab(inputId = "logarithmic", target = "8")
    }
  })


  # Downloadable csv of input data with imputed missing Cq values of reference genes ----
  observeEvent(qPCR.NF(), {
  output$download.qPCR1 <- downloadHandler(
    filename = function() {
      paste0(gsub(".csv", "", input$input.table), "_with_imputed_missing_values.csv")
    },
    content = function(file) {
      if (input$analyze > 0) {
      write.csv(qPCR.NF(), file)
      } else {
        write.csv(qPCR(), file)
      }
    }
  )
  }, ignoreNULL = FALSE)

  output$download.qPCR <- renderUI({
    req(input$input.table, qPCR())
    downloadButton("download.qPCR1", "Download table with imputed Cq values (if any)")
  })

  output$download.qPCR.b <- renderUI({
    req(input$input.table, qPCR())
    h6("")
    helpText("Download may start a few seconds after pressing the respective button.")
  })


  output$download.qPCR.c <- renderUI({
    req(input$input.table, qPCR())
    h6("")
    helpText("This is a table of preprocessed input data,
                    including imputed missing Cq values of reference genes
                    if the respective option was selected and if indeed the input
                    data contained missing values. The imputation is done
                    by the R package \"mice\" based on the available
                    Cq values of reference genes from replicates with
                    non-missing data using the following parameters:
                    imputation method = \"midastouch\" (weighted predictive
                    mean matching), number of multiple imputations = 20,
                    number of iterations = 20. The final number is the average of
                    the 20 imputations.")
  })


  # Downloadable csv of calculated reaction efficiencies ----
  output$download.eff.df1 <- downloadHandler(
    filename = function() {
      paste0(gsub(".csv", "", input$input.table), "_efficiency_table.csv")
    },
    content = function(file) {
      write.csv(eff.df(), file)
    }
  )

  output$download.eff.df <- renderUI({
    req(input$input.table, eff.df())
    downloadButton("download.eff.df1", "Download table with calculated efficiencies")
  })

  # Download a pdf file of standard curves ----
  output$download.stand.curves1 <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_standard_curves.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", stand.curves.output()$ml)
      }
    }
  )

  output$download.stand.curves <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, stand.curves.output()$ml)
      downloadButton("download.stand.curves1", "Download standard curves (PDF)")
    }
  })

  # Download a zip archive of standard curves ----
  output$download.stand.curves2 <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_standard_curves.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.1(ggplot.object = stand.curves.output()$stand.curves, fname = input$input.table$name, type.name = "standard_curves", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.stand.curves.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, stand.curves.output()$stand.curves)
      downloadButton("download.stand.curves2", "Download standard curves (ZIP archive of PNG files)")
    }
  })

  # Download a pdf file of Cq-Cq plots ----
  output$download.Cq.plots1 <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_Cq-Cq_plots.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", Cq.plots.output()$ml)
      }
    }
  )

  output$download.Cq.plots <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, Cq.plots.output()$ml)
      downloadButton("download.Cq.plots1", "Download Cq-Cq plots (PDF)")
    }
  })

  # Download a zip archive of Cq-Cq plots ----
  output$download.Cq.plots2 <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_Cq-Cq_plots.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.1(ggplot.object = Cq.plots.output()$Cq.plots, fname = input$input.table$name, type.name = "Cq-Cq_plots", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.Cq.plots.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, Cq.plots.output()$Cq.plots)
      downloadButton("download.Cq.plots2", "Download Cq-Cq plots (ZIP archive of PNG files)")
    }
  })

  # Download dotplots  in linear scale ----
  p1.results <- reactive({
    rd.plot.p1(rel.q.df = statistics.results()$rel.q.df, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
  })

  p2.results <- reactive({
    if (input$test.type == "parametric") {
    rd.plot.p2(rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
    }
  })

  p3.results <- reactive({
    if (input$test.type == "parametric") {
      rd.plot.p3(rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
    }
  })

  p2n.results <- reactive({
    if (input$test.type == "non-parametric") {
      rd.plot.p2n(rel.q.df = statistics.results()$rel.q.df, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
    }
  })

  output$download.p1.c <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p1.results())
      h6("")
      helpText("Download may start a few seconds after pressing the respective button.")
    }
  })



    # As a pdf file ----
  output$download.p1t <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_dotplot.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p1.results()$ml)
      }
    }
  )

  output$download.p1 <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p1.results()$ml)
      downloadButton("download.p1t", "Download dot plots (all points) in linear scale (PDF)")
    }
  })


  output$download.p2t <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_dotplot_CI.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p2.results()$ml)
      }
    }
  )

  output$download.p2 <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p2.results()$ml)
      downloadButton("download.p2t", "Download dot plots (means and confidence intervals) in linear scale (PDF)")
    }
  })

  output$download.p3t <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_bar_graph.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p3.results()$ml)
      }
    }
  )

  output$download.p3 <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p3.results()$ml)
      downloadButton("download.p3t", "Download bar graphs in linear scale (PDF)")
    }
  })

  output$download.p2nt <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_boxplot.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p2n.results()$ml)
      }
    }
  )

  output$download.p2n <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p2n.results()$ml)
      downloadButton("download.p2nt", "Download box plots in linear scale (PDF)")
    }
  })


    # As a zip file of png files ----
  output$download.p1z <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_dotplot.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p1.results()$p1, fname = input$input.table$name, type.name = "relative_expression_dotplot", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p1.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p1.results()$p1)
      downloadButton("download.p1z", "Download dot plots (all points) in linear scale (ZIP archive of PNG files)")
    }
  })

  output$download.p2z <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_dotplot_CI.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p2.results()$p2, fname = input$input.table$name, type.name = "relative_expression_dotplot_CI", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p2.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p2.results()$p2)
      downloadButton("download.p2z", "Download dot plots (means and confidence intervals) in linear scale (ZIP archive of PNG files)")
    }
  })

  output$download.p3z <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_bar_graph.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p3.results()$p3, fname = input$input.table$name, type.name = "relative_expression_bar_graph", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p3.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p3.results()$p3)
      downloadButton("download.p3z", "Download bar graphs in linear scale (ZIP archive of PNG files)")
    }
  })

  output$download.p2nz <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_expression_boxplot.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p2n.results()$p2n, fname = input$input.table$name, type.name = "relative_expression_boxplot", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p2n.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p2n.results()$p2n)
      downloadButton("download.p2nz", "Download box plots in linear scale (ZIP archive of PNG files)")
    }
  })

  # Download dotplots  in logarithmic scale ----
  p4.results <- reactive({
    rd.plot.p4(rel.q.log = statistics.results()$rel.q.log, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
  })

  p5.results <- reactive({
    if (input$test.type == "parametric") {
      rd.plot.p5(rel.q.mean.log = statistics.results()$rel.q.mean.log, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
    }
  })

  p6.results <- reactive({
    if (input$test.type == "parametric") {
      rd.plot.p6(rel.q.mean.log = statistics.results()$rel.q.mean.log, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
    }
  })

  p5n.results <- reactive({
    if (input$test.type == "non-parametric") {
      rd.plot.p5n(rel.q.log = statistics.results()$rel.q.log, rel.q.mean = statistics.results()$rel.q.mean, res.posthoc = statistics.results()$res.posthoc, ref.sample = statistics.results()$ref.sample, GOIs = inp.data()$GOIs, statistics = statistics.results()$statistics, posthoc = input$posthoc, sign.repr = input$sign.repr, p = input$p, stat.test = statistics.results()$stat.test, font.size = input$font.size)
    }
  })

  output$download.p4.c <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p4.results())
      h6("")
      helpText("Download may start a few seconds after pressing the respective button.")
    }
  })



  # As a pdf file ----
  output$download.p4t <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_dotplot.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p4.results()$ml)
      }
    }
  )

  output$download.p4 <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p4.results()$ml)
      downloadButton("download.p4t", "Download dot plots (all points) in logarithmic scale (PDF)")
    }
  })


  output$download.p5t <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_dotplot_SD.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p5.results()$ml)
      }
    }
  )

  output$download.p5 <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p5.results()$ml)
      downloadButton("download.p5t", "Download dot plots (means and standard deviations) in logarithmic scale (PDF)")
    }
  })

  output$download.p6t <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_bar_graph.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p6.results()$ml)
      }
    }
  )

  output$download.p6 <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p6.results()$ml)
      downloadButton("download.p6t", "Download bar graphs in logarithmic scale (PDF)")
    }
  })

  output$download.p5nt <- downloadHandler(
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_boxplot.pdf")
      }
    },
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      content = function(file) {
        ggplot2::ggsave(file, width = 210, height = 297, units = "mm", p5n.results()$ml)
      }
    }
  )

  output$download.p5n <- renderUI({
    if (input$plot.format == "PDF" | input$plot.format == "both") {
      req(input$input.table, p5n.results()$ml)
      downloadButton("download.p5nt", "Download box plots in logarithmic scale (PDF)")
    }
  })

  # As a zip file of png files ----
  output$download.p4z <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_dotplot.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p4.results()$p4, fname = input$input.table$name, type.name = "relative_log_expression_dotplot", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p4.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p4.results()$p4)
      downloadButton("download.p4z", "Download dot plots (all points) in logarithmic scale (ZIP archive of PNG files)")
    }
  })

  output$download.p5z <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_dotplot_SD.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p5.results()$p5, fname = input$input.table$name, type.name = "relative_log_expression_dotplot_SD", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p5.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p5.results()$p5)
      downloadButton("download.p5z", "Download dot plots (means and standard deviations) in logarithmic scale (ZIP archive of PNG files)")
    }
  })

  output$download.p6z <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_bar_graph.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p6.results()$p6, fname = input$input.table$name, type.name = "relative_log_expression_bar_graph", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p6.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p6.results()$p6)
      downloadButton("download.p6z", "Download bar graphs in logarithmic scale (ZIP archive of PNG files)")
    }
  })

  output$download.p5nz <- downloadHandler(
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      filename = function() {
        paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_boxplot.zip")
      }
    },
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      content = function(file) {
        fs <- png.plot.2(ggplot.object = p5n.results()$p5n, fname = input$input.table$name, type.name = "relative_log_expression_boxplot", png.size = c(input$png.width, input$png.height), png.dpi = input$png_dpi)
        zip(zipfile = file, files = fs, flags = "-j")
      }
    }, contentType = "application/zip"
  )

  output$download.p5n.zip <- renderUI({
    if (input$plot.format == "PNG" | input$plot.format == "both") {
      req(input$input.table, p5n.results()$p5n)
      downloadButton("download.p5nz", "Download box plots in logarithmic scale (ZIP archive of PNG files)")
    }
  })

  ## Downloadable tables of results ----
  save.tables <- reactive({
    rd.save.tables(input.table = input$input.table$name, rel.q.detailed = rel.q.norm.results()$rel.q.detailed, rel.q.detailed.log = rel.q.norm.results()$rel.q.detailed.log, rel.q.mean = statistics.results()$rel.q.mean, rel.q.mean.log = statistics.results()$rel.q.mean.log, p = input$p)
  })

  output$download.rel.q.detailed.c <- renderUI({
    req(input$input.table, rel.q.norm.results())
    h6("")
    helpText("Download may start a few seconds after pressing the respective button.")
  })

  # Downloadable csv of calculated relative expression values for each replicate ----
  output$download.rel.q.detailed1 <- downloadHandler(
    filename = function() {
      paste0(gsub(".csv", "", input$input.table), "_relative_expression_replicates.csv")
    },
    content = function(file) {
      write.csv(rel.q.norm.results()$rel.q.detailed, file)
    }
  )

  output$download.rel.q.detailed <- renderUI({
    req(input$input.table, rel.q.norm.results()$rel.q.detailed, save.tables())
    downloadButton("download.rel.q.detailed1", "Download table with calculated relative expression values for each replicate (linear scale)")
  })


  # Downloadable csv of calculated logarithmic relative expression values for each replicate ----
  output$download.rel.q.detailed.log1 <- downloadHandler(
    filename = function() {
      paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_replicates.csv")
    },
    content = function(file) {
      write.csv(save.tables()$rel.q.detailed.log, file)
    }
  )

  output$download.rel.q.detailed.log <- renderUI({
    req(input$input.table, save.tables()$rel.q.detailed.log)
    downloadButton("download.rel.q.detailed.log1", "Download table with calculated relative expression values for each replicate (logarithmic scale)")
  })

  # Downloadable csv of calculated mean relative expression values ----
  output$download.rel.q.mean1 <- downloadHandler(
    filename = function() {
      paste0(gsub(".csv", "", input$input.table), "_relative_expression_averaged.csv")
    },
    content = function(file) {
      write.csv(save.tables()$rel.q.mean, file, row.names = FALSE)
    }
  )

  output$download.rel.q.mean <- renderUI({
    req(input$input.table, save.tables()$rel.q.mean, )
    downloadButton("download.rel.q.mean1", "Download table with calculated relative mean expression values (linear scale)")
  })

  # Downloadable csv of calculated mean logarithmic relative expression values ----
  output$download.rel.q.mean.log1 <- downloadHandler(
    filename = function() {
      paste0(gsub(".csv", "", input$input.table), "_relative_log_expression_averaged.csv")
    },
    content = function(file) {
      write.csv(save.tables()$rel.q.mean.log, file, row.names = FALSE)
    }
  )

  output$download.rel.q.mean.log <- renderUI({
    req(input$input.table, save.tables()$rel.q.mean.log)
    downloadButton("download.rel.q.mean.log1", "Download table with calculated relative mean expression values (logarithmic scale)")
  })

  # Downloadable test data ----
  output$download.test.data.1 <- downloadHandler(
    filename = c("Test_data.csv"),
    content = function(file) {
      file.copy("Test_data.csv", file)
      }
  )

  output$download.test.data <- renderUI({
    downloadButton("download.test.data.1", "Download exemplary test data for use with the program (dilution-replicate setup)")
  })
  
  output$download.test.data.2 <- downloadHandler(
    filename = c("Test_data_precalc.csv"),
    content = function(file) {
      file.copy("Test_data_precalc.csv", file)
    }
  )

  output$download.test.data.precalc <- renderUI({
    downloadButton("download.test.data.2", "Download exemplary preprocessed test data for use with the program (independent of experimental setup)")
  })
  
  
  ## Counter of uploaded files
  output$count <- renderText({
    if (file.exists("counter.txt")) {
    cc <- read.table("counter.txt", stringsAsFactors = FALSE)
    cc[1,1]
    }
  })

  ## Pop-up window with terms of use
  shinyalert(
    title = "<h2 style='text-align:center;'><b>Data protection & Terms of Use / Datenschutzerklärung und Nutzungsbedingungen</b></h2>",
    text = tags$div(style = 'overflow:auto;max-height:50vh;text-align:justify;padding-right:10px;',
                    h5(style="color:Tomato;","Scroll down for the German version / Scrollen Sie bitte nach unten für die deutsche Version"),
                    h3(tags$b("Disclaimer and general terms of use")),
                    h5("The repDilPCR Shiny App is a tool for scientific analysis, visualization and interpretation of qPCR data. The tool has not been clinically validated; the use of data and results generated by the repDilPCR app is for scientific purposes only. Use for diagnostic or therapeutic purposes and/or implementation of the results in clinical studies is solely the responsibility of the user and the treating physician."),
                    h5("The following applies for the processing of all data by the repDilPCR Shiny App:",
                       tags$ul(tags$li("The upload of raw data is done exclusively by the user himself/herself, he/she alone is responsible for the technical quality (execution of laboratory methods, i.e. isolation of nucleic acids, reverse transcription and polymerase chain reaction) and correctness of the uploaded data."),
                               tags$li("The data uploaded by the user are stored temporarily for the duration of the processing only, there is no permanent data storage."),
                               tags$li("The results generated by the repDilPCR Shiny App can be viewed by the user within the app itself or downloaded as CSV, PDF and/or PNG files."),
                               tags$li("No data and no results are stored in the long term or passed on to third parties. Results are available exclusively for the user as a download. If the current session of the Shiny App is terminated (closing of the browser window, inactivity, timeout, ...), all uploaded and generated files of this session will be deleted from the server."),
                               tags$li("Each user is permitted to access and see only the data uploaded and processed by him/her in the respective session."),
                               tags$li("Access to data of other users/session IDs is not allowed."),
                               tags$li("The repDilPCR Shiny App is a purely scientific tool, which is solely meant to be used for processing, analysis and visualization of raw or pre-processed data for scientific purposes. Clinical validation has not been performed."),
                               tags$li("Further use of the results generated by the repDilPCR Shiny App, whether for basic research purposes, as information for treating physicians, for diagnostics, for therapy recommendation, or for the purpose of study design, is solely the responsibility of the user and/or the treating physician."))
                    ),
                    h3(tags$b("Imprint")),
                    h4(tags$b("Name and address")),
                    h5("The responsible party in the sense of the EU General Data Protection Regulation and other national data protection laws of the member states as well as other legal data protection regulations is:"),
                    h5("German Cancer Research Center - Foundation under Public Law"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Germany"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:kontakt@dkfz.de", "kontakt@dkfz.de")),
                    h5("Website: ", a(href = "https://www.dkfz.de", "www.dkfz.de")),
                    h2(""),
                    h4(tags$b("Name and address of the data protection officer")),
                    h5("Data Protection Officer"),
                    h5("German Cancer Research Center - Foundation under Public Law"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Germany"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:datenschutz@dkfz.de", "datenschutz@dkfz.de")),
                    h2(""),
                    h4(tags$b("Legal form")),
                    h5("The German Cancer Research Center is a foundation under public law of the State of Baden-Württemberg. The statutes of May 2018 are valid. Sales tax identification number: DE 143293537"),
                    h2(""),
                    h4(tags$b("Authorized representatives")),
                    h5("Prof. Dr. Michael Baumann (Chairman of the Board and Scientific Director)"),
                    h5("Ursula Weyrich (Commercial Director)"),
                    h5("(addresses as above)"),
                    h2(""),
                    h4(tags$b("Responsible for the content")),
                    h5("Address as above"),
                    h2(""),
                    h4(tags$b("Copyright notice")),
                    h5("All contents of the repDilPCR Shiny App are licensed under the GNU General Public License v3 or the GNU Affero General Public License v3 as applicable."),
                    h2(""),
                    h4(tags$b("Liability notice")),
                    h5("Despite careful control of the content, we do not assume any liability for the content as well as for the processed results. The use of the results generated by the repDilPCR Shiny App, be it for basic research purposes, as information for treating physicians, for diagnostics, for therapy recommendation or for the purpose of study designs etc. is solely the responsibility of the user and/or the treating physician."),
                    h2(""),
                    h4(tags$b("Scope")),
                    h5('These Terms of Use apply to your use of repDilPCR, including all content, features, services and rules set forth herein. We reserve the right to offer additional services on repDilPCR or to update repDilPCR on a regular basis. These Terms of Use are intended for both consumers and businesses. Both are grouped here under "Users". Consumers are natural persons who conclude a legal transaction for a purpose that can be attributed to neither their commercial nor their independent professional activity. Entrepreneurs are natural or legal persons or partnerships with legal capacity who, when concluding a legal transaction, act in the exercise of their commercial or independent professional activity.'),
                    h3(tags$b("Data protection")),
                    h4(tags$b("Name and address")),
                    h5("The responsible party in the sense of the EU General Data Protection Regulation and other national data protection laws of the member states as well as other legal data protection regulations is:"),
                    h5("German Cancer Research Center - Foundation under Public Law"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Germany"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:kontakt@dkfz.de", "kontakt@dkfz.de")),
                    h5("Website: ", a(href = "https://www.dkfz.de", "www.dkfz.de")),
                    h2(""),
                    h4(tags$b("Name and address of the data protection officer")),
                    h5("Data Protection Officer"),
                    h5("German Cancer Research Center - Foundation under Public Law"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Germany"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:datenschutz@dkfz.de", "datenschutz@dkfz.de")),
                    h2(""),
                    h4(tags$b("General information on data processing")),
                    h4(tags$i("Scope of the processing of personal data")),
                    h5("No information about the user is required. When opening the repDilPCR Shiny app, a session ID is generated, under which the upload of the scientific data, the processing as well as the download of the generated results takes place. There is no data bundling. The uploaded data is only stored temporarily for the duration of the processing and deleted after the session ends (see data deletion below). There is no long-term storage. No customer data or personal data is entered or stored."),
                    h2(""),
                    h5(tags$i("Scientific data")),
                    h5("For data originating from human source material applies:"),
                    h5("Data uploaded by the user to the repDilPCR Shiny App must be fully anonymized or pseudonymized. The uploaded data are stored only temporarily in the cache of the Shiny App to generate the results. Identification of persons based on the uploaded documents is not allowed. The user agrees to upload only anonymized or fully pseudonymized data to repDilPCR. No conclusion can be drawn about the origin of the data."),
                    h5("If data from human primary material are used, the user is responsible for obtaining patient or donor consent. This must allow the use of the primary material and further analysis/processing of the data for scientific research purposes."),
                    h5("An exception applies in those cases in which it is not possible to obtain prior consent for factual reasons and the processing of the data is permitted by legal regulations. The user of the Shiny app is solely responsible for this. The user agrees to these terms by giving his/her consent in the pop-up window that appears when opening the Shiny App."),
                    h5("The repDilPCR Shiny App is a tool exclusively for scientific analysis, visualization and interpretation of the processed data. The app has not been clinically validated; the use of the data and results generated by the repDilPCR App is for scientific purposes only. Use for diagnostic or therapeutic purposes and/or implementation of the results in clinical studies is solely the responsibility of the user and the treating physician."),
                    h5("When opening the Shiny app, the user is informed of the applicable terms of use as well as of this privacy policy via a pop-up window."),
                    h5("The pop-up window contains the explicit notice that in the case of data generated from human source material (patient material), these may only be uploaded and processed in anonymized or pseudonymized form and a declaration of consent must be available from the patient. The declaration must allow the use of the primary material as well as further analysis/processing of the data for research purposes."),
                    h5("The user accepts these conditions as well as the privacy policy by agreeing to the terms and conditions by clicking on the corresponding field within the pop-up window. Only then can data be uploaded and processed."),
                    h2(""),
                    h4(tags$i("Legal basis for the processing of personal data")),
                    h5("Insofar as data are processed that originate from human source material, the user is responsible for obtaining the declaration of consent from the respective patient(s)/donor(s) (see above). Art. 6 para. 1 point a of the EU General Data Protection Regulation (GDPR) serves as the legal basis for the processing of personal data."),
                    h5("When opening the app, the user agrees to process only anonymized or pseudonymized data."),
                    h2(""),
                    h4(tags$i("Data deletion and storage duration")),
                    h5("No data or results will be stored in the long term or shared with third parties. The data uploaded by the user in repDilPCR will be deleted as soon as the purpose of storage ceases to apply, i.e. with the expiration/termination of the session. The deletion of the data takes place immediately after the end of the session. Storage is only temporary for the time of processing. The user has access to the results within the Shiny App as well as the possibility to download results as CSV, PDF or PNG files. If the current session of the Shiny App is terminated (closing of the browser window, inactivity, timeout, ...), all uploaded and generated files of this session will be deleted."),
                    h2(""),
                    h4(tags$b("Provision of the website and creation of log files")),
                    h4(tags$i("Description and scope of data processing")),
                    h5("Each time the repDilPCR website is called up, our system automatically collects data and information from the computer system of the calling computer. The following data is transmitted in this process:",
                       tags$ul(tags$li("Information about the browser type and version used"),
                               tags$li("Date and time of access"),
                               tags$li("The Internet service provider of the user"))),
                    h5("The files uploaded by the user are processed as described above."),
                    h2(""),
                    h4(tags$i("Legal basis for the data processing")),
                    h5("The legal basis for the temporary storage of the data and the log files is Art. 6 (1) point f of the GDPR."),
                    h2(""),
                    h4(tags$i("Purpose of the data processing")),
                    h5("The temporary storage of the IP address by the system is necessary to enable delivery of the website to the user's computer. For this purpose, the user's IP address must remain stored for the duration of the session."),
                    h5("The storage in log files is done to ensure the functionality of the website. In addition, we use the data to ensure the security of our information technology systems. An evaluation of the data for marketing purposes does not take place."),
                    h5("These purposes are also covered by our legitimate interest in data processing according to Art. 6 para. 1 point f of the GDPR."),
                    h2(""),
                    h4(tags$i("Duration of storage")),
                    h5("The data is deleted as soon as it is no longer required to achieve the purpose for which it was collected. In respect to the collection of data for the provision of the website, this is the case when the respective session has ended."),
                    h2(""),
                    h4(tags$i("Possibility of objection and erasure")),
                    h5("The collection of data for the provision of the website and the storage of the data in log files are absolutely necessary for the operation of the website. Consequently, there is no possibility of objection on the part of the user."),
                    h2(""),
                    h4(tags$b("Use of cookies")),
                    h4(tags$i("Description and scope of data processing")),
                    h5("Our website does not use cookies."),
                    h2(""),
                    h4(tags$b("Personal settings and user account")),
                    h5("Personal accounts cannot be created on our website. The repDilPCR Shiny app is used without a user account via a session ID created in each case."),
                    h2(""),
                    h4(tags$b("E-mail contact")),
                    h4(tags$i("Description and scope of data processing")),
                    h5("It is possible to contact us via the e-mail address provided. In this case, the user's personal data transmitted with the e-mail will be stored within the e-mail conversation. The data will not be passed on to third parties and is used exclusively for processing the conversation."),
                    h2(""),
                    h4(tags$i("Purpose of data processing")),
                    h5("The case of contact by e-mail also constitutes the necessary legitimate interest in processing the data."),
                    h2(""),
                    h4(tags$i("Duration of storage")),
                    h5("The data are deleted as soon as they are no longer required to achieve the purpose for which they were collected. For personal data sent by e-mail, this is the case when the respective conversation with the user has ended. The conversation has ended when it can be inferred from the circumstances that the matter in question has been conclusively clarified."),
                    h2(""),
                    h4(tags$i("Possibility of objection and erasure")),
                    h5("If a user has engaged in an e-mail exchange with us, he/she can object to the storage of his/her personal data at any time. In such a case, the e-mail exchange cannot be continued and the data transmitted via e-mail will be deleted immediately."),
                    h5("You can revoke your consent at any time by sending an e-mail to ", a(href = "mailto:datenschutz@dkfz.de", "datenschutz@dkfz.de"), "."),
                    h5("In this case, all personal data stored in the course of the conversation with you will be deleted."),
                    h2(""),
                    h4(tags$b("Transfer of personal data to third parties")),
                    h5("No personal data will be passed on to third parties."),
                    hr(style = "height:4px;border:none;color:#333;background-color:#333;"),
                    h3(tags$b("Disclaimer und allgemeine Nutzungsbedingungen")),
                    h5("Die repDilPCR Shiny App ist ein Tool für wissenschaftliche Analyse, Visualisierung und Interpretation von qPCR Daten. Es liegt keine klinische Validierung vor, die Verwendung der durch die repDilPCR-App erzeugten Daten und Ergebnisse dient ausschließlich wissenschaftlichen Zwecken. Eine Nutzung zu diagnostischen oder therapeutischen Zwecken und/oder eine Implementierung der Ergebnisse in klinischen Studien unterliegt einzig und allein der Verantwortung des Nutzers/der Nutzerin (folgend genderneutral als „Nutzer“ definiert) und des behandelnden Arztes."),
                    h5("Für die Verarbeitung aller Daten innerhalb der repDilPCR Shiny App gilt:",
                       tags$ul(tags$li("Das Hochladen von Rohdaten erfolgt ausschließlich durch den Nutzer selbst, dieser allein ist für die technische Qualität (Durchführung der Labormethoden, d.h. Isolierung von Nukleinsäuren, reverse Transkription und Kultivierung und Polymerase-Kettenreaktion) und Richtigkeit der hochgeladenen Daten verantwortlich."),
                               tags$li("Die durch den Nutzer hochgeladenen Daten werden ausschließlich für die Dauer der Verarbeitung temporär gespeichert, es erfolgt keine dauerhafte Datenhaltung."),
                               tags$li("Die innerhalb der repDilPCR Shiny App erzeugten Ergebnisse können vom Nutzer entweder innerhalb der Shiny App angeschaut oder als CSV-, PDF- und/oder PNG-Dateien heruntergeladen werden."),
                               tags$li("Es werden keine Daten und keine Ergebnisse langfristig gespeichert oder weitergegeben. Diese stehen ausschließlich für den Nutzer als Download zur Verfügung. Wird die aktuelle Sitzung der Shiny App beendet (Schließen des Browserfensters, Inaktivität, Timeout, …), werden alle hochgeladenen und daraus erzeugte Dateien dieser Sitzung gelöscht."),
                               tags$li("Jeder Nutzer hat ausschließlich Zugriff und Einblick in die von ihm in der jeweiligen Sitzung (Session) hochgeladenen und prozessierten Daten."),
                               tags$li("Ein Zugriff auf Daten anderer Nutzer/Session IDs wird nicht erlaubt."),
                               tags$li("Die repDilPCR Shiny App ist ein rein wissenschaftliches Tool, welches allein der Prozessierung, Analyse und Visualisierung von Roh-oder präprozessierten Daten zu wissenschaftlichen Zwecken dient. Es liegt keine klinische Validierung vor."),
                               tags$li("Die weitere Nutzung der durch die repDilPCR Shiny App erzeugten Ergebnisse, sei es für Zwecke der Grundlagenforschung, Information für behandelnde Ärzte, Diagnostik, Therapieempfehlung oder Verwendung der Ergebnisse zum Zwecke von Studiendesigns, obliegt einzig und allein der Verantwortung des Nutzers und/oder des behandelten Arztes."))
                    ),
                    h3(tags$b("Impressum")),
                    h4(tags$b("Name und Anschrift")),
                    h5("Der Verantwortliche im Sinne der Datenschutz-Grundverordnung und anderer nationaler Datenschutzgesetze der Mitgliedsstaaten sowie sonstiger datenschutzrechtlicher Bestimmungen ist:"),
                    h5("Deutsches Krebsforschungszentrum - Stiftung des öffentlichen Rechts"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Deutschland"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:kontakt@dkfz.de", "kontakt@dkfz.de")),
                    h5("Website: ", a(href = "https://www.dkfz.de", "www.dkfz.de")),
                    h2(""),
                    h4(tags$b("Name und Anschrift des Datenschutzbeauftragten")),
                    h5("Datenschutzbeauftragter"),
                    h5("Deutsches Krebsforschungszentrum - Stiftung des öffentlichen Rechts"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Deutschland"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:datenschutz@dkfz.de", "datenschutz@dkfz.de")),
                    h2(""),
                    h4(tags$b("Rechtsform")),
                    h5("Das Deutsche Krebsforschungszentrum ist eine Stiftung des öffentlichen Rechts des Landes BadenWürttemberg. Gültig ist die Satzung vom Mai 2018. Umsatzsteuer-Identifikationsnummer: DE 143293537"),
                    h2(""),
                    h4(tags$b("Vertretungsberechtigte Personen")),
                    h5("Prof. Dr. Michael Baumann (Vorstandsvorsitzender und Wissenschaftlicher Vorstand)"),
                    h5("Ursula Weyrich (Kaufmännischer Vorstand)"),
                    h5("(Anschriften wie oben)"),
                    h2(""),
                    h4(tags$b("Inhaltlich verantwortlich")),
                    h5("Anschrift wie oben"),
                    h2(""),
                    h4(tags$b("Hinweise zum Urheberrecht")),
                    h5("Alle Inhalte der repDilPCR Shiny App sind unter der GNU General Public License v3 bzw. der GNU Affero General Public License v3 lizenziert."),
                    h2(""),
                    h4(tags$b("Haftungshinweis")),
                    h5("Trotz sorgfältiger inhaltlicher Kontrolle übernehmen wir keine Haftung für die Inhalte sowie für die prozessierten Ergebnisse. Die Nutzung der durch die repDilPCR Shiny App erzeugten Ergebnisse, sei es für Zwecke der Grundlagenforschung, Information für behandelnde Ärzte, Diagnostik, Therapieempfehlung oder Verwendung der Ergebnisse zum Zwecke von Studiendesigns oder Weiteres obliegt einzig und allein in der Verantwortung des Nutzers und/oder des behandelten Arztes."),
                    h2(""),
                    h4(tags$b("Geltungsbereich")),
                    h5('Diese Nutzungsbedingungen gelten für Ihre Nutzung von repDilPCR mit sämtlichen Inhalten, Funktionen, Diensten und hier aufgeführten Regeln. Wir behalten uns vor, auf repDilPCR zusätzliche Dienste anzubieten beziehungsweise repDilPCR regelmäßig zu aktualisieren. Diese Nutzungsbedingungen richten sich sowohl an Verbraucher als auch Unternehmen. Beide sind hier unter "Nutzer" zusammengefasst. Verbraucher sind natürliche Personen, die ein Rechtsgeschäft zu einem Zwecke abschließen, welches weder ihrer gewerblichen noch ihrer selbständigen beruflichen Tätigkeit zugerechnet werden kann. Unternehmer sind natürliche oder juristische Personen oder rechtsfähige Personengesellschaften, die bei Abschluss eines Rechtsgeschäfts in Ausübung ihrer gewerblichen oder selbständigen beruflichen Tätigkeit handeln.'),
                    h3(tags$b("Datenschutz")),
                    h4(tags$b("Name und Anschrift")),
                    h5("Der Verantwortliche im Sinne der Datenschutz-Grundverordnung und anderer nationaler Datenschutzgesetze der Mitgliedsstaaten sowie sonstiger datenschutzrechtlicher Bestimmungen ist:"),
                    h5("Deutsches Krebsforschungszentrum - Stiftung des öffentlichen Rechts"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Deutschland"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:kontakt@dkfz.de", "kontakt@dkfz.de")),
                    h5("Website: ", a(href = "https://www.dkfz.de", "www.dkfz.de")),
                    h2(""),
                    h4(tags$b("Name und Anschrift des Datenschutzbeauftragten")),
                    h5("Datenschutzbeauftragter"),
                    h5("Deutsches Krebsforschungszentrum - Stiftung des öffentlichen Rechts"),
                    h5("Im Neuenheimer Feld 280"),
                    h5("69120 Heidelberg"),
                    h5("Deutschland"),
                    h5("Tel.: +49 (0)6221 420"),
                    h5("E-Mail: ", a(href = "mailto:datenschutz@dkfz.de", "datenschutz@dkfz.de")),
                    h2(""),
                    h4(tags$b("Allgemeines zur Datenverarbeitung")),
                    h4(tags$i("Umfang der Verarbeitung personenbezogener Daten")),
                    h5("Es sind keine Angaben zum Nutzer notwendig. Beim Öffnen der repDilPCR Shiny App wird eine Sitzungskennung (Session ID) generiert, unter welcher das Hochladen der wissenschaftlichen Daten, das Prozessieren sowie das Herunterladen der generierten Ergebnisse erfolgt. Es findet keine Datenbündelung statt. Die hochgeladenen Daten werden nur temporär für die Dauer des Prozessierens gespeichert und nach Beenden der Sitzung gelöscht (siehe unten Datenlöschung). Es erfolgt keine langfristige Speicherung. Es werden keine Kundendaten oder personenbezogene Daten eingegeben oder gespeichert."),
                    h2(""),
                    h5(tags$i("Wissenschaftliche Daten")),
                    h5("Für Daten stammend von humanen Ausgangsmaterial gilt:"),
                    h5("Daten, die durch den Nutzer in die repDilPCR Shiny App hochgeladen werden, müssen vollständig anonymisiert oder pseudonymisiert sein. Die hochgeladenen Daten sind nur temporär im Zwischenspeicher zur Generierung der Resultate mit der Shiny App. Ein Personenbezug anhand der hochgeladenen Dokumente ist nicht erlaubt. Der Nutzer willigt ein, ausschließlich anonymisierte oder vollständig pseudonymisierte Daten in repDilPCR hochzuladen. Es kann kein Rückschluss auf den Ursprung der Daten gezogen werden."),
                    h5("Werden vom Nutzer Daten aus humanem Ausgangsmaterial verwendet, ist der Nutzer für die Einholung einer Patienten- oder Spender-Einwilligung verantwortlich. Diese muss die Verwendung des primären Materials und weiterer Analyse/Verarbeitung der Daten für wissenschaftliche Forschungszwecke erlauben."),
                    h5("Eine Ausnahme gilt in solchen Fällen, in denen eine vorherige Einholung einer Einwilligung aus tatsächlichen Gründen nicht möglich ist und die Verarbeitung der Daten durch gesetzliche Vorschriften gestattet ist. Hierfür ist allein der Nutzer der Shiny App verantwortlich. Der Nutzer willigt durch Bestätigung des Popup-Fensters beim Öffnen der Shiny App in diese Bedingungen ein."),
                    h5("Die repDilPCR Shiny App ist ein Tool ausschließlich für wissenschaftliche Analyse, Visualisierung und Interpretation der prozessierten Daten. Es liegt keine klinische Validierung vor, die Verwendung der durch die repDilPCR App erzeugten Daten und Ergebnisse dient ausschließlich wissenschaftlichen Zwecken. Eine Nutzung zu diagnostischen oder therapeutischen Zwecken und/oder eine Implementierung der Ergebnisse in klinischen Studien unterliegt einzig und allein in der Verantwortung des Nutzers und des behandelnden Arztes."),
                    h5("Der Nutzer wird beim Öffnen der Shiny App vor dem Hochladen der Files auf entsprechende Nutzungsregeln mittels Popup-Fenster hingewiesen, sowie auf diese Datenschutzerklärung hingewiesen."),
                    h5("Das Popup-Fenster beinhaltet den ausdrücklichen Hinweis, dass im Falle von Daten die aus humanen Ausgangsmaterial (Patientenmaterial) generiert wurden, diese ausschließlich anonymisiert oder pseudonymisiert hochgeladen und verarbeitet werden dürfen und eine Einwilligungserklärung des Patienten vorliegen muss, die die Verwendung des primären Materials zu Foschungszwecken und weitere Analyse/Verarbeitung der Daten zu Forschungszwecken erlaubt."),
                    h5("Der Nutzer akzeptiert diese Bedingungen sowie die Datenschutzerklärung durch Einwilligung in die AGBs durch Anklicken des entsprechenden Feldes innerhalb des Popup-Fensters. Erst dann kann ein Hochladen und weitere Verarbeitung der Rohdaten erfolgen."),
                    h2(""),
                    h4(tags$i("Rechtsgrundlage für die Verarbeitung personenbezogener Daten")),
                    h5("Soweit Daten prozessiert werden, die von humanen Ausgangsmaterial stammen ist der Nutzer für die Einholung der Einwilligungserklärung des Patienten/Spenders verantwortlich (siehe oben). Für die Verarbeitungsvorgänge personenbezogener Daten, dient Art. 6 Abs. 1 lit. a EU-Datenschutzgrundverordnung (DSGVO) als Rechtsgrundlage."),
                    h5("Der Nutzer willigt beim Öffnen der App ein, ausschließlich anonymisierte oder pseudonymisierte Daten zu prozessieren."),
                    h2(""),
                    h4(tags$i("Datenlöschung und Speicherdauer")),
                    h5("Es werden keine Daten und keine Ergebnisse langfristig gespeichert oder weitergegeben. Die durch den Nutzer in repDilPCR hochgeladenen Daten werden gelöscht, sobald der Zweck der Speicherung entfällt, mit Ablauf/Beendigung der Sitzung. Das Löschen der Daten erfolgt unmittelbar nach Sitzungsende. Eine Speicherung erfolgt nur temporär für die Zeit der Prozessierung. Der Nutzer hat innerhalb der Shiny App Einsicht auf die Ergebnisse sowie die Möglichkeit, Ergebnisse als CSV-, PDF- oder PNG-Dateien herunterzuladen. Wird die aktuelle Sitzung der Shiny App beendet (Schließen des Browserfensters, Inaktivität, Timeout, …), werden alle hochgeladenen und daraus erzeugte Dateien dieser Sitzung gelöscht."),
                    h2(""),
                    h4(tags$b("Bereitstellung der Website und Erstellung von Logfiles")),
                    h4(tags$i("Beschreibung und Umfang der Datenverarbeitung")),
                    h5("Bei jedem Aufruf der repDilPCR Website erfasst unser System automatisiert Daten und Informationen vom Computersystem des aufrufenden Rechners. Folgende Daten werden hierbei übermittelt:",
                       tags$ul(tags$li("Informationen über den Browsertyp und die verwendete Version"),
                               tags$li("Datum und Uhrzeit des Zugriffs"),
                               tags$li("Den Internet-Service-Provider des Nutzers"))),
                    h5("Die vom Nutzer hochgeladenen Dateien zur Prozessierung werden wie oben beschrieben verarbeitet."),
                    h2(""),
                    h4(tags$i("Rechtsgrundlage für die Datenverarbeitung")),
                    h5("Rechtsgrundlage für die vorübergehende Speicherung der Daten und der Logfiles ist Art. 6 Abs. 1 lit. f DSGVO."),
                    h2(""),
                    h4(tags$i("Zweck der Datenverarbeitung")),
                    h5("Die vorübergehende Speicherung der IP-Adresse durch das System ist notwendig, um eine Auslieferung der Website an den Rechner des Nutzers zu ermöglichen. Hierfür muss die IP-Adresse des Nutzers für die Dauer der Sitzung gespeichert bleiben."),
                    h5("Die Speicherung in Logfiles erfolgt, um die Funktionsfähigkeit der Website sicherzustellen. Zudem dienen uns die Daten zur Sicherstellung der Sicherheit unserer informationstechnischen Systeme. Eine Auswertung der Daten zu Marketingzwecken findet in diesem Zusammenhang nicht statt."),
                    h5("In diesen Zwecken liegt auch unser berechtigtes Interesse an der Datenverarbeitung nach Art. 6 Abs. 1 lit. f DSGVO."),
                    h2(""),
                    h4(tags$i("Dauer der Speicherung")),
                    h5("Die Daten werden gelöscht, sobald sie für die Erreichung des Zweckes ihrer Erhebung nicht mehr erforderlich sind. Im Falle der Erfassung der Daten zur Bereitstellung der Website ist dies der Fall, wenn die jeweilige Sitzung beendet ist."),
                    h2(""),
                    h4(tags$i("Widerspruchs- und Beseitigungsmöglichkeit")),
                    h5("Die Erfassung der Daten zur Bereitstellung der Website und die Speicherung der Daten in Logfiles ist für den Betrieb der Internetseite zwingend erforderlich. Es besteht folglich seitens des Nutzers keine Widerspruchsmöglichkeit."),
                    h2(""),
                    h4(tags$b("Verwendung von Cookies")),
                    h4(tags$i("Beschreibung und Umfang der Datenverarbeitung")),
                    h5("Unsere Webseite verwendet keine Cookies."),
                    h2(""),
                    h4(tags$b("Persönliche Einstellungen und angelegtes Profil")),
                    h5("Innerhalb unserer Website wird kein persönliches Profil angelegt. Die Nutzung der repDilPCR Shiny App erfolgt ohne Benutzerkonto über eine jeweils angelegte Sitzungskennung."),
                    h2(""),
                    h4(tags$b("E-Mail-Kontakt")),
                    h4(tags$i("Beschreibung und Umfang der Datenverarbeitung")),
                    h5("Eine Kontaktaufnahme über die bereitgestellte E-Mail-Adresse ist möglich. In diesem Fall werden die mit der E-Mail übermittelten personenbezogenen Daten des Nutzers innerhalb der E-Mail-Kontakt Konversation gespeichert. Es erfolgt in diesem Zusammenhang keine Weitergabe der Daten an Dritte. Die Daten werden ausschließlich für die Verarbeitung der Konversation verwendet."),
                    h2(""),
                    h4(tags$i("Zweck der Datenverarbeitung")),
                    h5("Im Falle einer Kontaktaufnahme per E-Mail liegt hieran auch das erforderliche berechtigte Interesse an der Verarbeitung der Daten."),
                    h2(""),
                    h4(tags$i("Dauer der Speicherung")),
                    h5("Die Daten werden gelöscht, sobald sie für die Erreichung des Zweckes ihrer Erhebung nicht mehr erforderlich sind. Für die personenbezogenen Daten die per E-Mail übersandt wurden, ist dies dann der Fall, wenn die jeweilige Konversation mit dem Nutzer beendet ist. Beendet ist die Konversation dann, wenn sich aus den Umständen entnehmen lässt, dass der betroffene Sachverhalt abschließend geklärt ist."),
                    h2(""),
                    h4(tags$i("Widerspruchs- und Beseitigungsmöglichkeit")),
                    h5("Nimmt der Nutzer per E-Mail Kontakt mit uns auf, so kann er der Speicherung seiner personenbezogenen Daten jederzeit widersprechen. In einem solchen Fall kann die Konversation nicht fortgeführt werden. Die über Email übermittelten Daten werden in diesem Fall unverzüglich gelöscht."),
                    h5("Sie können Ihre Einwilligung jederzeit per E-Mail an ", a(href = "mailto:datenschutz@dkfz.de", "datenschutz@dkfz.de"), " widerrufen."),
                    h5("Alle personenbezogenen Daten, die im Zuge der Kontaktaufnahme gespeichert wurden, werden in diesem Fall gelöscht."),
                    h2(""),
                    h4(tags$b("Weitergabe personenbezogener Daten an Dritte")),
                    h5("Es werden keine personenbezogenen Daten an Dritte weitergegeben.")
    ),
    size = "m", 
    closeOnEsc = FALSE,
    closeOnClickOutside = FALSE,
    html = TRUE,
    showConfirmButton = TRUE,
    showCancelButton = TRUE,
    confirmButtonText = "Accept / Zustimmen",
    confirmButtonCol = "#5eba7d",
    cancelButtonText = "Reject / Ablehnen",
    inputId = "accept",
    callbackR = function(x) {if(x == FALSE) shinyalert(
      text = "You have to accept the terms of use if you want to use the program.",
      size = "m",
      closeOnEsc = FALSE,
      closeOnClickOutside = FALSE,
      showConfirmButton = TRUE,
      confirmButtonText = "OK",
      callbackJS = "function(y) {if(y == true) {location.reload();}}"
    )}
  )
  

  ## Print warning messages if any
  warnings <- reactive({
    rd.warn(ref.sample = statistics.results()$ref.sample, rel.q.mean = statistics.results()$rel.q.mean, noref.warn = "A valid name of a sample to be used as baseline reference was not provided! Calculated relative quantities are not normalized to a particular sample.", statistics = input$statistics, posthoc = input$posthoc, nostatref.warn = "A reference group for the Dunnett post-hoc test has not been chosen. Please check your input.", frw = statistics.results()$frw, few.repl.warn = statistics.results()$few.repl.warn)
  })

  observeEvent(warnings()$noref.warn, {
    shinyalert("Warning!", warnings()$noref.warn, type = "warning")
  })

  observeEvent(warnings()$nostatref.warn, {
    shinyalert("Warning!", warnings()$nostatref.warn, type = "warning")
  })

  observeEvent(warnings()$few.repl.warn, {
    shinyalert("Warning!", warnings()$few.repl.warn, type = "warning")
  })

  ## Remove temporary files if any
  session$onSessionEnded(function() {
    delfil <- isolate(gsub(".{4}$", "*", input$input.table$name))
    if (length(Sys.glob(delfil)) != 0) {
      if (file.exists(Sys.glob(delfil))) {
        delfil.list <- Sys.glob(delfil)
        delfil.list <- delfil.list[!grepl(paste0(c(".R$", "Test_data.csv", "Test_data_precalc.csv", "counter.txt", "LICENSE", "README.md"), collapse = "|"), delfil.list)]
        file.remove(delfil.list)
      }
    }
  })

}

# Run app ----
shinyApp(ui = ui, server = server)