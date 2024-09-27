# options(shiny.autoreload = TRUE)
# Step1: Global setting ---------------------------------------------------------
options(shiny.fullstacktrace=TRUE)
HemaScopeR.runMode <- getOption("HemaScopeR.runMode", default = "client")
message("Run mode: ", HemaScopeR.runMode)
# 'client' for personal user, 'server' for running on server.

# Path for storing dataset files
if (is.null(getOption("HemaScopeR.cacheDir"))) {
  options(HemaScopeR.cacheDir = switch(HemaScopeR.runMode,
                                 client = file.path(tempdir(), "HemaScopeRShiny"), 
                                 server = "~/.HemaScopeRshiny"
  ))
}
HemaScopeR_DEST <- path.expand(file.path(getOption("HemaScopeR.cacheDir"), "datasets"))
if (!dir.exists(HemaScopeR_DEST)) {
  dir.create(HemaScopeR_DEST, recursive = TRUE)
}

options(shiny.maxRequestSize=1024*1024^2) #maximum file size upload

# Step2: Load necessary packages & data & functions ----------------------------------
message("Checking dependencies...")
source(system.file("shinyapp/utils_pkgs.R", package = "HemaScopeRShiny"))
source(system.file("shinyapp/utils_func.R", package = "HemaScopeRShiny"))
#source(system.file("shinyapp/utils_plot.R", package = "HemaScopeRShiny"))

# Put modules here --------------------------------------------------------
message("Loading modules and UIs...")
modules_path <- system.file("shinyapp", "modules", package = "HemaScopeRShiny", mustWork = TRUE)
modules_file <- dir(modules_path, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
sapply(modules_file, function(x, y) source(x, local = y), y = environment())

# Put page UIs here -----------------------------------------------------
pages_path <- system.file("shinyapp", "ui", package = "HemaScopeRShiny", mustWork = TRUE)
pages_file <- dir(pages_path, pattern = "\\.R$", full.names = TRUE, recursive = TRUE)
sapply(pages_file, function(x, y) source(x, local = y), y = environment())

# Step4: Run APP
# UI part ----------------------------------------------------------------------
message("Starting...")
ui <- tagList(
  tags$head(
    tags$title("HemaScopeRShiny"),
    tags$link(rel = "stylesheet", type = "text/css", href = "./css/global.css")
  ),

  useWaiter(), 
  waiterPreloader(html = tagList(
    spin_fading_circles(), 
    br(), br(),
    h1(strong("Welcome to use HemaScopeRShiny application!")),
    br(),
    p("A Specialized Bioinformatics Toolkit Designed for Analyzing both Single-cell and Spatial Transcriptome Sequencing Data from Hematopoietic Cells.",
      style = "font-size: 25px;"),
    br(),br(),
    p("Notes:", "(1) The initiation could take about 10 seconds. (2) Please zoom in or up screen for better representation.",
      style = "font-size: 16px;")
  ), color = "#2C3E50"),

  shinyjs::useShinyjs(),
  autoWaiter(html = spin_loader(), color = transparent(0.5)), # change style https://shiny.john-coene.com/waiter/

  navbarPage(
    id = "navbar",
    title = "HemaScopeRShiny",
    windowTitle = "HemaScopeRShiny",
    # inst/shinyapp/ui
    ui.page_home(),
    ui.page_scRNA_seq_pipeline(),
    ui.page_st_seq_pipeline(),
    ui.page_download(),
    ui.page_help(),
    ui.page_developers(),
    footer = ui.footer(),
    collapsible = TRUE,
    theme = tryCatch(shinythemes::shinytheme("flatly"),
                     error = function(e) {
                       "Theme 'flatly' is not available, use default."
                       NULL
                     })
  )
)

# Server Part ---------------------------------------------------------------
server <- function(input, output, session) {
  message("Shiny app run successfully! Enjoy it!\n")
  message("               --  HemaScopeR shiny team\n")
  # Stop warn
  storeWarn <- getOption("warn")
  options(warn = -1)
  # observe(print(input$navbar))
  # inst/shinyapp/server
  source(server_file("home.R"), local = TRUE)
  source(server_file("modules.R"), local = TRUE)
  source(server_file("general-analysis.R"), local = TRUE)
  observe_helpers(help_dir = system.file("shinyapp", "helper", package = "HemaScopeRShiny"))
}

# Run web app -------------------------------------------------------------
shiny::shinyApp(
  ui = ui,
  server = server
)