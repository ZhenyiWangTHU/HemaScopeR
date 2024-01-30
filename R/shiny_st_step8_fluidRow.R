step8_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 8. Deconvolution."),
        p("Please input the parameters for deconvolution."),
        textInput("cell2loc.sc.h5ad.dir", "cell2loc.sc.h5ad.dir (default: 'NULL'):", value = 'NULL'),
        numericInput("cell2loc.sc.max.epoch", "cell2loc.sc.max.epoch (default: 1000):", value = 1000),
        numericInput("cell2loc.st.max.epoch", "cell2loc.st.max.epoch (default: 10000):", value = 10000),
        selectInput("cell2loc.use.gpu", "cell2loc.use.gpu (default: FALSE):", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),

        actionButton("RunStep8", "Run Step8"),
        div(class = "spacer"), 
        uiOutput("runningStep8"),
        div(class = "spacer"),  
        uiOutput("step8_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 8. Deconvolution:"),
            uiOutput("Step8.st.file_list"))
    )