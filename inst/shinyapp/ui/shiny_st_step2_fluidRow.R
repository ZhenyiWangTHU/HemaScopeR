step2_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 2. Quality Control."),
        p("Please input the parameters for quality control and preprocessing."),
        numericInput("min.gene", "min.gene (default: 200):", value = 200),
        numericInput("min.nUMI", "min.nUMI (default: 500):", value = 500),
        numericInput("max.gene", "max.gene (default: Inf):", value = Inf),
        numericInput("max.nUMI", "max.nUMI (default: Inf):", value = Inf),
        numericInput("min.spot", "min.spot (default: 0):", value = 0),
        selectInput("bool.remove.mito", "bool.remove.mito:", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
        selectInput("species", "species:", choices = c("mouse" = "mouse", "human" = "human")),

        actionButton("RunStep2", "Run Step2"),
        div(class = "spacer"), 
        uiOutput("runningStep2"),
        div(class = "spacer"),  
        uiOutput("step2_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 2. Quality Control:"),
            uiOutput("Step2.st.file_list"))
    )