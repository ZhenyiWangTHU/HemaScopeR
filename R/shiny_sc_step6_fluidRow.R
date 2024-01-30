step6_fluidRow <- fluidRow(
          style = "margin-left: 10px;",
          column(
            6, align = "left", h3("Step 6. Find DEGs."),
            p("Please input the parameters for finding DEGs."),

            numericInput("min.pct", "The minimum percentage for a gene to be considered in differential gene detection(default: 0.25):", value = 0.25),
            numericInput("logfc.threshold", "The log-fold change threshold for differential gene detection(Default: 0.25):", value = 0.25),  
              
            actionButton("RunStep6", "Run Step6"),
            div(class = "spacer"), 
            uiOutput("runningStep6"),
            div(class = "spacer"),  
            uiOutput("step6_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 6. Find DEGs:"),
            uiOutput("Step6.file_list"))
        )