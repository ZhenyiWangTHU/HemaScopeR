step11_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 11. GSVA."),
        p("Please input the parameters for GSVA."),
        selectInput("Step11_GSVA.identify.cellType.features", "Whether to identify cell type-specific GSVA terms:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
        selectInput("Step11_GSVA.identify.diff.features", "Whether to identify differential GSVA terms:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
        textInput("Step11_GSVA.comparison.design", 
                  HTML('The comparison design for GSVA<br>
                       e.g. list(
                            list(c("sample1","sample2"),
                                 c("sample3","sample4")),<br>
                            list(c("sample5","sample6"),
                                 c("sample7","sample8")))<br>
                            which indicates "sample1 and sample2 vs sample 3 and sample4",<br>
                             and "sample5 and sample6 vs sample 7 and sample8" :', 
                  value = "NULL")),
        actionButton("RunStep11", "Run Step11"),
        div(class = "spacer"), 
        uiOutput("runningStep11"),
        div(class = "spacer"),  
        uiOutput("step11_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 11. GSVA:"),
            uiOutput("Step11.file_list"))
    )