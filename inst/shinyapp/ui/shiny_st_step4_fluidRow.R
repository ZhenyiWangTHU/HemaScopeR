step4_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 4. Differentially expressed genes."),
        p("Please input the parameters for finding differentially expressed genes in each cluster."),
        selectInput("only.pos", "only.pos:", choices = c("FALSE" = FALSE, "TRUE" = TRUE)),
        numericInput("min.pct", "min.pct (default: 0.25):", value = 0.25),
        numericInput("logfc.threshold", "logfc.threshold (default: 0.25):", value = 0.25),
        textInput("test.use", "test.use (default: 'wilcox'):", value = 'wilcox'),
        actionButton("RunStep4", "Run Step4"),
        div(class = "spacer"), 
        uiOutput("runningStep4"),
        div(class = "spacer"),  
        uiOutput("step4_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 4. Differentially expressed genes:"),
            uiOutput("Step4.st.file_list"))
    )