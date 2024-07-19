step4_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 4. Identify cell types automatically."),
        p("Please input the parameters for identifying cell types automatically."),
        selectInput("Org", "The organism for analysis (e.g., 'mmu' for mouse, 'hsa' for human):", choices = c("mmu" = "mmu", "hsa" = "hsa")),
        selectInput("Step4_Use_Which_Labels", "Step4_Use_Which_Labels:", choices = c('clustering' = 'clustering',
                                                                                     'abcCellmap.1' = 'abcCellmap.1',
                                                                                     'abcCellmap.2' = 'abcCellmap.2',
                                                                                     'abcCellmap.3' = 'abcCellmap.3',
                                                                                     'abcCellmap.4' = 'abcCellmap.4',
                                                                                     'HematoMap' = 'HematoMap')),
        selectInput("Step4_run_sc_CNV", "Step4_run_sc_CNV:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
        numericInput("ncores", "The number of CPU cores to use for parallel processing. (Default: 1):", value = 1),
        actionButton("RunStep4", "Run Step4"),
        div(class = "spacer"),
        uiOutput("runningStep4"),
        div(class = "spacer"),
        uiOutput("step4_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 4. Identify cell types automatically:"),
            uiOutput("Step4.file_list"))
    )