step5_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 5. Spatially variable features."),
        p("Please input the parameters for finding spatially variable features."),
        textInput("selection.method", "selection.method (default: 'moransi'):", value = 'moransi'),
        numericInput("n.top.show", "n.top.show (default: 10):", value = 10),
        numericInput("n.col.show", "n.col.show (default: 5):", value = 5),

        actionButton("RunStep5", "Run Step5"),
        div(class = "spacer"), 
        uiOutput("runningStep5"),
        div(class = "spacer"),  
        uiOutput("step5_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 5. Spatially variable features:"),
            uiOutput("Step5.st.file_list"))
    )