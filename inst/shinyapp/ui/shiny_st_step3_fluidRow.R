step3_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 3. Clustering."),
        p("Please input the parameters for normalization, PCA and clustering."),
        textInput("normalization.method", "normalization.method (default: 'SCTransform'):", value = "SCTransform"),
        numericInput("npcs", "npcs (default: 50):", value = 50),
        textInput("pcs.used", "pcs.used (default: 1:10):", value = "1:10"),
        numericInput("resolution", "resolution (default: 0.8):", value = 0.8),

        actionButton("RunStep3", "Run Step3"),
        div(class = "spacer"), 
        uiOutput("runningStep3"),
        div(class = "spacer"),  
        uiOutput("step3_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 3. Clustering:"),
            uiOutput("Step3.st.file_list"))
    )