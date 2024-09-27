step3_fluidRow <- fluidRow(
          style = "margin-left: 10px;",
          column(
            6, align = "left", h3("Step 3. Clustering"),
            p("Please input the parameters for clustering."),
            textInput("PCs.clustering", "PCs for clustering (default: 1:20):", value = "1:20"),
            numericInput("n.neighbors", "n.neighbors for clustering (default: 50):", value = 50),
            numericInput("resolution", "resolution for clustering (default: 0.4):", value = 0.4),
            actionButton("RunStep3", "Run Step3"),
            div(class = "spacer"),
            uiOutput("runningStep3"),
            div(class = "spacer"),
            uiOutput("step3_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 3. Clustering:"),
            uiOutput("Step3.file_list"))
        )