step5_fluidRow <- fluidRow(
          style = "margin-left: 10px;",
          column(
            6, align = "left", h3("Step 5. Visualization."),
            p("Please input the parameters for visualization."),
            numericInput("phate.knn", "The number of nearest neighbors for PhateR analysis(default: 50):", value = 50),
            numericInput("phate.npca", "The number of principal components for PhateR analysis(default: 20):", value = 20),
            numericInput("phate.t", "The t parameter for PhateR analysis(default: 10):", value = 10),
            numericInput("phate.ndim", "The number of dimensions for PhateR analysis(default: 2):", value = 2),  
            actionButton("RunStep5", "Run Step5"),
            div(class = "spacer"), 
            uiOutput("runningStep5"),
            div(class = "spacer"),  
            uiOutput("step5_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 5. Visualization:"),
            uiOutput("Step5.file_list"))
        )