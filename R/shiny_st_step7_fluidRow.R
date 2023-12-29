step7_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 7. CNV analysis."),
        p("Please input the parameters for CNV analysis."),
        textInput("copykat.genome", "copykat.genome (default: 'NULL'):", value = 'NULL'),
        numericInput("copykat.LOW.DR", "copykat.LOW.DR (default: 0.05):", value = 0.05),
        numericInput("copykat.UP.DR", "copykat.UP.DR (default: 0.1):", value = 0.1),
        numericInput("copykat.win.size", "copykat.win.size (default: 25):", value = 25),
        textInput("copykat.distance", "copykat.distance (default: 'euclidean'):", value = 'euclidean'),
        numericInput("copykat.n.cores", "copykat.n.cores (default: 1):", value = 1),

        actionButton("RunStep7", "Run Step7"),
        div(class = "spacer"), 
        uiOutput("runningStep7"),
        div(class = "spacer"),  
        uiOutput("step7_completed"))
    )