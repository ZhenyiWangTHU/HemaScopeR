step10_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 10. Niche analysis."),
        p("Please input the parameters for niche analysis (Please run step 8  deconvolution first)."),
        numericInput("Nich.cluster.n", "Nich.cluster.n (default: 4):", value = 4),
        actionButton("RunStep10", "Run Step10"),
        div(class = "spacer"), 
        uiOutput("runningStep10"),
        div(class = "spacer"),  
        uiOutput("step10_completed"))
    )