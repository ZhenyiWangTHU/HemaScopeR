step8_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 8. Calculate Heterogeneity."),
        p("Please set the parameter needs for calculating heterogeneity."),
        textInput("ViolinPlot.cellTypeOrders", "Set the orders of cell types (seperate by ','):", value = "NULL"),
        actionButton("RunStep8", "Run Step8"),
        div(class = "spacer"), 
        uiOutput("runningStep8"),
        div(class = "spacer"),  
        uiOutput("step8_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 8. Calculate Heterogeneity:"),
            uiOutput("Step8.file_list"))
    )