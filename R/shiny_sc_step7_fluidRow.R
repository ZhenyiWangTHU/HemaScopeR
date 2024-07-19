step7_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 7. Assign Cell Cycles."),
        p("Please input the parameters for assigning cell cycles."),
        textInput("cellcycleCutoff", "Set the cutoff for cell cycle scores(default: NULL):", value = "NULL"),
        actionButton("RunStep7", "Run Step7"),
        div(class = "spacer"), 
        uiOutput("runningStep7"),
        div(class = "spacer"),  
        uiOutput("step7_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 7. Assign Cell Cycles:"),
            uiOutput("Step7.file_list"))
    )