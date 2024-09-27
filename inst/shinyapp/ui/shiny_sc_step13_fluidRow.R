step13_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 13. TF Analysis."),
        p("Please input the parameters for TF analysis."),
        textInput("Step13_TF_Analysis.cellTypes_colors", "Set the hexadecimal codes of colors for cell types (seperate by ','):", value = "NULL"),
        textInput("Step13_TF_Analysis.groups_colors", "Set the hexadecimal codes of colors for groups (seperate by ','):", value = "NULL"),
        actionButton("RunStep13", "Run Step13"),
        div(class = "spacer"), 
        uiOutput("runningStep13"),
        div(class = "spacer"),  
        uiOutput("step13_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 13. TF Analysis:"),
            uiOutput("Step13.file_list"))
    )