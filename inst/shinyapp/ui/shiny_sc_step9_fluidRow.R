step9_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 9. Violin Plot for Marker Genes."),
        p("Set parameters for violin plot of marker genes."),
        textInput("marker.genes", "Enter marker genes for violin plot (seperate by ','):", value = "NULL"),
        textInput("ViolinPlot.cellTypeColors", "Set the hexadecimal codes of colors for cell types (seperate by ','):", value = "NULL"),
        actionButton("RunStep9", "Run Step9"),
        div(class = "spacer"), 
        uiOutput("runningStep9"),
        div(class = "spacer"),  
        uiOutput("step9_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 9. Violin Plot for Marker Genes:"),
            uiOutput("Step9.file_list"))
    )