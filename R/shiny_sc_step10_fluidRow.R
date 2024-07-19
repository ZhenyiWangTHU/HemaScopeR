step10_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 10. Calculate Lineage Scores."),
        p("Please input the parameters for calculating lineage scores."),
        textInput("lineage.genelist", HTML('The gene sets for calculating lineage scores.<br>Please enclose the gene sets with double quotation marks and seperate them by ";". <br>(e.g. "gene1,gene2,gene3";"gene4,gene5,gene6"):'), value = "NULL"),
        textInput("lineage.names", HTML('The names for the lineages. Please use "," to seperate them.<br>(e.g. lineage1,lineage2):'), value = "NULL"),
        textInput("groups_colors", HTML('The hexadecimal codes of colors for groups. Please use "," to seperate them.<br>(e.g. #FF0000,#0000FF):'), value = "NULL"),
        actionButton("RunStep10", "Run Step10"),
        div(class = "spacer"), 
        uiOutput("runningStep10"),
        div(class = "spacer"),  
        uiOutput("step10_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 10. Calculate Lineage Scores:"),
            uiOutput("Step10.file_list"))
    )