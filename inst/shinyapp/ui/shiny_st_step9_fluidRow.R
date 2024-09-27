step9_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 9. Cell cycle analysis."),
        p("Please input the parameters for cell cycle analysis."),
        textInput("s.features", HTML('The gene sets for calculating S phase scores(e.g. "gene1,gene2,gene3"):'), value = "NULL"),
        textInput("g2m.features", HTML('The gene sets for calculating G2M phase scores(e.g. "gene1,gene2,gene3"):'), value = "NULL"),
        actionButton("RunStep9", "Run Step9"),
        div(class = "spacer"), 
        uiOutput("runningStep9"),
        div(class = "spacer"),  
        uiOutput("step9_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 9. Cell cycle analysis:"),
            uiOutput("Step9.st.file_list"))
    )