step6_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 6. Spatial interaction."),
        p("Please input the parameters for analyzing spatial interaction."),
        textInput("commot.signaling_type", "commot.signaling_type (default: 'Secreted Signaling'):", value = 'Secreted Signaling'),
        textInput("commot.database", "commot.database (default: 'CellChat'):", value = 'CellChat'),
        numericInput("commot.min_cell_pct", "commot.min_cell_pct (default: 0.05):", value = 0.05),
        numericInput("commot.dis_thr", "commot.dis_thr (default: 500):", value = 500),
        numericInput("commot.n_permutations", "commot.n_permutations (default: 100):", value = 100),

        actionButton("RunStep6", "Run Step6"),
        div(class = "spacer"), 
        uiOutput("runningStep6"),
        div(class = "spacer"),  
        uiOutput("step6_completed")),
      column(
            6, align = "left",
            h3("Browse files in Step 6. Spatial interaction:"),
            uiOutput("Step6.st.file_list"))
    )