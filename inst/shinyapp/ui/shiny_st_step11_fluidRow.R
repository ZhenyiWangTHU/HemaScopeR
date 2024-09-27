step11_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 11. Generate the Report."),
        actionButton("RunStep11", "Run Step11"),
        div(class = "spacer"), 
        uiOutput("runningStep11"),
        div(class = "spacer"),  
        uiOutput("step11_completed"))
    )