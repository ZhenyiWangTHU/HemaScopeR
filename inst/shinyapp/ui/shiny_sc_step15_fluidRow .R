step15_fluidRow <- fluidRow(
          style = "margin-left: 10px;",
          column(
            6, align = "left", h3("Step 15. Generate the Report."),
            p("Please input the parameters for generating the report."),

            actionButton("RunStep15", "Run Step15"),
            div(class = "spacer"), 
            uiOutput("runningStep15"),
            div(class = "spacer"),  
            uiOutput("step15_completed"))
        )