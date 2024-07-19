step14_fluidRow <- fluidRow(
          style = "margin-left: 10px;",
          column(
            6, align = "left", h3("Step 14. Cell-Cell Interaction."),
            p("Please input the parameters for cell-cell interaction analysis."),
            selectInput("sorting", "The cell groups were sorted?:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
            actionButton("RunStep14", "Run Step14"),
            div(class = "spacer"), 
            uiOutput("runningStep14"),
            div(class = "spacer"),  
            uiOutput("step14_completed")),
          column(
            6, align = "left",
            h3("Browse files in Step 14. Cell-Cell Interection:"),
            uiOutput("Step14.file_list"))
          )
