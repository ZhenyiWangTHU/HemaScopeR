step12_fluidRow <- fluidRow(
          style = "margin-left: 10px;",
          column(
            6, align = "left", h3("Step 12. Construct Trajectories."),
            textInput("Step12_Construct_Trajectories.clusters", "Set the cell types for constructing trajectories(seperate by ',' and default is 'all'):", value = "all"),
            # monocle
            p("Set the parameters for monocle2."),
            selectInput("Step12_Construct_Trajectories.monocle", "Whether to run monocle2:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
            # slingshot
            p("Set the parameters for slingshot."),
            selectInput("Step12_Construct_Trajectories.slingshot", "Whether to run slingshot:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
            textInput("slingshot.start.clus", "Set the root clusters:", value = "NULL"),
            textInput("slingshot.end.clus", "Set the tip clusters:", value = "NULL"),
            textInput("slingshot.colors", "Set the hexadecimal codes of colors for clusters (seperate by ','):", value = "NULL"),
            # scVelo
            p("Set the parameters for scVelo."),
            selectInput("Step12_Construct_Trajectories.scVelo", "Whether to run scVelo:", choices = c("TRUE" = TRUE, "FALSE" = FALSE)),
            textInput("loom.files.path", "Enter the paths of loom files (seperate by ';'):", value = "NULL"),
            actionButton("RunStep12", "Run Step12"),
            div(class = "spacer"),
            uiOutput("runningStep12"),
            div(class = "spacer"),
            uiOutput("Step12_completed")),
            column(
                    6, align = "left",
                    h3("Browse files in Step 12. Construct Trajectories:"),
                    uiOutput("Step12.file_list"))
        )