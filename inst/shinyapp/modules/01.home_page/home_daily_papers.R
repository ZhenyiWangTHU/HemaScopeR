ui.home_daily_papers <- function(id) {
  ns <- NS(id)
  tagList(
    # h3("Daily papers", icon = "dice", align = "center"),
    fluidRow(
      column(12, align = "center",
        uiOutput(ns("papers_latest")),

      ),
    ),
  )
}

server.home_daily_papers <- function(input, output, session) {
  ns <- session$ns
  papers_latest = reactiveValues(papers="")

  papers_latest = reactive({

  })

  output$papers_papers_latest = renderUI({

  })
}

