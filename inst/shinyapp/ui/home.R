home_text.list = list(
  intro = paste0(
    "Thank you for using HemaScopeRShiny",packageVersion("HemaScopeR")," based on ","HemaScopeR", packageVersion("HemaScopeR"),
    ". Our web tool aims to povide a user-friendly platform to explore both single-cell and spatial transcriptome sequencing data from hematopoietic cells, including myeloid and lymphoid lineages.",
    " If you have any questions during use, please do not hesitate to contact us via Github issues, Zhenyi Wang at wangzy17@tsinghua.org.cn, or Yuxin Miao at miaoyx21@mails.tsinghua.edu.cn.",
    " If the tool has faciliated your research, welcome to cite our work. :)"
  )
)


ui.page_home <- function() {
  tabPanel(
    title = "Home",
    icon = icon("home"), # create icon http://shiny.rstudio.com/reference/shiny/latest/icon.html
    fluidPage(
      useShinydashboard(),
      tags$head(
        tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap4.css")
      ),
      fluidRow(
        column(4,
          tags$div(
            column(
              12,
              tags$h2("Welcome to HemaScopeRShiny!"),
              tags$p(home_text.list$intro, style = "font-size: 19px;"),
              fluidRow(
                column(4, align = "center",
                  actionBttn("bt01","Github",
                    style = "bordered", color = "primary", icon = icon("github"),
                    onclick=paste0("window.open('https://github.com/ZhenyiWangTHU/HemaScopeR','_blank')"))
                ),
                column(4, align = "center",
                  actionBttn("bt02","Tutorial",
                    style = "bordered", color = "primary", icon = icon("book"),
                    onclick=paste0("window.open('https://ZhenyiWangTHU.github.io/HemaScopeRShiny_Book','_blank')"))
                ),
                column(4, align = "center",
                  actionBttn("bt03","Article",
                    style = "bordered", color = "primary", icon = icon("newspaper"))
                ),
              ),
              uiOutput("citation"),
              tags$hr(),
            )
          )
        ),
        column(4,
          wellPanel(
            style = "background: #b3cde3",
            h3("Daily papers", icon = "dice", align = "center"),
            ui.home_daily_papers("homepage_daily_papers"),
          )
        ),
        column(4,
          wellPanel(
            style = "background: #a6cee3",
            h3("Marker Genes Query", align = "center"),
            # br(),
            ui.home_Marker_Genes_Query_box("homepage_marker_genes_query_box"),
          )
        ),
      ),
      tags$br(),
      h2(strong("※ Shiny Page Gallery")),
      tags$hr(style = "border:none; border-top:5px solid #5E81AC;"),
      fluidRow(
        column(
          12,
          slickR::slickROutput("slick_output", width = "90%", height = "700px")
        )
      ),
      br(),
      h2(strong("※ Latest significant release notes")),
      tags$hr(style = "border:none; border-top:5px solid #5E81AC;"),
      fluidRow(
        column(12, #offset = 1,
          tags$ul(
            tags$li("2024-05-1: During the major version upgrade of HemaScopeR.",style = "font-size: 20px;"),
            tags$li("2023-12-30: Open source HemaScopeR on GitHub.",style = "font-size: 20px;"),
            tags$li("See more update logs in our", 
              a("Github", href="https://github.com/ZhenyiWangTHU/HemaScopeR/"), ".",
              style = "font-size: 20px;")
          )  
        )
      )
    )
  )
}