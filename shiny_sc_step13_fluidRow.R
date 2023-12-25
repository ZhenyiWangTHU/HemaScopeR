step13_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        6, align = "left", h3("Step 13. TF Analysis."),
        p("Please input the parameters for TF analysis."),
        textInput("Step13_TF_Analysis.cellTypes_colors", "Set the hexadecimal codes of colors for cell types (seperate by ','):", value = "NULL"),
        textInput("Step13_TF_Analysis.groups_colors", "Set the hexadecimal codes of colors for groups (seperate by ','):", value = "NULL"),
        actionButton("RunStep13", "Run Step13"),
        div(class = "spacer"), 
        uiOutput("runningStep13"),
        div(class = "spacer"),  
        uiOutput("step13_completed")),
      column(
          6, align = "left",
          div(
            style = "border: 2px solid #ddd; padding: 10px; text-align: center; min-height: 300px; min-width: 400px;",
            conditionalPanel(
              condition = "!output.vln_plot",
              div(
                h4("Graphical Display Area"),
                p("Graph will be displayed here after clicking 'Run'.")
              )
            ),
            plotOutput("vln_plot")  # 创建一个输出区域来显示 VlnPlot
          ),
          fluidRow(  
            column(6,
              textInput("color", "Change Vln Color (Hexadecimal):", value = "#FF5733"),
              actionButton("edit_button", "Edit the Figure")
            )
          ), 
          fluidRow(
            column(
              6, # 控制PDF Width输入框的布局
              numericInput("pdf_width", "PDF Width (default: 8):", value = 8)
            ),        
            column(
              6, # 控制PNG Width输入框的布局
              numericInput("png_width", "PNG Width (default: 800):", value = 800)
            )
          ),
          fluidRow(  
            column(
              6, # 控制PDF Height输入框的布局
              numericInput("pdf_height", "PDF Height (default: 6):", value = 6)
            ),
            column(
              6, # 控制PNG Height输入框的布局
              numericInput("png_height", "PNG Height (default: 600):", value = 600)
            )
          ),
          fluidRow(  
            column(
              6, # Download PDF按钮的布局
              downloadButton("download_pdf_button", "Download PDF")
            ),
            column(
              6, # Download PNG按钮的布局
              downloadButton("download_png_button", "Download PNG")
            )
          )
        )
    )