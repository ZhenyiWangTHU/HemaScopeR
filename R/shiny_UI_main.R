# ui--------------------------------------------------------------------------------------------------------------------------------------------
ui <- fluidPage(
  shinyjs::useShinyjs(),  # shinyjs

  # ui1
  div(id = "ui1", style = "display: flex; flex-direction: column; align-items: center; justify-content: center; height: 70vh;",
      fluidRow(),
      fluidRow(
        column(3, align = "center", imageOutput('logo'))
      ),
      fluidRow(
        column(12, align = "center", h1("HemaScopeR: A Specialized Bioinformatics Toolkit Designed for Analyzing both Single-cell and Spatial Transcriptome Sequencing Data from Hematopoietic Cells", class = "h1-font"))
      ),
      fluidRow(div(class = "spacer")),  # empty line
      fluidRow(div(class = "spacer")),  # empty line
      fluidRow(
        column(4, align = "center", actionButton("start_button", "Start scRNA-seq pipeline", style = "font-family: 'arial'; font-size: 18pt;"))
      ),
      fluidRow(div(class = "spacer")),  # empty line
      fluidRow(
        column(4, align = "center", actionButton("start_button_st", "Start st-seq pipeline", style = "font-family: 'arial'; font-size: 18pt;"))
      ),
      uiOutput("ui_styles")
  ),

  # ui2
  div(id = "ui2", style = "display: none;",
      titlePanel("scRNA-seq pipeline"),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          style = "width: 45%",
          fluidRow(
            column(12, align = "left",
                   actionButton("back_button", "Back to Home", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step1", "Step 1. Input Data", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step2", "Step 2. Quality Control", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step3", "Step 3. Clustering", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step4", "Step 4. Identify Cell Types", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step5", "Step 5. Visualization", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step6", "Step 6. Find DEGs", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step7", "Step 7. Assign Cell Cycles", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step8", "Step 8. Calculate Heterogeneity", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step9", "Step 9. Violin Plot for Marker Genes", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step10", "Step 10. Calculate Lineage Scores", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step11", "Step 11. GSVA", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step12", "Step 12. Construct Trajectories", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step13", "Step 13. TF Analysis", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step14", "Step 14. Cell-Cell Interaction", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step15", "Step 15. Generate the Report", style = "width: 100%; text-align: left;"))
          )
        ),
        mainPanel = mainPanel(
          style = "width: 55%; margin-left: -400px; margin-top: -10px;",
          uiOutput("stepContent")
       )
      )
     ),

  # ui3
  div(id = "ui3", style = "display: none;",
      titlePanel("st-seq pipeline"),
      sidebarLayout(
        sidebarPanel = sidebarPanel(
          style = "width: 45%",
          fluidRow(
            column(12, align = "left",
                   actionButton("back_button_st", "Back to Home", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step1_st", "Step 1. Input Data", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step2_st", "Step 2. Quality Control", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step3_st", "Step 3. Clustering", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step4_st", "Step 4. Find Differentially Expressed Genes", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step5_st", "Step 5. Spatially Variable Features", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step6_st", "Step 6. Spatial Interaction", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step7_st", "Step 7. CNV Analysis", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step8_st", "Step 8. Deconvolution", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step9_st", "Step 9. Cell Cycle Analysis", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step10_st", "Step 10. Niche Analysis", style = "width: 100%; text-align: left;")),
            column(12, align = "left",
                   actionButton("step11_st", "Step 11. Generate the Report", style = "width: 100%; text-align: left;"))
          )
        ),
        mainPanel = mainPanel(
          style = "width: 55%; margin-left: -400px; margin-top: -10px;",
          uiOutput("stepContent_st")
       )
      )
     )
)