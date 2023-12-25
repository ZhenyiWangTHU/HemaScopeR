step1_fluidRow <- fluidRow(
      style = "margin-left: 10px;",
      column(
        3, align = "left", h3("Step 1. Input Data"),
        p("Please select the input data."),
        textInput("input.data.dirs", HTML("Enter data path and please use ';' to seperate different files<br>(e.g. /path1/file1/data1;/path2/file2/data2;path3/file3/data3):"), value = "NULL"),
        textInput("project.names", "Enter project name:", value = "project.names"),
        textInput("output.dir", "Enter output path:", value = ""),
        textInput("pythonPath", "Enter the path of Python:", value = "NULL"),
        selectInput("Step1_Input_Data.type", "Select Data Type:", choices = c("cellranger-count", "Seurat", "Matrix")),
        numericInput("gene.column", "Gene Column (default: 2):", value = 2),
        numericInput("min.cells", "Minimum Cells (default: 10):", value = 10),
        numericInput("min.feature", "Minimum Features (default: 200):", value = 200),
        textInput("mt.pattern", "Mt Pattern (default: '^MT-'):", value = '^MT-'),
        actionButton("load_data_button", "Load Data"),
        div(class = "spacer"),
        uiOutput("loadingData"),
        div(class = "spacer"),
        uiOutput("data_dim_output")
      )
    )