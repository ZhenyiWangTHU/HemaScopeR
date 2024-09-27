step1_fluidRow_st <- fluidRow(
      style = "margin-left: 10px;",
      column(
        3, align = "left", h3("Step 1. Input Data"),
        p("Please select the input data."),
        textInput("input.data.dir", HTML("Enter data path:"), value = "NULL"),
        textInput("sampleName", "Enter sample name:", value = "Hema_ST"),
        textInput("output.dir", "Enter output path:", value = ""),
        textInput("pythonPath", "Enter the path of Python:", value = "NULL"),
        actionButton("load_data_button", "Load Data"),
        div(class = "spacer"),
        uiOutput("loadingData"),
        uiOutput("step1_completed"))
    )