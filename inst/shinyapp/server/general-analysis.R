custom_file <- reactiveValues()

# # Upload custom feature/phenotype data
# observeEvent(input$ga_input_feature_file,{
#   req(input$ga_input_feature_file)
  
#   # ext <- tools::file_ext(input$ga_input_feature_file$name)
#   # shiny::validate(need(ext %in% c("csv","tsv", "gz"), "Please upload a csv/tsv file"))
  
#   inFile <- input$ga_input_feature_file
#   if (is.null(inFile))
#     return(NULL)
#   df <- data.table::fread(inFile$datapath, header = TRUE)
#   message("Saving custom feature data to temp directory.")
#   saveRDS(df, file = file.path(tempdir(), "custom_feature_data.rds"))
  
#   custom_file$fData <- df
# })

# # Individual analysis pages -------------------------------------------------------------
# observeEvent(req(input$navbar=="General Dataset Analysis"),{
#   callModule(
#     server.modules_ga_scatter_correlation, "module_ga_scatter_correlation",
#     selected_database_rm_phenotype, selected_database_add_url_and_phenotype,
#     custom_file
#   )
#   callModule(
#     server.modules_ga_matrix_correlation, "module_ga_matrix_correlation",
#     selected_database_rm_phenotype, selected_database_add_url_and_phenotype,
#     custom_file
#   )
#   callModule(
#     server.modules_ga_group_comparison, "module_ga_group_comparison",
#     selected_database_rm_phenotype, selected_database_add_url_and_phenotype,
#     custom_file
#   )
#   callModule(
#     server.modules_ga_surv_analysis, "module_ga_surv_analysis",
#     selected_database_rm_phenotype, selected_database_add_url_and_phenotype,
#     custom_file
#   )
#   callModule(
#     server.modules_ga_dim_distribution, "module_ga_dim_distribution",
#     selected_database_rm_phenotype, selected_database_add_url_and_phenotype,
#     custom_file
#   )
# }, once = TRUE)  

# # Show use alert ----------------------------------------------------------

# observeEvent(input$use_ga_page, {
#   shinyalert(
#     title = "General Analysis Usage",
#     text = paste(
#       "Firstly, select datasets from Repository page, the datasets and corresonding clinical datasets will be automatically loaded here (You can also upload your own data with the upload button).",
#       "Secondly, use any analysis feature below by clicking the tab.",
#       "Lastly, control how to analyze from left panel and filter samples from right panel. The result plot should be shown at the middle.",
#       sep = "\n\n"
#     ),
#     type = "info",
#     timer = 0,
#     confirmButtonCol = "#202324"
#   )
# })

# observeEvent(input$ga_drop_button, {
#   if (is.null(selected_database_add_url_and_phenotype())) {
#     sendSweetAlert(session,
#       title = "Warning!", text = "Please select datasets from Repository page firstly",
#       type = "warning"
#     )
#   }
# })