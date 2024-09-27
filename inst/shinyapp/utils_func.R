# Obtain path to individual server code parts ----------------------------
server_file <- function(x) {
  server_path <- system.file("shinyapp", "server",
    package = "HemaScopeRShiny", mustWork = TRUE
  )
  file.path(server_path, x)
}


# Set utility functions ---------------------------------------------------
