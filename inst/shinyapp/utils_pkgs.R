if (!requireNamespace("pacman")) {
  install.packages("pacman", repos = "http://cran.r-project.org")
}
library(pacman)

# if (!requireNamespace("gganatogram")) {
#   pacman::p_load(remotes)
#   tryCatch(
#     remotes::install_github("jespermaag/gganatogram"),
#     error = function(e) {
#       remotes::install_git("https://gitee.com/XenaShiny/gganatogram")
#     }
#   )
# }


# pacman::p_load(

# )

# if (!requireNamespace("ggradar")) {
#   pacman::p_load(remotes)
#   tryCatch(
#     remotes::install_github("ricardo-bion/ggradar"),
#     error = function(e) {
#       remotes::install_git("https://gitee.com/XenaShiny/ggradar")
#     }
#   )
# }

# if (packageVersion("UCSCXenaTools") < "1.4.4") {
#   tryCatch(
#     install.packages("UCSCXenaTools", repos = "http://cran.r-project.org"),
#     error = function(e) {
#       warning("UCSCXenaTools <1.4.4, this shiny has a known issue (the download button cannot be used) to work with it. Please upate this package!",
#         immediate. = TRUE
#       )
#     }
#   )
# }