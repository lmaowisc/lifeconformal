#' Release call-local native owners before the package DLL is unloaded
#' @param libpath Installed package library path.
#' @noRd
.onUnload <- function(libpath) {
  library.dynam.unload("lifeconformal", libpath)
}
