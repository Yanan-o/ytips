#' @title Find matching paths in OneDrive
#' @description Searches detected OneDrive folders, not the full computer, for
#' matching paths.
#'
#' @param x File or folder name, or a path suffix using `/`.
#' @param contents Whether matching folders return their first-level contents.
#' @return All matching paths, or folder contents when `contents = TRUE`.
#'
#' @examples
#' \dontrun{
#' onedrive("Y.tips.R")
#' onedrive("Y.R")
#' onedrive("Y.R", contents = TRUE)
#' onedrive("Y.R/Y.tips.R")
#' }
#'
#' @family utilities
#' @export
#'

onedrive <- function(x, contents = FALSE){
  if(missing(x) || !is.character(x) || length(x) != 1L || !nzchar(x)){stop("x must be one file or folder name")}
  yy_check_pkg("fs")

  roots <- c(Sys.getenv(c("OneDrive", "OneDriveConsumer", "OneDriveCommercial"), unset = ""), Sys.glob(path.expand("~/OneDrive*")),
             Sys.glob(path.expand("~/Library/CloudStorage/OneDrive*")))
  roots <- roots[nzchar(roots) & fs::dir_exists(roots)]
  roots <- unique(normalizePath(roots, winslash = "/", mustWork = TRUE)); roots <- roots[grepl("(^|/)onedrive", tolower(roots))]
  character() -> res
  for(i in roots){
    if(fs::dir_exists(i)){c(res, as.character(fs::dir_ls(i, recurse = TRUE, glob = paste0("*", x), fail = FALSE))) -> res}
  }
  res <- unique(res)
  if(length(res) == 0 || !isTRUE(contents)){return(res)}
  unlist(lapply(res, function(i){
    if(fs::dir_exists(i)){as.character(fs::dir_ls(i, type = "any", fail = FALSE))} else{i}
  }), use.names = FALSE) %>% unique()
}
