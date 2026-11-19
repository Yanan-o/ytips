#' @title ytips: reusable R helpers for ecology and data analysis
#' @description `ytips` contains plotting, model-summary, piecewiseSEM, ordination and
#' co-occurrence-network helpers migrated from the author's analysis scripts.
#'
#' @keywords internal
#' @importFrom ggplot2 aes annotate arrow coord_flip coord_radial
#'   element_blank element_rect element_text expansion geom_bar geom_boxplot
#'   geom_col geom_contour geom_errorbar geom_errorbarh geom_hline geom_jitter
#'   geom_line geom_path geom_point geom_pointrange geom_rect geom_ribbon
#'   geom_segment geom_smooth geom_text geom_tile geom_violin geom_vline ggplot
#'   guide_colorbar guide_legend guides labs margin position_dodge
#'   position_nudge position_stack scale_color_manual scale_colour_manual
#'   scale_fill_gradientn scale_fill_manual scale_linetype_manual
#'   scale_shape_manual scale_y_discrete theme theme_bw theme_classic
#'   theme_minimal theme_void xlim ylim
#' @importFrom dplyr across all_of arrange bind_cols bind_rows contains count
#'   dense_rank desc distinct everything group_by if_else matches mutate pull
#'   recode reframe rename rename_with select slice slice_head summarise ungroup
#'   where
#' @importFrom purrr imap map set_names
#' @importFrom tibble as_tibble column_to_rownames remove_rownames
#'   rownames_to_column
#' @importFrom tidyr pivot_longer separate
#' @importFrom stringr str_remove_all str_replace str_replace_all str_split
#' @importFrom lme4 fixef getData VarCorr
#' @importFrom lmerTest lmer
#' @importFrom stats anova as.formula density formula getCall lm median qqline
#'   qqnorm quantile sd setNames shapiro.test var
#' @importFrom graphics abline barplot boxplot hist legend par title
#' @importFrom utils capture.output combn head
#' @importFrom grDevices colorRampPalette rgb
#' @importFrom grid unit
#' @importFrom scales alpha
#' @importFrom magrittr %>%
"_PACKAGE"

.onAttach <- function(libname, pkgname){
  if(nzchar(Sys.getenv("_R_CHECK_PACKAGE_NAME_"))){return(invisible())}
  packages <- c("openxlsx2", "here", "ggplot2", "dplyr", "tibble", "tidyr", "stringr", "purrr", "lmerTest")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if(length(missing) > 0){
    packageStartupMessage("Installing ytips dependencies: ", paste(missing, collapse = ", "))
    utils::install.packages(missing, dependencies = TRUE)
  }
  for(i in packages){
    if(requireNamespace(i, quietly = TRUE) && !paste0("package:", i) %in% search()){suppressPackageStartupMessages(base::attachNamespace(i))}
  }
  invisible()
}


yy_pkg_available <- function(pkg){
  requireNamespace(pkg, quietly = TRUE)
}



#' @title Check and report missing packages
#' @description Checks whether required packages are installed. In an
#' interactive R session, missing packages trigger a prompt to install them or
#' stop the current call. Non-interactive sessions stop with the installation
#' command.
#'
#' @param pkg Package name or character vector of package names.
#' @param install.hint Optional installation command.
#' @param ask Whether to ask before installing. Defaults to `interactive()`.
#' @return Invisibly returns `TRUE` when all packages are available.
#' @family utilities
#' @keywords internal
#'

yy_check_pkg <- function(pkg, install.hint = NULL, ask = interactive()){
  missing <- pkg[!vapply(pkg, yy_pkg_available, logical(1))]
  if(length(missing) > 0){
    hint <- if(is.null(install.hint)){
      paste0('install.packages(c(', paste0('"', missing, '"', collapse = ', '), '))')
    } else{
      install.hint
    }
    if(isTRUE(ask)){
      choice <- utils::menu(
        c("Install now", "Do not install"),
        title = paste0("Package(s) not installed: ", paste0(missing, collapse = ", "))
      )
      if(choice == 1L){
        message("Installing package dependencies with:\n", hint)
        tryCatch(
          eval(parse(text = hint), envir = parent.frame()),
          error = function(e){
            rlang::abort(c(
              "Package installation failed.",
              x = conditionMessage(e),
              i = paste0("Install manually with: ", hint)
            ))
          }
        )
        missing <- pkg[!vapply(pkg, yy_pkg_available, logical(1))]
        if(length(missing) == 0){return(invisible(TRUE))}
        rlang::abort(c(
          "Package installation did not make all required packages available.",
          x = paste0("Still missing: ", paste0(missing, collapse = ", ")),
          i = paste0("Install manually with: ", hint)
        ))
      }
    }
    rlang::abort(c(
      paste0("Package(s) required but not installed: ", paste0(missing, collapse = ", ")),
      i = paste0("Install with: ", hint)
    ))
  }
  invisible(TRUE)
}

yy_check_phyloseq <- function(){
  yy_check_pkg("phyloseq", 'pak::pak("bioc::phyloseq")')
}

yy_check_network_dependencies <- function(spiec.easi = FALSE, taxa.shuffle = FALSE){
  yy_check_phyloseq()
  yy_check_pkg("igraph")
  if(isTRUE(spiec.easi)){
    yy_check_pkg("SpiecEasi", 'pak::pak("zdk123/SpiecEasi")')
  } else if(!isTRUE(taxa.shuffle)){
    yy_check_pkg("Hmisc")
  }
  invisible(TRUE)
}

yy_check_sem_dependencies <- function(){
  yy_check_pkg(
    c("piecewiseSEM", "DiagrammeR"),
    'pak::pak(c("piecewiseSEM", "DiagrammeR"))'
  )
}

yy_check_random_forest_dependencies <- function(){
  yy_check_pkg(
    c("rfPermute", "randomForest"),
    'install.packages(c("rfPermute", "randomForest"))'
  )
}
