#' @title Color palette helper
#' @description Generates a color vector from one of three preset palettes.
#'
#' @param n Number of colors wanted, or an index vector.
#' @param colbar Palette identifier: 1, 2, or 3.
#' @return A character vector of colors.
#' @examples
#' color_y(4)
#' color_y(4:5, colbar = 2)
#' @family plot-tools
#' @export
#'

color_y <- function(n, colbar = 3){
  if(colbar == 1){c("black", "steelblue2", "#FF1493") -> cols}
  if(colbar == 2){c("#A9D179", "#F5AE6B", "#4387B5", "#CC79A7", "#BCB8D3", "#84CAC0", "#fc9272") -> cols}
  if(colbar == 3){c("#15983DFF", "#149BEDFF", "#EC579AFF", "#FA6B09FF", "#8F2F8BFF", "#EE0011FF", "#F9B90AFF") -> cols}
  if(length(n) == 1){if(n > 7){grDevices::colorRampPalette(c("steelblue2", "white", "#FF1493"))(n) -> cols};
    if(n == 1){"black" -> cols}; if(n == 2){cols[c(2, 3)] -> cols}; if(n > 2 & n < 8){cols[1:n] -> cols}
  } else{cols[n] -> cols}; return(cols) }



#' @title Default clean ggplot theme
#' @description Creates the package's default white-background theme.
#'
#' @param base.size Base font size. Defaults to 11.
#' @return A ggplot2 theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_y()
#' @family plot-tools
#' @export
#'

theme_y <- function(base.size = 11){
  ggplot2::theme_bw(base_size = base.size) +
  ggplot2::theme(panel.grid = ggplot2::element_blank()) + ggplot2::theme(legend.background = ggplot2::element_blank()) +
  ggplot2::theme(axis.text = ggplot2::element_text(colour = 1), plot.title = ggplot2::element_text(hjust = .5, face = "bold", size = 9))
  }


no_legend <- ggplot2::theme(legend.position = "none")
no_axis_y <- ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank(), axis.title.y = ggplot2::element_blank())
no_axis_x <- ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank())
no_axis_y_title <- ggplot2::theme(axis.title.y = ggplot2::element_blank())
no_axis_x_title <- ggplot2::theme(axis.title.x = ggplot2::element_blank())



#' @title Save or print a ggplot or patchwork figure
#' @description Saves or prints ggplot and patchwork figures.
#'
#' @param p A ggplot object, or a list of ggplot objects.
#' @param file Output path without the `.tiff` extension.
#' @param nrow Number of rows.
#' @param width Figure width in inches.
#' @param height Figure height in inches.
#' @param spacer Integer positions at which to insert patchwork spacers.
#' @param xlim,ylim Optional shared axis limits.
#' @param x.axis,y.axis Whether to retain x/y axis text and ticks.
#' @param x.lab,y.lab Whether to retain x/y axis titles.
#' @param lab.collect Axis-title collection option.
#' @param title.p Optional panel titles.
#' @param save Logical; save using `file` when `TRUE`, otherwise print the plot.
#' @param tag.levels Optional panel labels, such as `"a"`.
#' @param dpi Output resolution in dots per inch. Defaults to 300.
#' @param show.legend One of `"collect"`, `"last"`, `"no"`, `"all"`, or "keep"`.
#' @param legend.position Legend position.
#' @param ... Additional theme elements passed to [ggplot2::theme()].
#' @return The output path invisibly when saving; otherwise the printed plot.
#' @examples
#' library(ggplot2)
#' p <- ggplot(mtcars, aes(wt, mpg)) + geom_point()
#' y_ggsave(p, save = FALSE)
#' @family plot-tools
#' @export
#'

y_ggsave <- function(p, file = "fig/0", nrow = 1, width = 4, height = 3*nrow, spacer = NULL, ylim = NULL, xlim = NULL,
                     x.axis = TRUE, y.axis = TRUE, x.lab = TRUE, y.lab = TRUE, lab.collect = NULL, title.p = NULL, save = TRUE,
                     tag.levels = NULL, dpi = 300, show.legend = c("collect", "last", "no", "all", "keep"), legend.position = "right", ...){

  if(inherits(p, "list") == FALSE){if(!is.null(title.p)){p + ggplot2::labs(title = title.p) -> p}} else{
    if(!is.null(spacer)){for(i in rev(spacer)){append(p, patchwork::plot_spacer(), after = i - 1) -> p}}
    length(p) -> n; if(n/nrow < 1){1 -> nrow}; ceiling(n/nrow) -> ncol;
    if(!x.axis & nrow > 1){purrr::map(p[1:((nrow - 1)*ncol)], ~.x + no_axis_x) -> p[1:((nrow - 1)*ncol)]}
    if(!y.axis){purrr::map(p[-seq(1, n, by = ncol)], ~.x + no_axis_y) -> p[-seq(1, n, by = ncol)]}
    if(!x.lab & nrow > 1){purrr::map(p[1:((nrow - 1)*ncol)], ~.x + no_axis_x_title) -> p[1:((nrow - 1)*ncol)]}
    if(!y.lab){purrr::map(p[-seq(1, n, by = ncol)], ~.x + no_axis_y_title) -> p[-seq(1, n, by = ncol)]}
    if(show.legend[1] == "last"){c(purrr::map(p[1:n - 1], ~.x + no_legend), p[[n]] + ggplot2::theme(legend.position = legend.position)) -> p}
    if(!is.null(title.p)){length(title.p) -> k; if(k < n){purrr::map(1:k, function(i){p[[i]] + ggplot2::labs(title = title.p[i])}) -> p[1:k]
    } else{purrr::map(1:n, function(i){p[[i]] + ggplot2::labs(title = title.p[i])}) -> p}}; patchwork::wrap_plots(p, ncol = ncol, axis_titles = lab.collect) -> p }

   p + patchwork::plot_annotation(tag_levels = tag.levels) -> p;
  if(!is.null(xlim)){p & ggplot2::xlim(xlim) -> p}; if(!is.null(ylim)){p & ggplot2::ylim(ylim) -> p};
  if(show.legend[1] == "collect" & length(legend.position) == 1){p + patchwork::plot_layout(guides = "collect") -> p};
  if(show.legend[1] %in% c("collect", "all")){p & ggplot2::theme(legend.position = legend.position) -> p}; if(show.legend[1] == "no"){p & no_legend -> p}
  p & ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 11), plot.title = ggplot2::element_text(hjust = .5, face = "bold", size = 9), legend.background = ggplot2::element_blank(), ...) -> p
  if(save){dir.create(dirname(paste0(file, ".tiff")), recursive = TRUE, showWarnings = FALSE); ggplot2::ggsave(paste0(file, ".tiff"), p, width = width, height = height, dpi = dpi, units = "in")} else{print(p)}
}
