#' @title Plot a numeric response by groups
#' @description Plots a numeric response by groups.
#'
#' @param df Data frame whose first two columns are the response and grouping
#' variable; an optional third column defines colors.
#' @param plot Plot type or an existing ggplot object.
#' @param cols Optional color values.
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param error.bar Error bars to use: standard deviation or standard error.
#' @param hline Optional horizontal reference line or "ck" .
#' @param vline Optional vertical reference line.
#' @param point.size Point size.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_classplot <- function(df, plot = c('box', 'bar', 'vio', 'viobox', 'line', 'linebox'), cols = NULL,
                        fig.name = NULL, error.bar = c('sd', 'se'), hline = FALSE, vline = FALSE, point.size = 1.5){
    # df must be y, x, group
    if(plot[1] == 'vio'){yy_check_pkg("ggbeeswarm")}; if(plot[1] == 'viobox'){yy_check_pkg("ggdist")}

    names(df)[1] -> name_y; if(!is.null(fig.name)){fig.name[match(name_y, fig.name[,1]), 2] -> name_y}
    if(!(c('y') %in% names(df))){rename(df, 'y' = 1) -> df}; if(!(c('x') %in% names(df))){rename(df, 'x' = 2) -> df};
    if(!('group' %in% names(df))){"NA" -> df$group}; df %>% mutate(group = factor(group, levels = unique(group))) -> df;
    if(is.null(cols)){color_y(n = length(unique(df$group))) -> cols}; df %>% mutate(x = factor(x, levels = unique(x))) -> df;
    df %>% group_by(x, group) %>% reframe(median = median(y), sd = stats::sd(y), se = sd/sqrt(length(y)), y = mean(y)) -> df_m
    if(error.bar[1] == 'se'){df_m$se -> df_m$sd; df_m$se -> df_m$sd} ; ggplot(data = df, aes(x = x, y = y)) -> p

    if(plot[1] == 'box'){ p + geom_boxplot(aes(fill = group), color = 'black', alpha = .2, width = .4, outlier.size = 1, outlier.shape = 19) -> p}
    if(plot[1] == 'boxx'){ p + geom_boxplot(aes(fill = group, color = group), alpha = .2, width = .4, outlier.size = 1, outlier.shape = 19) -> p}
    if(plot[1] == 'bar'){ p + geom_col(data = df_m, aes(fill = group), width = .7, position = position_dodge(width = .8)) +
      geom_errorbar(data = df_m, aes(ymin = y-sd, ymax = y+sd, color = group), width = .3, size = 1, position = position_dodge(width = .8)) -> p}
    if(plot[1] == 'vio'){ p + geom_violin(aes(fill = group), color = NA, alpha = .2) +
       ggbeeswarm::geom_quasirandom(aes(color = group), width = .32, pch = 19, size = point.size, alpha = .5) +
       geom_point(data = df_m, aes(color = group), pch = 19, size = 2.8) +
        geom_errorbar(data = df_m, aes(ymin = y-sd, ymax = y+sd, color = group), width = 0, size = 1) -> p}
    if(plot[1] == 'viobox'){ p +
       geom_boxplot(aes(fill = group, color = group), alpha = .2, width = .2, outlier.shape = 16, position = position_nudge(x = -.18, y = 0)) +
       ggdist::stat_halfeye(aes(fill = group), adjust = .7, width = .35, .width = 0, alpha = .2, justification = -.18, point_colour = NA) +
       geom_jitter(aes(color = group), size = point.size/2, width = .03, alpha = .5) -> p}
    if(plot[1] == 'line'){ p +
       geom_line(data = df_m, aes(group = group, color = group), linewidth = .5, position = position_dodge(width = .5)) +
       geom_point(data = df_m, aes(color = group), pch = 19, size = 1.5, position = position_dodge(width = .5)) +
       geom_errorbar(data = df_m, aes(ymin = y-sd, ymax = y+sd, color = group), width = 0, size = .5, position = position_dodge(width = .5)) -> p}
    if(plot[1] == 'linebox'){ p +
       geom_line(data = df_m, aes(y = median, group = group, color = group), linewidth = .5, position = position_dodge(width = .5)) +
       geom_boxplot(aes(fill = group, color = group), alpha = .2, width = .39, outlier.size = 1, outlier.shape = 19, position = position_dodge(width = .5)) -> p}

    if(hline == 'ck'){ifelse(grepl(plot[1], 'box'), df_m$median[1], df_m$y[1]) ->  hline };
    if(!isFALSE(hline)){p + geom_hline(yintercept = hline, lty = 2, size = .3) -> p}; if(!isFALSE(vline)){p + geom_vline(xintercept = vline, lty = 2, size =.3) -> p}
    p + labs(x = '', y = name_y) + scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
      theme_y() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) %>% return()
    }



#' @title Plot fitted-model predictions by group
#' @description Plots fitted-model predictions by group.
#'
#' @param fit A fitted model with exactly one categorical predictor.
#' @param plot Plot type or an existing ggplot object.
#' @param pairwise Whether to perform all pairwise comparisons.
#' @param p.adjust P-value adjustment method.
#' @param group Grouping variable, column name, or grouping table.
#' @param ylim Optional y-axis limits.
#' @param hline Whether to add a horizontal guide below the observed points.
#' @param rr.inc Numeric controls for label spacing.
#' @param ... Additional arguments passed to [y_classplot()], including `cols`,
#' `fig.name`, `error.bar`, `vline`, or `point.size`.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_fit_classplot <- function(fit, plot = 'box', pairwise = TRUE, p.adjust = 'default', group = NULL, ylim = NULL, hline = FALSE, rr.inc = c(.08, 0), ...){

     # p.adjust = 'none', p-value is equal to y_fit_summary()
     # only one character x permits in fit; group is a df(x, g) or 'name' in fit; using to coloring point
     yy_check_pkg("car"); if(isTRUE(pairwise)){yy_check_pkg("rcompanion")}

     fit %>% y_fit_summary() -> fit.res;
     fit.res$rr.txt -> rr.txt; fit.res$name_x[1] -> name_x; fit.res$df %>% arrange(across(all_of(name_x))) -> df.0

     if(!is.null(group)){if(is.character(group)){df.0[[group]] -> df.0$'group'} else{group[match(df.0[[name_x]], group[,1]), 2] -> df.0$'group'}}
     rename(df.0, y = all_of(fit.res$name_y), x = all_of(name_x)) -> df; group_by(df, x) %>% reframe(max.y = max(y), mean.y = mean(y)) -> df_m
     if(is.null(ylim)){max(df$y) + .3*(max(df$y) - min(df$y)) -> ymax; min(df$y) -> ymin} else{ylim[1] -> ymin; ylim[2] -> ymax}
     if(!isFALSE(hline)){if(is.null(ylim)){min(df$y) - .03*(max(df$y) - min(df$y)) -> hline.value; min(ymin, hline.value) -> ymin} else{ymin + .01*(ymax - ymin) -> hline.value}} else{NULL -> hline.value}

     # p.sig & model P
      if(p.adjust == 'default'){ifelse(pairwise, 'tukey', 'dunnett') -> p.adjust}
      eval(parse(text = paste0('emmeans::emmeans(fit, ~', name_x, ')') )) -> emm

      if(pairwise){
        emmeans::contrast(emm,  method = "pairwise", adjust = p.adjust) %>% as.data.frame() %>%
        rcompanion::cldList(p.value ~ contrast, reversed = FALSE, data = .) %>% mutate(p.sig = Letter) -> est.p;
       if(length(unique(est.p$p.sig)) == 1){' ' -> est.p$p.sig} } else{
        emmeans::contrast(emm,  method = "trt.vs.ctrl", adjust = p.adjust) %>% as.data.frame() %>%
        mutate(p.sig = ifelse(p.value < .01, "**", ifelse(p.value < .05, "*", " ")), x = gsub(' -.*', '',contrast)) -> est.p }

      car::Anova(fit) %>% data.frame() %>% .[1,] %>%  select(contains('Pr')) %>% rename_with(~'P') %>%
      mutate(P = ifelse(P < 0.01,  "~`<`~0.01", paste0("~`=`~", formatC(P, format = "f", digits = 2))) %>% paste0("italic(p)", .) ) -> P.txt

    # rr.txt & plot
      if(pairwise){ rbind(paste0(rr.txt[,1], collapse = '~~'), P.txt$P) %>% data.frame(rr.txt = .) -> rr.txt
        rr.txt %>% mutate(xposition = .7, yposition = ymax - rr.inc[1]*c(0:(nrow(.)-1)+rr.inc[2])*(ymax-ymin)) -> rr.txt
        est.p %>% cbind(df_m) %>% mutate(yposition = max.y + .05 * (ymax - ymin)) -> Est.txt
       } else{ merge(est.p, df_m, by = 'x', sort = FALSE) %>%
          mutate(p.sig = p.sig) %>% mutate(xstart = 1, xend = x, yposition = 0) %>% dplyr::filter(p.value < 0.05) -> Est.txt
          if(nrow(Est.txt) > 0){for(j in 1:nrow(Est.txt)){Est.txt[, "yposition"][j] = max(df$y) + j*(ymax - max(df$y))/nrow(Est.txt)}}
        rr.txt %>% mutate(xposition = 0.88*length(unique(df$x)), yposition = ymax-rr.inc[1]*c(1:(nrow(.))+rr.inc[2])*(ymax-ymin)) -> rr.txt }

    c(fit.res$name_y, name_x, if("group" %in% names(df.0)){"group"}) -> plot_vars
    y_classplot(df.0[, plot_vars, drop = FALSE], plot = plot, ...) +
      geom_text(data = rr.txt, aes(x = xposition, y = yposition), label = rr.txt[, "rr.txt"], hjust = "left", parse = TRUE, size = 2.5) -> p

    if(pairwise){p + geom_text(data = Est.txt, aes(x = x, y = yposition, label = p.sig), hjust = .5, size = 3) -> p } else {
      if(nrow(Est.txt) > 0){yy_check_pkg("ggpubr"); p + ggpubr::geom_bracket(data = Est.txt, aes(xmin = xstart, xmax = xend, y.position = yposition,
                            label = p.sig), lty = 1, size = .3, label.size = 3.5, tip.length = .0, vjust = .7) -> p }}
  if(!is.null(hline.value)){p + geom_hline(yintercept = hline.value, lty = 2, linewidth = .3) -> p}
  p + ylim(ymin, ymax) %>% return()
  }



#' @title Plot fitted lines and confidence bands
#' @description Plots one or more fitted lines.
#'
#' @param fit.list A fitted model, or named list of models, each with one
#' numeric predictor.
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param cols Optional color values.
#' @param point.size Point size.
#' @param line.size Line width.
#' @param point.alpha Point transparency.
#' @param line.alpha Line transparency.
#' @param ylim Optional y-axis limits.
#' @param show.p Significance threshold used for display.
#' @param show.p.lty Line type for non-significant fitted lines.
#' @param show.se Whether to draw or add confidence bands.
#' @param show.smooth Whether to add a smoothed fitted line.
#' @param rr.loc Optional label position.
#' @param rr.inc Numeric controls for label spacing.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_fit_linesplot <- function(fit.list, fig.name = NULL, cols = NULL, point.size = 2, line.size = 1, point.alpha = .3, line.alpha = 1, ylim = NULL,
                            show.p = 0.05, show.p.lty = 2, show.se = FALSE, show.smooth = FALSE, rr.loc = NULL, rr.inc = c(.08, 0)){
   # fit for one line is OK; fits are same to y_fit_summary, but only one x permits;
   # fit.list: a list includes multiple fit and each fit can be crafted as a line (name of sub-list is the line_name)

    if(class(fit.list)[1] == 'list'){
      fit.list %>% map(y_fit_summary) -> fit.res; c(fit.res[[1]]$name_y, fit.res[[1]]$name_x[1]) -> name_vars

      purrr::map2(fit.res, fit.list, ~{fit_res <- .x;
        fit_res$df %>% select(y = all_of(fit_res$name_y), x = all_of(fit_res$name_x[1]))}) %>% bind_rows(.id = 'line') -> df_fit
      purrr::map2(fit.list, fit.res, ~{name_x <- .y$name_x[1];
        y_fit_predict(.x, newdata = 'grid') %>% select(all_of(c(name_x, 'ypre', 'ypre.Q2.5', 'ypre.Q97.5'))) %>%
        rename(x = all_of(name_x))}) %>% bind_rows(.id = 'line') -> df_pre

      map(fit.res, "Est.txt") %>% imap(~data.frame(line = .y, P = .x$p, P.txt = .x$p.txt)) %>% bind_rows() -> P
      map(fit.res, "rr.txt") %>% imap(~data.frame(line = .y, txt = paste0(.x[,1], collapse = '~~'))) %>%
        bind_rows() %>% dplyr::left_join(P, by = 'line') %>% mutate(txt = paste0(txt, '~~', P.txt)) -> df_txt0
    } else{
      fit.list %>% y_fit_summary() -> fit.res; c(fit.res$name_y, fit.res$name_x[1]) -> name_vars

      fit.res$df %>% select(y = all_of(name_vars[1]), x = all_of(name_vars[2])) %>% mutate(line = 'line') -> df_fit
      y_fit_predict(fit.list, newdata = 'grid') %>% select(all_of(c(name_vars[2], 'ypre', 'ypre.Q2.5', 'ypre.Q97.5'))) %>%
        rename(x = all_of(name_vars[2])) %>% mutate(line = 'line') -> df_pre
      data.frame(txt = rbind(paste0(fit.res$rr.txt[,1], collapse = '~~'), fit.res$Est.txt$p.txt), P = fit.res$Est.txt$p) %>% mutate(line = "line") -> df_txt0 }

    # plot
    df_fit %>% mutate(line = factor(line, levels = unique(line))) -> df_fit
    df_pre %>% mutate(line = factor(line, levels = unique(line))) -> df_pre
    if(is.null(fig.name)){name_vars -> fig_name} else{fig.name[match(name_vars, fig.name[,1]), 2] -> fig_name}
    if(is.null(ylim)){.9*min(df_fit$y, df_pre$ypre.Q2.5) -> ymin; 1.1*max(df_fit$y, df_pre$ypre.Q97.5) -> ymax} else{ylim[1] -> ymin; ylim[2] -> ymax}
    if(is.null(cols)){color_y(n = length(fit.list)) -> cols}; if(is.null(rr.loc)){c(min(df_pre$x), ymax) -> rr.loc}
    df_txt0 %>% mutate(linetype = ifelse(P<show.p, 1, show.p.lty), yposition = rr.loc[2]-rr.inc[1]*c(0:(nrow(.)-1)+rr.inc[2])*(ymax-ymin), xposition = rr.loc[1]) -> df_txt

    ggplot() + geom_point(data = df_fit, aes(x = x, y = y, color = line), pch = 19, alpha = point.alpha, size = point.size) +
      geom_line(data = df_pre, aes(x = x, y = ypre, color = line, lty = line),  size = line.size, alpha = line.alpha) +
      geom_text(data = df_txt, aes(x = xposition, y = yposition, color = line, label = txt), hjust = "left", parse = TRUE, alpha = line.alpha, size = 2.8, show.legend = FALSE) +
      scale_color_manual("Group", values = cols) + scale_linetype_manual(values = df_txt[,"linetype"], guide = "none") + labs(x = fig_name[[2]], y = fig_name[[1]]) -> p

    if(length(fit.list) == 1){p + geom_text(data = df_txt, aes(x = xposition, y = yposition, label = txt), hjust = "left", parse = TRUE,  size = 2.8, show.legend = FALSE) -> p}
    if(show.smooth){p + geom_smooth(data = df_fit, aes(x = x, y = y, color = line, fill = line), size = .5, alpha = .1, lty = 1, se = show.se) + scale_fill_manual("Group", values = cols) -> p}
    if(show.se){p + geom_ribbon(data = df_pre, aes(x = x, ymin = ypre.Q2.5, ymax = ypre.Q97.5, fill = line), alpha = .1) + scale_fill_manual("Group", values = cols) -> p}
    p +  theme_y() + ylim(ymin, ymax) %>% return()
  }



#' @title Add a fitted line to an existing plot
#' @description Adds a fitted curve to an existing ggplot.
#'
#' @param plot Plot type or an existing ggplot object.
#' @param fit A fitted model with one numeric predictor.
#' @param cols Optional color values.
#' @param size Line width.
#' @param alpha Transparency.
#' @param show.p Significance threshold used for display.
#' @param lty.p Line type for non-significant fitted lines.
#' @param show.se Whether to draw or add confidence bands.
#' @param txt Optional label text.
#' @param txt.cols Label color.
#' @param txt.size Label size.
#' @param txt.position Label position.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_fit_addline <- function(plot, fit, cols = '#FF0000', size = 1, alpha = 1, show.p = 0.05, lty.p = 2,
                          show.se = FALSE, txt = NULL, txt.cols = '#FF0000', txt.size = 8, txt.position = NULL){

    fit %>% y_fit_summary() -> fit.res
    fit %>% y_fit_predict(newdata = 'grid') %>% rename(x = all_of(fit.res$name_x[1])) -> df_pre
    fit.res$Est.txt %>% select(Estimate, p, p.sig) -> txt0

    if(is.null(txt)){txt0$p.sig -> txt} else{ txt -> txt}
    if(is.null(txt.position)){c(max(df_pre$x), y = ifelse(txt0$Estimate > 0, max(df_pre$ypre), min(df_pre$ypre))) -> txt.position}
    if(show.se){plot + geom_ribbon(data = df_pre, aes(x = x, ymin = ypre.Q2.5, ymax = ypre.Q97.5), fill = cols, alpha = .1) -> plot}

    plot + geom_line(data = df_pre, aes(x = x, y = ypre),  color = cols, size = size, alpha = alpha, lty = if_else(txt0$p < show.p, 1, lty.p)) +
        annotate("text", x = txt.position[1], y = txt.position[2], label = txt, size = txt.size, color = txt.cols, alpha = alpha) %>% return()
  }



#' @title Plot a two-predictor fitted surface
#' @description Evaluates a fitted model over two numeric predictors.
#'
#' @param fit A fitted model with one response and exactly two numeric predictors.
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param fill.lim Optional limits for the fill scale.
#' @param cols Optional color values.
#' @return A named list containing the two-predictor surface plot (`xxp`)
#' and the combined partial plots (`partial.p`).
#'
#' @family plotting
#' @export
#'

y_fit_surface_plot <- function(fit, fig.name = NULL, fill.lim = NULL, cols = c("steelblue2", "white", "#FF1493")){
    # fits are same to y_fit_summary, but only two x permits

    fit %>% y_fit_summary() -> fit.res; fit.res$df -> df
    fit.res$name_y -> y; fit.res$name_x[1] -> x1; fit.res$name_x[2] -> x2
    if(is.null(fig.name)){c(y, x1, x2) -> fig_name} else{fig.name[match(c(y, x1, x2), fig.name[,1]), 2] -> fig_name}

    y_fit_predict(fit, newdata = 'grid') %>% rename(x1 = all_of(x1), x2 = all_of(x2)) -> df_pre.xx

    ggplot(df_pre.xx, aes(x = x1, y = x2)) + geom_tile(aes(fill = ypre)) +
      geom_contour(aes(z = ypre), color = "white", bins = 30, linewidth = .1) +
      scale_fill_gradientn(colours = colorRampPalette(cols)(8), limits = fill.lim) +
      annotate("text", x = min(df_pre.xx$x1), y = max(df_pre.xx$x2) - .1*c(0,1)*(max(df_pre.xx$x2) - min(df_pre.xx$x2)),
                label = fit.res$rr.txt[[1]], hjust = "left", parse = TRUE,  size = 2.8) +
      labs(x = fig_name[2], y = fig_name[3], fill = fig_name[1]) + theme_y() +
      theme(legend.position = c(.68, .73), legend.title = element_text(hjust = 1)) -> xxplot # ypre of x1 & x2

    expand.grid(x1 = df[[x1]], x2 = median(df[[x2]])) %>% rename_with(~c(x1, x2)) -> new.x1
    expand.grid(x1 = median(df[[x1]]), x2 = df[[x2]]) %>% rename_with(~c(x1, x2)) -> new.x2

    y_fit_predict(fit) %>% rename(y = all_of(y), x1 = all_of(x1), x2 = all_of(x2)) -> df_pre
    y_fit_predict(fit, newdata = new.x1) %>% rename(x1 = all_of(x1), x2 = all_of(x2)) -> df_pre.x1
    y_fit_predict(fit, newdata = new.x2) %>% rename(x1 = all_of(x1), x2 = all_of(x2)) -> df_pre.x2

    ggplot() + geom_point(data = df_pre, aes(x = x1, y = y), alpha = 0.6) + labs(x = fig_name[2], y = fig_name[1], title = 'x2 = median') +
      geom_line(data = df_pre.x1, aes(x = x1, y = ypre), color = "red", linewidth = 1) + theme_classic() -> p1

    ggplot() + geom_point(data = df_pre, aes(x = x2, y = y), alpha = 0.6) + labs(x = fig_name[3], y = fig_name[1], title = 'x1 = median') +
      geom_line(data = df_pre.x2, aes(x = x2, y = ypre), color = "red", linewidth = 1) + theme_classic() -> p2
    p1 + p2 -> p # partial ypre of x1 or x2 (median)

  list("partial.p" = p, "xxp" = xxplot) %>% return()
  }



#' @title Plot a partial residual relationship from a fitted model
#' @description Plots a partial residual relationship from a fitted model.
#'
#' @param fit A fitted `lm` or `lmer` model with at least one numeric
#' fixed-effect predictor.
#' @param x.fixed A single character string naming the numeric fixed-effect
#' predictor whose partial relationship is plotted.
#' @param fig.name Optional two-column data frame: input labels in column 1 and
#' display labels in column 2. Applies to the response and `x.fixed`.
#' @param ... Additional arguments passed to [y_fit_linesplot()], such as
#' `show.p`, `show.se`, `point.size`, or `line.size`.
#' @return A named list with `df_resid`, the model data and residual columns,
#' and `plot`, a ggplot object.
#'
#' @family plotting
#' @export
#'

y_fit_partial_plot <- function(fit, x.fixed = "", fig.name = NULL, ...){

  fit %>% y_fit_summary() -> fit.res; fit.res$df -> df; fit.res$name_y -> y
  stats::update(fit, formula. = stats::as.formula(paste0('. ~ . - ', x.fixed))) -> fit1
  stats::update(fit1, formula. = stats::as.formula(paste0(x.fixed, ' ~ .'))) -> fit2

  df %>% mutate(residual.y = stats::resid(fit1), residual.x.fixed = stats::resid(fit2)) -> df_residuals

  # plot
  lm(residual.y ~ residual.x.fixed, data = df_residuals) -> fit_residual # 此处不需要再加random
  c(y, x.fixed) -> name.fig; if(!is.null(fig.name)){fig.name[match(name.fig, fig.name[,1]), 2] -> name.fig}

  fit_residual %>% y_fit_linesplot(fig.name = NULL, ...) +
    labs(y = paste0("Residuals of ", name.fig[1]), x = paste0("Residuals of ", name.fig[2])) -> p
  return(list("df_resid" = df_residuals, "plot" = p))
}



#' @title Decompose explanatory variation
#' @description Compares explanatory variation for two response components.
#'
#' @param fit A fitted `lm` or `lmer` model.
#' @param y1,y2 Single character strings naming the two numeric response
#' components. Both columns must be present in the data used to fit `fit`.
#' @param name.yy Character vector of two display labels for `y1` and `y2`.
#' @return A named list with `df_p`, containing the relative explanatory power
#' and labels, and `plot`, a ggplot object.
#'
#' @family plotting
#' @export
#'

y_fit_decomposition_plot <- function(fit, y1 = "", y2 = "", name.yy = c(y1, y2)){

  var_pred <- function(model){
    if(inherits(model, "merMod")){stats::predict(model, re.form = NA) -> p} else{stats::predict(model) -> p}
    var(p) }

  fit %>% y_fit_summary() -> fit.res
  tryCatch(eval(getCall(fit)$data, environment(stats::formula(fit))), error = function(e) NULL) -> df
  if(!is.data.frame(df)){fit.res$df -> df}
  stats::update(fit, formula. = stats::as.formula(paste0(y1, ' ~ .')), data = df) -> fit_y1
  stats::update(fit, formula. = stats::as.formula(paste0(y2, ' ~ .')), data = df) -> fit_y2
  var_total <- var_pred(fit) + if(inherits(fit, "merMod")){sum(as.data.frame(VarCorr(fit))$vcov)} else{stats::sigma(fit)^2}
  var_y <- vapply(list(fit_y1, fit_y2), var_pred, numeric(1))
  data.frame(explanatory.power = var_y/sum(var_y), name = name.yy) -> df_p

  ggplot(df_p, aes(x = .3, y = explanatory.power, fill = rev(name))) +
    geom_bar(stat = "identity", width = .5, alpha = .5) +
    geom_text(aes(x = .7, label = rev(name)), hjust = .1, size = 2.5) +
    scale_fill_manual("", values = rev(c("darkgreen", "darkorange2"))) +
    coord_flip() +  labs(y = "Relative explanatory power (%)") + xlim(c(0, .8)) +
    theme_classic() + theme(legend.position = "none") + theme(plot.margin = margin(-2, 0, 0, 0)) +
    theme(axis.line.y = element_blank(), plot.background = element_blank()) + no_axis_y -> p
  list("df_p" = df_p, "plot" = p) %>% return()
}



#' @title Draw a correlation bubble plot
#' @description Draws a correlation bubble matrix.
#'
#' @param df Data frame with x labels, y labels, correlations, and optional
#' significance labels.
#' @param fill.lim Optional limits for correlation colors.
#' @param point.size Bubble-size multiplier.
#' @param fig.name Optional two-column data frame: input labels in column 1 and
#' display labels in column 2. Applies to both axes.
#' @param legend.title Character string used as the legend title.
#' @param base.size Base text size.
#' @param legend.direction Legend direction, `"v"` or `"h"`.
#' @param cols Colors for low, middle, and high correlations.
#' @return A ggplot object.
#'
#' @examples
#' cor <- data.frame(
#'   x = c("a", "b", "a", "b"),
#'   y = c("a", "a", "b", "b"),
#'   r = c(1, .4, .4, 1),
#'   p = c("", "*", "*", "")
#' )
#' y_corrplot(cor)
#'
#' @family plotting
#' @export
#'

y_corrplot <- function(df, fill.lim = NULL, point.size = .8, fig.name = NULL, legend.title = ' ',
                       base.size = 10, legend.direction = 'v', cols = NULL){

  if(is.null(cols)){colorRampPalette(c("steelblue2", "white", "#FF1493"))(8) -> cols}
  if(ncol(df) == 3){mutate(df, p.sig = ' ') -> df}; df %>% rename_with(~c('name_x', 'name_y', 'value', 'p.sig')) -> df
  if(!is.null(fig.name)){fig.name[match(df$name_x, fig.name[,1]), 2] -> df$name_x; fig.name[match(df$name_y, fig.name[,1]), 2] -> df$name_y;}
  df %>% mutate(name_x = factor(name_x, levels = unique(name_x)), name_y = factor(name_y, levels = rev(unique(name_y))) ) -> df

  if(is.null(fill.lim)){c(min(df$value), max(df$value)) -> fill.lim}; seq(0, max(abs(fill.lim)), length.out = 8) -> limits.cut
  df %>% mutate(value_cut = cut(abs(value), breaks = c(limits.cut, Inf), labels = FALSE, include.lowest = TRUE)) -> df.p

  ggplot(df.p, aes(x = name_x, y = name_y)) +
    geom_tile(color = 1, fill = "white", linewidth = .05) +
    geom_point(aes(fill = value), pch = 22,  colour = 'white', size = point.size * df.p$value_cut) +
    geom_text(aes(label = p.sig), color = 1, size = 5, nudge_y = -.18, hjust = .5) + labs(x = "", y = "") -> p

  if(legend.direction[1] == 'v'){guide_colorbar() -> legend.guide} else{
      guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = .5 , label.position = "bottom",
                     label.theme = element_text(angle = 70, size = 8,  hjust = 1, vjust = 1.02)) -> legend.guide}

    p + scale_fill_gradientn(name = legend.title, colours = cols, limits = fill.lim, guide = legend.guide) +
      theme_minimal(base.size) + theme(panel.grid = element_blank()) + theme(plot.background = element_rect(fill = "white", linetype = 0)) +
      theme(axis.text = element_text(color = 1)) + theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1.02)) %>%
    return()
}
