#' @title Plot a numeric response by groups
#' @description Supports box, bar, violin, line and combined group plots. The
#' first two data columns are treated as response and grouping variables. If
#' present, the third column is used for color groups.
#'
#' @param df Input data frame.
#' @param plot Plot type or an existing ggplot object.
#' @param cols Optional color values.
#' @param fig.name Optional two-column lookup table mapping input labels to
#' display labels.
#' @param error.bar Error bars to use: standard deviation or standard error.
#' @param hline Optional horizontal reference line or "ck" .
#' @param vline Optional vertical reference line.
#' @param point.size Point size.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_classplot <- function(df, plot = c('box', 'bar', 'vio', 'viobox', 'line', 'linebox'), cols = NULL,
                        fig.name = NULL, error.bar = c('sd', 'se'), hline = FALSE, vline = FALSE, point.size = 1.5, ...){
    # df must be y, x, group
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
#' @description Combines fitted predictions with class-level plotting. The
#' fitted model should contain one categorical predictor. Significance letters
#' or brackets are added according to pairwise. Significance letters or
#' brackets are added according to pairwise.
#'
#' @param fit A fitted model object supported by the function.
#' @param plot Plot type or an existing ggplot object.
#' @param pairwise Whether to perform all pairwise comparisons.
#' @param p.adjust P-value adjustment method.
#' @param group Grouping variable, column name, or grouping table.
#' @param ylim Optional y-axis limits.
#' @param hline Whether to add a horizontal guide below the observed points.
#' @param rr.inc Numeric controls for label spacing.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_fit_classplot <- function(fit, plot = 'box', pairwise = TRUE, p.adjust = 'default', group = NULL, ylim = NULL, hline = FALSE, rr.inc = c(.08, 0), ...){

     # p.adjust = 'none', p-value is equal to y_fit_summary()
     # only one character x permits in fit; group is a df(x, g) or 'name' in fit; using to coloring point
     fit %>% y_fit_summary(output = 'all', diagnostics = FALSE) -> fit.res; fit.res$rr.txt -> rr.txt; fit.res$df %>% arrange(across(2)) -> df.0
     if(!is.null(group)){if(is.character(group)){df.0[[group]] -> df.0$'group'} else{group[match(df.0[[2]], group[,1]), 2] -> df.0$'group'}}
     rename(df.0,'y' = 1, 'x' = 2) -> df; group_by(df, x) %>% reframe(max.y = max(y), mean.y = mean(y)) -> df_m
     if(is.null(ylim)){max(df$y) + .3*(max(df$y) - min(df$y)) -> ymax; min(df$y) -> ymin} else{ylim[1] -> ymin; ylim[2] -> ymax}
     if(!isFALSE(hline)){if(is.null(ylim)){min(df$y) - .03*(max(df$y) - min(df$y)) -> hline.value; min(ymin, hline.value) -> ymin} else{ymin + .01*(ymax - ymin) -> hline.value}} else{NULL -> hline.value}

     # p.sig & model P
      if(p.adjust == 'default'){ifelse(pairwise, 'tukey', 'dunnett') -> p.adjust}
      eval(parse(text = paste0('emmeans::emmeans(fit, ~', names(df.0)[2], ')') )) -> emm

      if(pairwise){
        emmeans::contrast(emm,  method = "pairwise", adjust = p.adjust) %>% as.data.frame() %>%
        rcompanion::cldList(p.value ~ contrast, reversed = FALSE, data = .) %>% mutate(p.sig = Letter) -> est.p;
       if(length(unique(est.p$p.sig)) == 1){' ' -> est.p$p.sig} } else{
        emmeans::contrast(emm,  method = "trt.vs.ctrl", adjust = p.adjust) %>% as.data.frame() %>%
        mutate(p.sig = ifelse(p.value < .01, "**", ifelse(p.value < .05, '*', " ")), x = gsub(' -.*', '',contrast)) -> est.p }

      car::Anova(fit) %>% data.frame() %>% .[1,] %>%  select(contains('Pr')) %>% rename_with(~'P') %>%
      mutate(P = ifelse(P < 0.01,  "~`<`~0.01", paste0("~`=`~", round(P, 2))) %>% paste0("italic(p)", .) ) -> P.txt

    # rr.txt & plot
      if(pairwise){ rbind(paste0(rr.txt[,1], collapse = '~~'), P.txt$P) %>% data.frame(rr.txt = .) -> rr.txt
        rr.txt %>% mutate(xposition = .7, yposition = ymax - rr.inc[1]*c(0:(nrow(.)-1)+rr.inc[2])*(ymax-ymin)) -> rr.txt
        est.p %>% cbind(df_m) %>% mutate(yposition = max.y + .05 * (ymax - ymin)) -> Est.txt
       } else{ merge(est.p, df_m, by = 'x', sort = FALSE) %>%
          mutate(p.sig = p.sig) %>% mutate(xstart = 1, xend = x, yposition = 0) %>% dplyr::filter(p.value < 0.05) -> Est.txt
          if(nrow(Est.txt) > 0){for(j in 1:nrow(Est.txt)){Est.txt[, "yposition"][j] = max(df$y) + j*(ymax - max(df$y))/nrow(Est.txt)}}
        rr.txt %>% mutate(xposition = 0.88*length(unique(df$x)), yposition = ymax-rr.inc[1]*c(1:(nrow(.))+rr.inc[2])*(ymax-ymin)) -> rr.txt }

    y_classplot(df.0, plot = plot, ...) +
      geom_text(data = rr.txt, aes(x = xposition, y = yposition), label = rr.txt[, "rr.txt"], hjust = "left", parse = TRUE, size = 2.5) -> p

    if(pairwise){p + geom_text(data = Est.txt, aes(x = x, y = yposition, label = p.sig), hjust = .5, size = 3) -> p } else {
      if(nrow(Est.txt) > 0){p + ggpubr::geom_bracket(data = Est.txt, aes(xmin = xstart, xmax = xend, y.position = yposition,
                                label = p.sig), lty = 1, size = .3, label.size = 3.5, tip.length = .0, vjust = .7) -> p }}
  if(!is.null(hline.value)){p + geom_hline(yintercept = hline.value, lty = 2, linewidth = .3) -> p}
  p + ylim(ymin, ymax) %>% return()
  }



#' @title Plot fitted lines and confidence bands
#' @description Plots one or more fitted lines.
#'
#' @param fit.list A fitted model or a named list of fitted models.
#' @param fig.name Optional two-column lookup table mapping input labels to
#' display labels.
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
   # fit for one line is OK; fits are same to y_fit_summary, but only one x permits; show.p.lty = NA: lines p>.05 not show
   # lines.list: a list includes multiple fit and each fit can be crafted as a line (name of sub-list is the line_name)

    if(class(fit.list)[1] == 'list'){
      fit.list %>% map(~{y_fit_summary(.x, output = 'all', diagnostics = FALSE)}) -> fit.res; fit.res[[1]] %>% .$df %>% names() -> name_vars
      fit.list %>% map(~{.x %>% y_fit_predict() %>% rename('y' = 1, 'x' = 2)}) %>%
        bind_rows(.id = "line") %>% mutate(line = factor(line, levels = unique(line)))-> df_pre
        map(fit.res, "Est.txt") %>% map(~data.frame(P = .x$p, P.txt = .x$p.txt)) %>% bind_rows(.id = 'line') -> P
        map(fit.res, "rr.txt") %>% map(~data.frame(txt = paste0(.x[,1], collapse = '~~') )) %>% bind_rows() %>%
        cbind(., P) %>% mutate(txt = paste0(txt, '~~', P.txt)) -> df_txt0
      } else{
      fit.list %>% y_fit_summary(output = 'all', diagnostics = FALSE) -> fit.res; fit.res$df -> df; names(df) -> name_vars
      fit.list %>% y_fit_predict() %>% rename('y' = 1, 'x' = 2) %>% mutate(line = "line") -> df_pre
      data.frame(txt = rbind(paste0(fit.res$rr.txt[,1], collapse = '~~'), fit.res$Est.txt$p.txt), P = fit.res$Est.txt$p) %>% mutate(line = "line") -> df_txt0 }

    # plot
    if(is.null(fig.name)){name_vars -> fig_name} else{fig.name[match(name_vars, fig.name[,1]), 2] -> fig_name}
    if(is.null(ylim)){.9*min(df_pre$y, df_pre$ypre.Q2.5) -> ymin; 1.1*max(df_pre$y, df_pre$ypre.Q97.5) -> ymax} else{ylim[1] -> ymin; ylim[2] -> ymax}
    if(is.null(cols)){color_y(n = length(fit.list)) -> cols}; if(is.null(rr.loc)){c(min(df_pre$x), ymax) -> rr.loc}
    df_txt0 %>% mutate(linetype = ifelse(P<show.p, 1, show.p.lty), yposition = rr.loc[2]-rr.inc[1]*c(0:(nrow(.)-1)+rr.inc[2])*(ymax-ymin), xposition = rr.loc[1]) -> df_txt

    ggplot(df_pre) + geom_point(aes(x = x, y = y, color = line), pch = 19, alpha = point.alpha, size = point.size) +
      geom_line(aes(x = x, y = ypre, color = line, lty = line),  size = line.size, alpha = line.alpha) +
      geom_text(data = df_txt, aes(x = xposition, y = yposition, color = line, label = txt), hjust = "left", parse = TRUE, alpha = line.alpha, size = 2.8, show.legend = FALSE) +
      scale_color_manual("Group", values = cols) + scale_linetype_manual(values = df_txt[,"linetype"], guide = "none") + labs(x = fig_name[[2]], y = fig_name[[1]]) -> p

    if(length(fit.list) == 1){p + geom_text(data = df_txt, aes(x = xposition, y = yposition, label = txt), hjust = "left", parse = TRUE,  size = 2.8, show.legend = FALSE) -> p}
    if(show.smooth){p + geom_smooth(aes(x = x, y = y, color = line, fill = line), size = .5, alpha = .1, lty = 1, se = show.se) + scale_fill_manual("Group", values = cols) -> p}
    if(show.se){p + geom_ribbon(aes(x = x, ymin = ypre.Q2.5, ymax = ypre.Q97.5, fill = line), alpha = .1) + scale_fill_manual("Group", values = cols) -> p}
    p +  theme_y() + ylim(ymin, ymax) %>% return()
  }



#' @title Add a fitted line to an existing plot
#' @description Adds a fitted curve to an existing ggplot.
#'
#' @param plot Plot type or an existing ggplot object.
#' @param fit A fitted model object supported by the function.
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

    fit %>% y_fit_predict() %>% rename('y' = 1, 'x' = 2) -> df_pre
    fit %>% y_fit_summary(output = 'all', diagnostics = FALSE) %>% .$Est.txt %>% select(Estimate, p, p.sig) -> txt0

    if(is.null(txt)){txt0$p.sig -> txt} else{ txt -> txt}
    if(is.null(txt.position)){c(max(df_pre$x), y = ifelse(txt0$Estimate > 0, max(df_pre$ypre), min(df_pre$ypre))) -> txt.position}
    if(show.se){plot + geom_ribbon(data = df_pre, aes(x = x, ymin = ypre.Q2.5, ymax = ypre.Q97.5), fill = cols, alpha = .1) -> plot}

    plot + geom_line(data = df_pre, aes(x = x, y = ypre),  color = cols, size = size, alpha = alpha, lty = if_else(txt0$p < show.p, 1, lty.p)) +
        annotate("text", x = txt.position[1], y = txt.position[2], label = txt, size = txt.size, color = txt.cols, alpha = alpha) %>% return()
  }



#' @title Plot a two-predictor fitted surface
#' @description Evaluates a fitted model over two numeric predictors.
#'
#' @param fit A fitted model object supported by the function.
#' @param fig.name Optional two-column lookup table mapping input labels to
#' display labels.
#' @param fill.lim Optional limits for the fill scale.
#' @param cols Optional color values.
#' @param output Output mode; see the function description.
#' @return A ggplot object or a list containing partial plots.
#'
#' @family plotting
#' @export
#'

y_fit_surface_plot <- function(fit, fig.name = NULL, fill.lim = NULL, cols = c("steelblue2", "white", "#FF1493"), output = ""){
    # fits are same to y_fit_summary, but only two x permits

    fit %>% y_fit_summary(output = 'all', diagnostics = FALSE) -> fit.res; fit.res$df -> df
    names(df)[1] -> y; names(df)[2] -> x1;  names(df)[3] -> x2
    if(is.null(fig.name)){names(df) -> fig_name} else{fig.name[match(names(df), fig.name[,1]), 2] -> fig_name}

    expand.grid(
      x1 = seq(min(df[[x1]]), max(df[[x1]]), length = 100),
      x2 = seq(min(df[[x2]]), max(df[[x2]]), length = 100)) %>% rename_with(~c(x1, x2))-> new.xx

    y_fit_predict(fit, newdata = new.xx) %>% rename('x1' = 1, 'x2' = 2) -> df_pre.xx

    ggplot(df_pre.xx, aes(x = x1, y = x2)) + geom_tile(aes(fill = ypre)) +
      geom_contour(aes(z = ypre), color = "white", bins = 30, linewidth = .1) +
      scale_fill_gradientn(colours = colorRampPalette(cols)(8), limits = fill.lim) +
      annotate("text", x = min(df_pre.xx$x1), y = max(df_pre.xx$x2) - .1*c(0,1)*(max(df_pre.xx$x2) - min(df_pre.xx$x2)),
                label = fit.res$rr.txt[[1]], hjust = "left", parse = TRUE,  size = 2.8) +
      labs(x = fig_name[2], y = fig_name[3], fill = fig_name[1]) + theme_y() +
      theme(legend.position = c(.68, .73), legend.title = element_text(hjust = 1)) -> xxplot # ypre of x1 & x2

    expand.grid(x1 = df[[x1]], x2 = median(df[[x2]])) %>% rename_with(~c(x1, x2)) -> new.x1
    expand.grid(x1 = median(df[[x1]]), x2 = df[[x2]]) %>% rename_with(~c(x1, x2)) -> new.x2

    y_fit_predict(fit) %>% rename('y' = 1, 'x1' = 2, 'x2' = 3) -> df_pre
    y_fit_predict(fit, newdata = new.x1) %>% rename('x1' = 1, 'x2' = 2) -> df_pre.x1
    y_fit_predict(fit, newdata = new.x2) %>% rename('x1' = 1, 'x2' = 2) -> df_pre.x2

    ggplot() + geom_point(data = df_pre, aes(x = x1, y = y), alpha = 0.6) + labs(x = fig_name[2], y = fig_name[1], title = 'x2 = median') +
      geom_line(data = df_pre.x1, aes(x = x1, y = ypre), color = "red", linewidth = 1) + theme_classic() -> p1

    ggplot() + geom_point(data = df_pre, aes(x = x2, y = y), alpha = 0.6) + labs(x = fig_name[3], y = fig_name[1], title = 'x1 = median') +
      geom_line(data = df_pre.x2, aes(x = x2, y = ypre), color = "red", linewidth = 1) + theme_classic() -> p2
    p1 + p2 -> p # partial ypre of x1 or x2 (median)

    if(output == ""){xxplot %>% return()} else{list("partial.p" = p, "xxp" = xxplot) %>% return()}
  }



#' @title Plot a partial residual relationship
#' @description Plots residualised relationships. The focal and response
#' variables are residualized against the remaining predictors before their
#' relationship is plotted.
#'
#' @param df Input data frame.
#' @param y Response variable name or response values.
#' @param x Numeric or character input used by the function.
#' @param x.fix Predictor to retain as the focal partial relationship.
#' @param fit A fitted model object supported by the function.
#' @param show.p Significance threshold used for display.
#' @return A list containing residual data and a ggplot object.
#'
#' @family plotting
#' @export
#'

y_partial_plot <- function(df, y = "", x = "", x.fix = "", fit = "lmer", show.p = 0.05){

  paste0(y, " ~ ", paste0(x[x != x.fix], collapse = " + ")) -> fit_formula1
  paste0(x.fix, " ~ ", paste0(x[x != x.fix], collapse = " + ")) -> fit_formula2

  if(fit == "lmer"){paste0(fit_formula1, " + (1|block)") %>% as.formula() %>% lmer(., data = df) -> fit1;
    paste0(fit_formula2, " + (1|block)") %>% as.formula() %>% lmer(., data = df) -> fit2} else{
      fit_formula1 %>% as.formula %>% lm(., data = df) -> fit1;
      fit_formula2 %>% as.formula %>% lm(., data = df) -> fit2}

  data.frame(stats::resid(fit1), stats::resid(fit2)) %>% rename_with(~c(paste0("residuals.", y), paste0("residuals.", x.fix))) -> df_residuals
  if(fit == "lmer"){df_residuals %>% mutate(block = df$block) -> df_residuals}

  # plot
  lm(df_residuals[,1] ~ df_residuals[,2]) -> fit
  summary(fit)$r.squared %>% round(2) -> RR ; anova(fit)$"Pr(>F)"[1] %>% round(2) -> p
  data.frame(xmin = min(df_residuals[,2]), ymax = max(df_residuals[,1]),
             txt = paste0("R^2~`=`~", RR, "~~~italic(p)~`=`~", p)) -> lm_txt

  df_residuals %>%
    ggplot(., aes(x = .[[2]], y = .[[1]])) + geom_point(color = "grey20", size = 2) +
    geom_smooth(method = "lm", se = T, color = 1, fill = "grey", size = 0.8, lty =  if_else(p < show.p, 1, 2)) +
    geom_text(data = lm_txt, aes(x = xmin, y = ymax), label = lm_txt$txt, hjust = "left", parse = T,  size = 2.8) +
    labs(y = paste0("Residuals of ", y), x = paste0("Residuals of ", x.fix)) + theme_classic() -> resid.plot

  return(list("df_resid" = df_residuals, "plot" = resid.plot))
}



#' @title Decompose explanatory variation
#' @description Partitions explanatory power between two response components.
#' The data should contain block, x.fixed, y, and the two components named by
#' y1y2.name. The data should contain block, x.fixed, y, and the two
#' components named by y1y2.name.
#'
#' @param df Input data frame.
#' @param y1y2.name Names of the two decomposed response components.
#' @param output Output mode; see the function description.
#' @return A ggplot object or a list containing its data.
#'
#' @family plotting
#' @export
#'

y_decomposition_plot <- function(df, y1y2.name = c("y1", "y2"), output = "plot"){
  # df includes block, x.fixed, y, y1, y2
  # y ~ x.fixed + (1|block); decomposing y into y1, y2

  lmer(y ~ x.fixed + (1|block), data = df) -> fit1
  stats::model.matrix(fit1) %*% fixef(fit1) -> fit_pre1 # same to: stats::predict(fit1, re.form = NA)
  fit_pre1 %>% var() -> var_fixef; var_fixef

  VarCorr(fit1) %>% as.data.frame() -> var_rb; var_rb
  var_rb$vcov[1] + var_rb$vcov[2] -> var_Error; var_Error # (fit_pre1 - df_FD_Height$p.FD_Height.norm) %>% var()
  var_fixef + var_Error -> sum_var_FD; sum_var_FD

  lmer(y1 ~ x.fixed + (1|block), data = df) -> fit2
  lmer(y2 ~ x.fixed + (1|block), data = df) -> fit3
  stats::model.matrix(fit2) %*% fixef(fit2) %>% var() -> var_fit_pre2
  stats::model.matrix(fit3) %*% fixef(fit3) %>% var() -> var_fit_pre3

  expla_inter = var_fit_pre2/sum_var_FD
  expla_intra = var_fit_pre3/sum_var_FD

  rbind(
    expla_inter/(expla_inter + expla_intra), expla_intra/(expla_inter + expla_intra)) %>%
    data.frame(explanatory.power = .) %>% mutate(name = y1y2.name) -> df_p

  ggplot(df_p, aes(x = .3, y = explanatory.power, fill = rev(name))) +
    geom_bar(stat = "identity", width = .5, alpha = .5) +
    geom_text(aes(x = .7, label = rev(name)), hjust = .1, size = 2.5) +
    scale_fill_manual("", values = rev(c("darkgreen", "darkorange2"))) +
    coord_flip() +  labs(y = "Relative explanatory power (%)") + xlim(c(0, .8)) +
    theme_classic() + theme(legend.position = "none") + theme(plot.margin = margin(-2, 0, 0, 0)) +
    theme(axis.line.y = element_blank(), plot.background = element_blank()) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()) -> p

  if(output == "all"){list(df_p, p) %>% return()} else{p %>% return()}
}



#' @title Draw a correlation bubble plot
#' @description Draws a bubble matrix for correlations and significance. Input
#' columns are label-x, label-y, correlation value, and an optional
#' significance label.
#'
#' @param df Input data frame.
#' @param fill.lim Optional correlation color limits.
#' @param point.size Point size.
#' @param legend.title Legend title.
#' @param legend.direction Legend direction.
#' @param cols Optional color values.
#' @return A ggplot object.
#'
#' @family plotting
#' @export
#'

y_corrplot <- function(df, fill.lim = NULL, point.size = .66, legend.title = ' ',
                       legend.direction = 'v', cols = c("#FF1493", "white", "steelblue2")){

  # df includes name_x, name_y, value, p.sig

  if(ncol(df) == 3){mutate(df, p.sig = ' ') -> df}; df %>% rename_with(~c('name_x', 'name_y', 'value', 'p.sig')) %>%
    mutate(name_x = factor(name_x, levels = unique(name_x)), name_y = factor(name_y, levels = rev(unique(name_y))) ) -> df

  if(is.null(fill.lim)){c(min(df$value), max(df$value)) -> fill.lim}; if(fill.lim[1] < 0 && fill.lim[2] > 0){
    seq(0, max(abs(fill.lim)), length.out = 8) -> limits.cut} else{seq(fill.lim[1], fill.lim[2], length.out = 8) -> limits.cut}

  df %>% mutate(value_cut = cut(abs(value), breaks = c(limits.cut, Inf), labels = FALSE)) -> df.p

  ggplot(df.p, aes(x = name_x, y = name_y)) +
    geom_tile(color = 1, fill = "white", size = .05) +
    geom_point(aes(fill = value), pch = 22,  colour = alpha('white', 0), size = point.size * df.p$value_cut) +
    geom_text(aes(label = p.sig), color = 1, size = 5, nudge_y = -.18, nudge_x = 0) + labs(x = "", y = "") -> p

  if(legend.direction != "horizontal"){guide_colorbar() -> legend.guide} else{
      guide_colorbar(direction = "horizontal", title.position = "top", title.hjust = .5 , label.position = "bottom",
                     label.theme = element_text(angle = 70, size = 8,  hjust = 1, vjust = 1.02)) -> legend.guide}

    p + scale_fill_gradientn(name = legend.title, colours = colorRampPalette(cols)(8), limits = fill.lim, guide = legend.guide) +

      theme_minimal() + theme(panel.grid = element_blank()) + theme(plot.background = element_rect(fill = "white", linetype = 0)) +
      theme(axis.text = element_text(color = 1)) + theme(axis.text.x = element_text(angle = 70, hjust = 1, vjust = 1.02)) %>%
    return()
}
