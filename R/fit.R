#' @title Describe the distribution of a numeric vector
#' @description Draws distribution diagnostics.
#'
#' @param x Numeric vector of observed values.
#' @param discrete Logical; whether `x` contains discrete counts.
#' @return The distribution-diagnostic result.
#'
#' @family fitting
#' @export
#'

y_descdist <- function(x, discrete = FALSE){
  yy_check_pkg("fitdistrplus")

  shapiro.test(x) -> A ; par(mfrow = c(2, 2)) # if p > 0.1, normality, p > .05 is ok
  paste0(A$method, ":\nw = ", round(A$statistic, 2), ' p = ', round(A$p.value, 2)) -> norm.txt
  boxplot(x, horizontal = TRUE); abline(v = mean(x), col = "red"); title(norm.txt); hist(x);
  plot(density(x), col = "red"); qqnorm(x, col = "blue"); qqline(x); par(mfrow = c(1, 1));
  fitdistrplus::descdist(x, discrete = discrete, boot = 999, method = "sample", graph = TRUE)
}



#' @title Detect and optionally process outliers
#' @description Identifies values outside the interquartile range.
#'
#' @param x Numeric vector.
#' @param multiplier Interquartile-range multiplier.
#' @param action `"flag"` returns a logical vector identifying outliers,
#' `"na"` replaces outliers with `NA`, and `"winsorize"` replaces outliers with
#' the nearest boundary.
#' @return A logical or numeric vector, depending on `action`.
#'
#' @family fitting
#' @export
#'

y_is_outlier <- function(x, multiplier = 3, action = c("flag", "na", "winsorize")){

  if (sum(!is.na(x)) < 3) return(if(action[1] == "flag") rep(FALSE, length(x)) else x)

  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE); iqr <- q[2] - q[1]
  if (iqr == 0) return(if(action[1] == "flag") rep(FALSE, length(x)) else x)

  lower <- q[1] - multiplier * iqr; upper <- q[2] + multiplier * iqr
  outlier <- !is.na(x) & (x < lower | x > upper)
  if (action[1] == "flag") return(outlier)
  if (action[1] == "na") {x[outlier] <- NA; return(x)}
  if (action[1] == 'winsorize') {x[outlier & x < lower] <- lower; x[outlier & x > upper] <- upper; return(x)}
}



#' @title Summarise a fitted model
#' @description Returns model data, fixed-effect estimates and R-squared values.
#'
#' @param fit A fitted `lm`, `glm`, `lmer`, `glmer`, `glmmTMB` or `brmsfit`
#' model object.
#' @return A list with `df`, `name_y`, `name_x`, `Est.txt`, numeric R-squared
#' values in `rr`, and formatted labels in `rr.txt`.
#'
#' @examples
#' fit <- lm(mpg ~ wt + hp, data = mtcars)
#' y_fit_summary(fit)$Est.txt
#'
#' @family fitting
#' @export
#'

y_fit_summary <- function(fit){

    insight::get_data(fit) -> df
    insight::find_response(fit)[1] -> name_y
    insight::find_predictors(fit, effects = "fixed", flatten = TRUE) -> name_x

    if(inherits(fit, 'brmsfit')){
        brms::fixef(fit) %>% data.frame() %>% rownames_to_column("x") %>%
          rename_with(~ c("x", "Estimate", "std.Error", "Q2.5", "Q97.5")) %>%
          mutate(p = ifelse(Q2.5 <= 0 & Q97.5 >= 0, 1, 0)) -> Est.txt } else {

    if(inherits(fit, 'glmmTMB')){
      summary(fit)$coefficients$cond -> coefficients} else {
      if(inherits(fit, 'lmerMod')){ summary(lmerTest::as_lmerModLmerTest(fit))$coefficients -> coefficients
      } else {summary(fit)$coefficients -> coefficients} }
      coefficients %>% data.frame() %>% rownames_to_column("x") -> Est.txt
      names(Est.txt)[grepl("rror", names(Est.txt))] <- "std.Error"; names(Est.txt)[grepl("^Pr", names(Est.txt))] <- "p";
      Est.txt %>% mutate(Q2.5 = Estimate - 1.96*std.Error, Q97.5 = Estimate + 1.96*std.Error) -> Est.txt }

    # p_value
    df %>% select(where(~ !is.numeric(.))) %>% names() -> x.name.sub
    Est.txt %>% mutate(x = gsub(paste0(x.name.sub, collapse = "|"), "", x)) -> Est.txt
    Est.txt <- Est.txt[Est.txt$x != "(Intercept)", , drop = FALSE]; rownames(Est.txt) <- NULL
    if("p" %in% names(Est.txt)) {Est.txt %>% mutate(p.sig = cut(p, breaks = c(-Inf, .01, .05, Inf), labels = c("**", "*", " ")),
      p.txt = ifelse(p < .01, "~`<`~0.01", paste0("~`=`~", formatC(p, format = "f", digits = 2))) %>% paste0("italic(p)", .)) -> Est.txt
      } else { Est.txt %>% mutate(p.sig = NA, p.txt = NA) -> Est.txt }

    # RR
    performance::r2(fit) %>% unlist() -> RR; if(any(is.na(RR[1]))){MuMIn::r.squaredGLMM(fit)[2] -> RR[1]}; RR -> R2
    map(RR, ~{.x = ifelse(.x < 0.01, "~`<`~0.01", paste0("~`=`~", formatC(.x, format = "f", digits = 2)))}) %>% bind_cols() -> RR
    if(inherits(fit, 'lm')|length(RR) == 1){paste0("R^2", RR[[1]]) -> rr.txt; data.frame('R2' = R2[1]) -> R2} else{
      rbind(paste0("R[m]^2", RR[[2]]), paste0("R[c]^2", RR[[1]])) -> rr.txt; data.frame('Rm' = R2[2], 'Rc' = R2[1]) -> R2};

   list('df' = df, "name_y" = name_y, "name_x" = name_x, "Est.txt" = Est.txt, "rr" = R2, "rr.txt" = data.frame(rr.txt = rr.txt))
 }



#' @title Plot fitted-model coefficients
#' @description Plots fixed-effect estimates and confidence intervals.
#'
#' @param fit A fitted model accepted by [y_fit_summary()].
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param show.random.forest Whether to add random-forest variable importance.
#' @param rr.inc R-squared label spacing.
#' @param diagnostics Whether to draw model residual diagnostics.
#' @return A ggplot object.
#'
#' @family fitting
#' @export
#'

y_fit_estplot <- function(fit, fig.name = NULL, show.random.forest = FALSE, rr.inc = c(.08, 0), diagnostics = TRUE){

    if(diagnostics){yy_check_pkg("DHARMa")}
    if(show.random.forest){yy_check_random_forest_dependencies()}
    if(isTRUE(diagnostics)){DHARMa::simulateResiduals(fittedModel = fit, n = 300) %>% plot()}
    # Model-fit plot; black points should follow the red line.

    y_fit_summary(fit) -> fit.res
    fit.res$df -> df; fit.res$Est.txt -> Est.txt; fit.res$rr.txt -> rr.txt

    Est.txt %>% mutate(p.color = ifelse(p < .05, "red", "grey58")) -> df.p
    1.1*nrow(Est.txt) - rr.inc[1]*c(0:(nrow(rr.txt)- 1)+rr.inc[2])*(nrow(Est.txt)) -> RR.position

    max(df.p$Q97.5) + .1*(max(df.p$Q97.5) - min(df.p$Q2.5)) -> max.x; ifelse(max.x < .1, .1, max.x) -> max.x
    min(df.p$Q2.5) - .35*(max(df.p$Q97.5) - min(df.p$Q2.5)) -> min.x; ifelse(min.x > .1, -.1, min.x) -> min.x
    if(!is.null(fig.name)){fig.name[match(df.p$x, fig.name[,1]), 2] -> df.p$x}; factor(df.p$x, levels = rev(df.p$x)) -> df.p$x

    ggplot(df.p, aes(y = x, x = Estimate)) +
      geom_vline(xintercept = 0, linetype = "dotted", linewidth = .5) +
      geom_pointrange(aes(xmin = Q2.5, xmax = Q97.5), size = .5, color = df.p$p.color, pch = 15) +
      annotate("text", x = min.x, y = RR.position, hjust = "left", parse = TRUE, size = 2.8, label = rr.txt$rr.txt) +
      theme_y() + theme(axis.text.y = element_text(size = 9)) + xlim(min.x, max.x) + labs(x = "Effect size (Estimate)", y = "") +
      scale_y_discrete(expand = expansion(mult = c(.2, .2))) -> fit.plot

    if(show.random.forest){
      y_random_forest_plot(y = fit.res$name_y, x = fit.res$name_x, df = df, sort = FALSE, output = 'all') -> RF.res
      RF.res$RF.result %>% select(x = name_x, X.IncMSE) %>% merge(df.p, ., by = "x", sort = F, all.x = T) %>%
        mutate(X.IncMSE2 = X.IncMSE/max(X.IncMSE)*max.x) -> df.p.RF
      fit.plot + geom_rect(data = df.p.RF, aes(ymin = as.numeric(x) -.4, ymax = as.numeric(x) +.4, xmin = min.x, xmax = X.IncMSE2), fill = "grey", alpha = .2) +
        annotate('text', x = min.x, y = 0, hjust = "left", vjust = -.8, parse = F, size = 2.5, label = paste0('In.MSE (%): 0 - ', formatC(max(df.p.RF$X.IncMSE), format = "f", digits = 1))) -> fit.plot
    }
    return(fit.plot)
}



#' @title Build a wide fitted-model table
#' @description Builds a one-row coefficient table for a fitted model.
#'
#' @param fit A fitted model accepted by [y_fit_summary()].
#' @param digits Number of decimal places.
#' @param show `"est"` displays coefficient estimates with significance marks;
#' `"mse"` displays observed group means and standard errors.
#' @return A one-row wide data frame.
#'
#' @family fitting
#' @export
#'

y_fit_msetable <- function(fit, digits = 2, show = c("est", "mse")){

  fit %>% y_fit_summary() -> fit_res; fit_res$name_x[1] -> name_x
  names(attr(stats::model.matrix(fit), "contrasts")) -> contrast_terms
  is_categorical <- length(fit_res$name_x) == 1 && (is.factor(fit_res$df[[name_x]]) || is.character(fit_res$df[[name_x]]) ||
    any(vapply(contrast_terms, function(term){name_x %in% all.vars(stats::as.formula(paste0("~", term)))}, logical(1))))

  if(is_categorical){
    data.frame(y = fit_res$df[[fit_res$name_y]], x = fit_res$df[[name_x]]) %>%
      group_by(x) %>% reframe(mean = mean(y, na.rm = TRUE), se = sd(y, na.rm = TRUE)/sqrt(sum(!is.na(y)))) -> mse }

  if(show[1] == "mse"){
	 data.frame(x = as.character(mse$x), value = paste0(formatC(mse$mean, format = "f", digits = digits), " \u00b1 ", formatC(mse$se, format = "f", digits = digits))) -> df_table
  } else {
    as.character(fit_res$Est.txt$p.sig) -> sig; sig[is.na(sig)] <- ""
    sub("^factor\\([^)]*\\)", "", fit_res$Est.txt$x) -> est_x
    paste0(formatC(fit_res$Est.txt$Estimate, format = "f", digits = digits), sig) %>% data.frame(x = est_x, value = .) -> est
    if(!is_categorical){est -> df_table} else{
      data.frame(x = as.character(mse$x), mean = mse$mean, se = mse$se) %>% dplyr::left_join(est, by = "x") %>%
	      mutate(value = ifelse(is.na(value), paste0(formatC(mean, format = "f", digits = digits), " \u00b1 ", formatC(se, format = "f", digits = digits)), value)) -> df_table } }

  paste0(names(fit_res$rr), " = ", formatC(unlist(fit_res$rr), format = "f", digits = digits), collapse = "; ") -> R2
  as.data.frame(as.list(setNames(df_table$value, df_table$x)), check.names = FALSE) %>%
    mutate(R2 = R2) %>% return()
}



#' @title Generate fitted values and confidence intervals
#' @description Adds fitted values and 95% confidence intervals to model data.
#'
#' @param fit A fitted `lm`, `glm`, `lmer`, `glmer`, `glmmTMB` or `brmsfit` model.
#' @param newdata Model data, predictor values, or `"grid"` for a generated grid.
#' @param re_formula `NA` excludes random effects; `NULL` includes them.
#' @return A data frame containing fitted values and intervals.
#'
#' @family fitting
#' @export
#'

y_fit_predict <- function(fit, newdata = NULL, re_formula = NA){
    if(is.null(newdata)){insight::get_data(fit) -> df
    } else if(is.character(newdata)){
      insight::find_predictors(fit, effects = 'fixed', flatten = TRUE) -> vars.x; length(vars.x) -> n.x
      if(n.x <= 2){100 -> grid.length} else{max(5, floor(10000^(1/n.x))) -> grid.length}
      insight::get_datagrid(fit, length = grid.length, range = 'range', protect_integers = FALSE,
                            include_response = FALSE, include_random = FALSE) -> df} else {newdata -> df}

    insight::get_predicted(fit, data = df, predict = 'response', ci = .95,
      include_random = is.null(re_formula)) %>% as.data.frame() %>% select(1:4) %>%
      rename_with(~c("ypre", "ypre.se", "ypre.Q2.5", "ypre.Q97.5")) %>% select(!ypre.se) -> ypre
    df %>% cbind(ypre) %>% as_tibble() %>% return()
  }



#' @title Update a model formula while preserving random effects
#' @description Updates a model formula and keeps its random effects.
#'
#' @param fit A fitted `lm`, `glm`, `lmer`, `glmer` or `glmmTMB` model.
#' @param str Formula update text, such as `"+ x1 - x2"`.
#' @param x.order Optional term order.
#' @return An updated fitted model.
#'
#' @family fitting
#' @export
#'

y_fit_update <- function(fit, str = '', x.order = NULL){
  if(inherits(fit, "merMod")){getData(fit) -> df} else{eval(fit$call$data, environment(stats::formula(fit))) -> df}
  if(is.null(df)){stats::model.frame(fit) -> df}; if(is.null(x.order)){names(df) -> x.order}

  stats::update(stats::formula(fit), stats::as.formula(paste0('. ~ . ', str))) -> new.formula
  reformulas::nobars(new.formula) -> fixed.formula; stats::terms(fixed.formula) -> fixed.terms
  attr(fixed.terms, 'term.labels') -> terms.x
  vapply(terms.x, \(term){min(c(match(all.vars(stats::as.formula(paste0('~', term))), x.order), Inf), na.rm = TRUE)}, numeric(1)) -> term.rank
  order(grepl(':', terms.x, fixed = TRUE), term.rank) -> terms.id
  terms.x[terms.id] %>% paste0(collapse = ' + ') -> fixed.rhs
  if(length(terms.x) == 0){'1' -> fixed.rhs}; if(attr(fixed.terms, 'intercept') == 0){paste0('0 + ', fixed.rhs) -> fixed.rhs}

  reformulas::findbars(new.formula) -> random
  if(length(random) == 0){'' -> random} else{paste0(' + (', random, ')', collapse = '') -> random}
  paste0(paste(deparse(fixed.formula[[2]]), collapse = ''), ' ~ ', fixed.rhs, random) %>%
    stats::as.formula(env = environment(stats::formula(fit))) -> updated.formula
  getCall(fit) -> fit.call; updated.formula -> fit.call$formula
  eval(fit.call, environment(stats::formula(fit))) %>% return()
}



#' @title Rank fixed-effect model subsets
#' @description Ranks fixed-effect model subsets with [MuMIn::dredge()].
#'
#' @param fit A fitted model accepted by [MuMIn::dredge()], preferably fitted
#' with `na.action = na.fail`.
#' @param rank Ranking criterion accepted by [MuMIn::dredge()].
#' @param top Number of ranked subsets, excluding the full model.
#' @return A list with `best_fit` and a model table in `fit_rank`.
#'
#' @examples
#' fit <- lm(mpg ~ wt * hp + disp, data = mtcars, na.action = na.fail)
#' result <- y_fit_models_rank(fit, top = 5)
#' result$fit_rank
#'
#' @family fitting
#' @export
#'

y_fit_models_rank <- function(fit, rank = "AICc", top = 10){
  # AIC comparisons require ML for mixed models; the returned best model is restored to REML.
  is.reml <- inherits(fit, "merMod") && lme4::isREML(fit)
  if(is.reml){stats::update(fit, REML = FALSE) -> fit.select} else{fit -> fit.select}
  suppressMessages(MuMIn::dredge(fit.select, rank = rank)) -> selection

  min(as.integer(top), nrow(selection)) -> n.top
  MuMIn::get.models(selection, subset = seq_len(n.top)) -> top.models
  top.models[[1]] -> best_fit; if(is.reml){stats::update(best_fit, REML = TRUE) -> best_fit}

  c(list(fit.select), top.models) -> all.models
  all.models %>% map(y_fit_summary) -> fit.summaries
  all.models %>% map(~paste(deparse(stats::formula(.x)), collapse = "")) %>% unlist() %>% unname() -> formulas
  fit.summaries %>% map("rr") %>% map(~{
    if(all(c("Rm", "Rc") %in% names(.x))){data.frame(Rm = .x[["Rm"]], Rc = .x[["Rc"]])
    } else{data.frame(R2 = .x[[1]])}
  }) %>% bind_rows() -> RR

  tibble::tibble(rank = c(0L, seq_len(n.top)), formula = formulas, RR,
    AICc = vapply(all.models, MuMIn::AICc, numeric(1)),
    AIC = vapply(all.models, stats::AIC, numeric(1)),
    BIC = vapply(all.models, stats::BIC, numeric(1))) -> fit_rank

  list("best_fit" = best_fit, "fit_rank" = fit_rank) %>% return()
}



#' @title Plot model-averaged coefficients
#' @description Plots coefficients averaged across AICc models within
#' `delta < 3`.
#'
#' @param fit A fitted model accepted by [MuMIn::dredge()].
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @return A list with two ggplot objects: `sw` and `avg`.
#'
#' @family fitting
#' @export
#'

y_fit_model_avg <- function(fit, fig.name = NULL){
  MuMIn::dredge(fit, rank = "AICc") %>%
    MuMIn::model.avg(., subset = delta < 3) %>% summary() -> m.avg_res

  m.avg_res$coefmat.subset %>% data.frame() %>% rownames_to_column("x") %>% slice(-1) -> df.con.avg
  m.avg_res$sw %>% data.frame(sw = .) %>% rownames_to_column("x") %>% merge(df.con.avg, by = "x") %>%
    mutate(x = factor(x, levels = rev(names(m.avg_res$sw)))) %>% arrange(x) %>%
    mutate(color = ifelse(Pr...z.. < .05, "red", "grey58")) -> df.p

  if(!is.null(fig.name)){fig.name[match(df.p$x, fig.name[, 1]), 2] -> df.p$x}
  factor(df.p$x, levels = rev(df.p$x)) -> df.p$x

  ggplot(df.p, aes(x = sw, y = x)) + geom_col(width = .8, alpha = .6, fill = df.p$color) +
    labs(x = "Relative importance (Sum of Akaike weights)", y = "") -> p.sw
  ggplot(df.p, aes(x = Estimate, y = x)) +
    geom_vline(xintercept = 0, linetype = "dotted", linewidth = .5) +
    geom_pointrange(aes(xmin = Estimate - Adjusted.SE, xmax = Estimate + Adjusted.SE),
      linewidth = .5, color = df.p$color, shape = 15) +
    labs(x = "Model-averaged coefficients (conditional average \u00b1 SE)", y = "") -> p.avg

  list(sw = p.sw + theme_y(), avg = p.avg + theme_y()) %>% return()
}



#' @title Plot random-forest variable importance
#' @description Runs a permuted random forest and plots variable importance.
#'
#' @param df Data frame containing the response and all predictor columns.
#' @param x Character vector of predictor column names.
#' @param y Name of one response column.
#' @param x.group Optional data frame with predictor names in the first column
#' and display groups in the second.
#' @param sort Whether to sort variables by importance.
#' @param n.top Maximum number of variables to display.
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param cols Optional color values.
#' @param coord.flip Whether to flip plot coordinates.
#' @param output Output mode; see the function description.
#' @return A ggplot object or a list containing importance results.
#'
#' @family fitting
#' @export
#'

y_random_forest_plot <- function(df, x = "", y = "", x.group = NULL, sort = FALSE, n.top = 9,
                      fig.name = NULL, cols = NULL, coord.flip = TRUE, output = "p"){

    yy_check_random_forest_dependencies()

    # x = c('x1', 'x2', ...); x.group is a df same to fig.name
    if(is.null(cols)){color_y(n = 7, colbar = 1) -> cols}
    df %>% select(all_of(y), all_of(x)) -> df_RF; paste(y, '~', paste(x, collapse = '+')) %>% as.formula() -> RF_formula
    rfPermute::rfPermute(RF_formula, data = df_RF, ntree = 999, num.rep = 200, importance = TRUE) -> RF_fit
    randomForest::randomForest(RF_formula, data = df_RF, ntree = 999, importance = TRUE) %>% capture.output() %>% data.frame() %>% .[9, 1] -> RR

    RF_fit %>% rfPermute::importance() %>% round(3) %>% data.frame() %>% rownames_to_column("name_x") %>%
      mutate(p.sig = cut(X.IncMSE.pval, breaks = c(-Inf, .01, .05, Inf), labels = c("**", "*", " ")) ) %>%
      mutate(name_x = factor(name_x, levels = x)) %>% arrange(name_x) -> RF.result; RF.result -> df_p

    if(sort & is.numeric(n.top)){df_p %>% arrange(desc(X.IncMSE)) %>% slice(1:n.top) -> df_p}
    if(sort & coord.flip){df_p %>% arrange(X.IncMSE) -> df_p}; if(!sort & coord.flip){df_p %>% arrange(desc(name_x)) -> df_p}

    if(is.null(fig.name)){df_p$name_x -> df_p$fig_x} else{fig.name[match(df_p$name_x, fig.name[,1]), 2]-> df_p$fig_x}
    if(is.null(x.group)){mutate(df_p, group_x = ifelse(X.IncMSE.pval < .05, 'red', 'grey68')) -> df_p} else{x.group[match(df_p$name_x, x.group[,1]), 2] -> df_p$group_x}

    df_p %>% mutate(fig_x = factor(fig_x, levels = fig_x)) %>% ggplot(aes(x = fig_x, y = X.IncMSE)) -> p.0
    if(is.null(x.group)){p.0 + geom_col(fill = df_p$group_x, width = .8, alpha = .6) -> p} else{
      p.0 + geom_col(aes(fill = group_x), width = .8, alpha = .6) + geom_text(aes(y = X.IncMSE, label = p.sig), nudge_y = .5, size = 5) +
      scale_fill_manual('', values = cols) -> p};  p + labs(x = '', y = 'Increase in MSE (%)') + theme_y() -> p;

  if(coord.flip){p + coord_flip() -> p} else{p + theme(axis.text.x = element_text(angle = 45, hjust = 1)) -> p}
  if(grepl('p', output)){return(p)} else{return(list("RF.result" = RF.result, "RR" = RR, "RF.plot" = p))}
  }
