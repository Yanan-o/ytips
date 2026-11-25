#' @title Describe the distribution of a numeric vector
#' @description Plots Shapiro-Wilk diagnostics, distribution graphics and the
#' Cullen-Frey diagnostic. Combines a Shapiro-Wilk test, boxplot, histogram,
#' density and QQ plot with the Cullen-Frey distribution diagnostic.
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
#' @description Flags observations outside the interval defined by the first
#' and third quartiles plus `multiplier` times the interquartile range.
#'
#' @param x Numeric vector.
#' @param multiplier Numeric multiplier for the interquartile range.
#' @param action `"flag"` returns a logical vector identifying outliers,
#' `"na"` replaces outliers with `NA`, and `"winsorize"` replaces outliers with
#' the nearest boundary.
#' @return A logical vector for `"flag"`, otherwise a numeric vector with the
#' same length as `x`.
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
#' @description Extracts model data, response and fixed-effect predictor names,
#' fixed-effect coefficients, and variance explained. Supports `lm`, `glm`,
#' `lmer`, `glmer`, `glmmTMB` and `brmsfit` models. Use models fitted with
#' `lmerTest::lmer()` to retain mixed-model p-values.
#'
#' @param fit A fitted `lm`, `glm`, `lmer`, `glmer`, `glmmTMB` or `brmsfit`
#' model object.
#' @return A list containing the model data (`df`), response name (`name_y`),
#' fixed-effect predictor names (`name_x`), coefficient table (`Est.txt`) and
#' formatted R-squared values (`rr.txt`).
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
    if("p" %in% names(Est.txt)) {Est.txt %>% mutate(p.sig = cut(p, breaks = c(-Inf, .01, .05, Inf), labels = c("**", "* ", "  ")),
      p.txt = ifelse(p < .01, "~`<`~0.01", paste0("~`=`~", formatC(p, format = "f", digits = 2))) %>% paste0("italic(p)", .)) -> Est.txt
      } else { Est.txt %>% mutate(p.sig = NA, p.txt = NA) -> Est.txt }

    # RR
    performance::r2(fit) %>% unlist() -> RR; if(any(is.na(RR[1]))){MuMIn::r.squaredGLMM(fit)[2] -> RR[1]};
    map(RR, ~{.x = ifelse(.x < 0.01, "~`<`~0.01", paste0("~`=`~", formatC(.x, format = "f", digits = 2)))}) %>% bind_cols() -> RR
    if(inherits(fit, 'lm')|length(RR) == 1){paste0("R^2", RR[[1]]) -> rr.txt} else{
      rbind(paste0("R[m]^2", RR[[2]]), paste0("R[c]^2", RR[[1]])) -> rr.txt}; data.frame(rr.txt = rr.txt) -> rr.txt

   list('df' = df, "name_y" = name_y, "name_x" = name_x, "Est.txt" = Est.txt, "rr.txt" = rr.txt)
 }


#' @title Plot fitted-model coefficients
#' @description Draws fixed-effect estimates and confidence intervals. Residual
#' diagnostics and random-forest variable importance can be added.
#'
#' @param fit A fitted model accepted by [y_fit_summary()].
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param show.random.forest Whether to overlay random-forest variable
#' importance. Requires `rfPermute` and `randomForest`.
#' @param rr.inc Numeric controls for label spacing.
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
    # 检测模型（左黑点应落于红线，右黑线不水平则遗漏非线性项

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
#' @description Returns one table row for a single fitted model. Each
#' coefficient or categorical group forms a column, and R-squared is included
#' as the final column.
#'
#' @param fit A single fitted model accepted by [y_fit_summary()]. A model with
#' several predictors is supported when `show = "est"`. For `show = "mse"`,
#' the model must contain exactly one categorical predictor, such as `y ~ x`.
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
	 data.frame(x = as.character(mse$x), value = paste0(formatC(mse$mean, format = "f", digits = digits), " ± ", formatC(mse$se, format = "f", digits = digits))) -> df_table
  } else {
    as.character(fit_res$Est.txt$p.sig) -> sig; sig[is.na(sig)] <- ""
    sub("^factor\\([^)]*\\)", "", fit_res$Est.txt$x) -> est_x
    paste0(formatC(fit_res$Est.txt$Estimate, format = "f", digits = digits), sig) %>% data.frame(x = est_x, value = .) -> est
    if(!is_categorical){est -> df_table} else{
      data.frame(x = as.character(mse$x), mean = mse$mean, se = mse$se) %>% dplyr::left_join(est, by = "x") %>%
	      mutate(value = ifelse(is.na(value), paste0(formatC(mean, format = "f", digits = digits), " ± ", formatC(se, format = "f", digits = digits)), value)) -> df_table } }

  as.data.frame(as.list(setNames(df_table$value, df_table$x)), check.names = FALSE) %>%
    mutate(R2 = paste0(fit_res$rr.txt[[1]], collapse = " ")) %>% return()
}


#' @title Generate fitted values and confidence intervals
#' @description Returns observed data with fitted values and 95% confidence
#' intervals for supported models. By default, predictions use fixed effects
#' only; `re_formula` controls whether random effects are retained.
#'
#' @param fit A fitted `lm`, `glm`, `lmer`, `glmer`, `glmmTMB` or `brmsfit` model.
#' @param newdata `NULL` uses the model data; a data frame supplies predictor
#' values directly; a character value such as `"grid"` generates a balanced grid
#' across the observed predictor ranges.
#' @param re_formula Random-effect handling. `NA` gives population-level
#' predictions, while `NULL` retains all random effects.
#' @return A data frame containing fitted values and intervals.
#' @details Automatic grids include all fixed-effect predictors. Numeric
#' variables use evenly spaced values across their observed range, while
#' categorical variables retain their observed levels. Up to two fixed-effect
#' predictors use 100 numeric values each; with more predictors, the number of
#' values per predictor is reduced to keep the grid near 10,000 combinations.
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



#' @title Compare subsets of a fitted model
#' @description Fits predictor subsets and ranks them by AICc or predictor
#' p-values. A warning is issued when more than 255 candidate models are
#' requested. output = "best_fit" returns only the selected model.
#'
#' @param fit A fitted `lm`, `lmer` or `glmer` model object.
#' @param step Model-selection rule: AICc, p < 0.1, or p < 0.05.
#' @param output Output mode; see the function description.
#' @return The best model or a list of candidate-model results.
#'
#' @family fitting
#' @export
#'

y_fit_step <- function(fit, step = "aicc", output = "all"){
      # in output No.1, 2 is full_fit & best_fit
      # step could be '0.1' or '0.05' to catch differ methods basing p-values of x;

      all.vars(stats::formula(fit, fixed.only = T))[-1] -> vars.x
      unlist(lapply(1:length(vars.x), \(k) combn(vars.x, k, simplify = FALSE)), recursive = FALSE) %>% c("00", .) -> comb_x

        comb_x %>% imap(~{paste0('-', paste0(.x, collapse = '-'))}) %>%
                   imap(~{paste0('stats::update(stats::formula(fit), .~. ', .x, ")")}) %>%
                   imap(~{eval(parse(text = .x)) %>% as.formula()}) -> formula_all

        suppressMessages(suppressWarnings(formula_all %>% imap(~{stats::update(fit, formula = .x)}))) -> fit_all

        fit_all %>% imap(~{MuMIn::AICc(.x) %>% data.frame(AICc = .)}) %>% bind_rows() -> AICc_all
        fit_all %>% imap(~{stats::formula(.x) %>% deparse1() %>% data.frame(formula = .)}) %>% bind_rows(.id = "fit.id") -> formula_all
        fit_all %>% imap(~{MuMIn::r.squaredGLMM(.x) %>% head(1) %>% data.frame()}) %>% bind_rows() -> RR_all
        fit_all %>% imap(~{anova(fit, .x) %>% data.frame() %>% select(Pr..Chisq.) %>% .[2,1] %>% data.frame(ANOVA.P = .)}) %>% bind_rows() -> ANOVA.P
        fit_all %>% imap(~{summary(.x) %>% .[["coefficients"]] %>% data.frame() %>% rename_with(~ "p", matches("Pr")) %>% .[-1,]}) %>%
          imap(~{cbind(summarise(.x, 'x.p_0.05' = sum(p > 0.05)), summarise(.x, 'x.p_0.1' = sum(p > 0.1)))}) %>% bind_rows() -> Number.x.p

      # step_res: formula, R2m, R2c, AICc, ANOVA.Pr(>Chisq), number of x that p > 0.05 or 0.1
      cbind(formula_all, RR_all, AICc_all, ANOVA.P, Number.x.p) %>% as_tibble() %>% mutate(across(R2m:ANOVA.P, ~round(., 2))) -> step_res

      if(grepl("aicc", step)){step_res %>% arrange(AICc) -> step_res}                      # best_fit is the fit with min AICc
      if(grepl(".1", step)){step_res %>% arrange(x.p_0.1, desc(R2m), AICc) -> step_res}    # best_fit is all p of x < 0.1, and take the max R2m
      if(grepl(".05", step)){step_res %>% arrange(x.p_0.05, desc(R2m), AICc) -> step_res}  # best_fit is all p of x < 0.05, and take the max R2m

      if(step_res[1, 'R2m'] == 0){2 -> id.best_fit; warning('\n0_o| best_fit is y ~ 1\n-----------------------\n')} else{1 -> id.best_fit}
      step_res[id.best_fit, 'fit.id'] %>% as.numeric() %>% fit_all[[.]] -> best_fit

      # all_res: row = 1 is full_fit & row = 2 is best_fit & include all model
      rbind(step_res %>% dplyr::filter(fit.id == 1), step_res) %>% as_tibble() -> all_res
      if (output == "best_fit") {best_fit %>% return()} else{list("best_fit" = best_fit,  "all_res" = all_res) %>% return()}
    }



#' @title Update a model formula while preserving random effects
#' @description Updates a model formula with compact text while retaining the
#' model's fixed-effect terms and random effects. Fixed-effect terms are ordered
#' by their first appearance in `x.order`, and random effects are placed last.
#'
#' @param fit A fitted `lm`, `glm`, `lmer`, `glmer` or `glmmTMB` model.
#' @param str Formula update text, such as `"+ x1 - x2"`.
#' @param x.order Optional character vector giving term order. Defaults to the
#' column order in the model data.
#' @return An updated fitted model with one model fit.
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



#' @title Plot model-averaged coefficients
#' @description Averages models within three AICc units and plots relative
#' importance and coefficients. Models within three AICc units of the best
#' model are averaged. The first plot shows relative importance; the second
#' shows coefficients.
#'
#' @param fit A fitted model accepted by `MuMIn::dredge()`.
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @return A named list containing two ggplot objects.
#'
#' @family fitting
#' @export
#'

y_fit_model_avg <- function(fit, fig.name = NULL){
      # using MuMIn::model.avg to determine Model-averaged coefficients of delta-AICc < 3
      # sw: sum of Akaike weights, indicating its relative importance across models with ΔAICc < 3.

      all.vars(stats::formula(fit, fixed.only = T))[-1] -> vars.x
      MuMIn::dredge(fit, rank = 'AICc') %>% MuMIn::model.avg(., subset = delta < 3) %>% summary() -> m.avg_res

      m.avg_res$coefmat.subset %>% data.frame() %>% rownames_to_column('x') %>% slice(-1) -> df.con.avg
      m.avg_res$sw %>% data.frame(sw = .) %>% rownames_to_column('x') %>% merge(df.con.avg, by = 'x') %>%
        mutate(x = factor(x, levels = vars.x)) %>% arrange(x) %>% mutate(color = ifelse(Pr...z.. < .05, "red", "grey58"))-> df.p

      if(!is.null(fig.name)){fig.name[match(df.p$x, fig.name[, 1]), 2] -> df.p$x}; factor(df.p$x, levels = rev(df.p$x)) -> df.p$x

      ggplot(df.p, aes(x = sw, y = x)) + geom_col(width = .8, alpha = .6, fill = df.p$color) + labs(x = 'Relative importance (Sum of Akaike weights)', y = '') -> p.sw

      ggplot(df.p, aes(x = Estimate, y = x)) + geom_vline(xintercept = 0, linetype = "dotted", size = .5) +
        geom_pointrange(aes(xmin = Estimate - Adjusted.SE, xmax = Estimate + Adjusted.SE), size = .5, color = df.p$color, pch = 15) +
        labs(x = 'Model-averaged coefficients (conditional average +/- SE)', y = '') -> p.avg

      list(p.sw, p.avg) %>% map(~.x + theme_y() + theme(axis.text.y = element_text(size = 9))) %>% set_names(c('sw', 'avg')) %>% return()
    }



#' @title Plot random-forest variable importance
#' @description Runs a permuted random forest and plots variable importance.
#' Permutation importance controls the displayed significance marks.
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
      mutate(p.sig = cut(X.IncMSE.pval, breaks = c(-Inf, .01, .05, Inf), labels = c("**", "*", "")) ) %>%
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
