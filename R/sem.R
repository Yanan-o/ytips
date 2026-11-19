#' @title Build a piecewiseSEM formula from fitted models
#' @description Builds a piecewiseSEM formula from one or more fitted models.
#' The first model defines the main response. Additional models may be
#' supplied as a list. Predictor levels control the order of indirect paths.
#'
#' @param fit A fitted `lmer` model, or a list whose first element is a fitted
#' `lmer` model.
#' @param sem.f Formula-building mode, optionally including step or full .
#' @param sem.data Data used to fit the SEM.
#' @param x.level1 Predictors in the first path level.
#' @param x.level2 Predictors in the second path level.
#' @param x.level3 Predictors in the third path level.
#' @param corr Whether to include correlations among predictors in the same
#' level.
#' @return A list containing a formula string and model data.
#'
#' @examples
#' \dontrun{
#'  formula <- y_psem_formula(fit, x.level1 = "treatment")
#' }
#'
#' @family piecewiseSEM
#' @export
#'

y_psem_formula <- function(fit, sem.f = '', sem.data = NULL, x.level1 = NULL, x.level2 = NULL, x.level3 = NULL, corr = TRUE){
      # sem.f: way to setup model, see y_fit_step, using grepl to identity;
      # 'step_.1', 'step_.05', 'step_aicc' is OK; but 'step_.._full', add all links
      # x.level, distinguish the levels of x; 1 is must, all x but 1,2,3 are 4, max is 4
      # 3 maybe enough, e.g., 1-3 could be vars of plant, soil, microbe, and the y could be EMF
      # corr = T: all x in each level will be linked by %~~%

      # noting that only 'lmer' using for the sub_fit, but you can provide a fit.list include all fits and,
      # add the cov_fit by running this function; [1] is the main fit, others is the sub_fit, totally as you supply.

      if(is.list(fit)){fit[-1] -> sub_fit.0; fit[[1]] -> fit} else{sub_fit.0 = NULL}

      if(is.null(sem.data)){lme4::getData(fit) -> df} else{sem.data -> df}
      getCall(fit) %>% deparse1() %>%  str_replace(., "data = .*,", "data = df,") %>% paste0(",") -> main_fit
      reformulas::findbars(stats::formula(fit)) %>% as.character() %>% paste0(" + (", .,")") -> random

      all.vars(stats::formula(fit, fixed.only = T))[-1] -> vars.x
      vars.x[order(factor(vars.x, levels = vars.x))] -> x

      x[x %in% x.level1] -> x.sub1 # e.g., plant
      x[x %in% x.level2] -> x.sub2 # e.g., soil
      x[x %in% x.level3]  -> x.sub3 # e.g., microbe
      x[!x %in% c(x.sub1, x.sub2, x.sub3)] -> x.sub4 # others

        if(length(x.sub1) == 0){x.level1[[1]] -> x.sub1; cat("\n0_o| there is no x1 in fit & x.level1 is needed\n----------------\n")}
        if(length(x.sub2) == 0){NULL -> x.sub2}; if(length(x.sub3) == 0){NULL -> x.sub3}; if(length(x.sub4) == 0){NULL -> x.sub4}

        # default, only x in main fit
          x.sub2 %>% imap(~{paste0("lmer(", .x, " ~ ", paste0(x.sub1, collapse = " + "), random, ", data = df),")}) -> sub_fit2
          x.sub3 %>% imap(~{paste0("lmer(", .x, " ~ ", paste0(c(x.sub1, x.sub2), collapse = " + "), random, ", data = df),")}) -> sub_fit3
          x.sub4 %>% imap(~{paste0("lmer(", .x, " ~ ", paste0(c(x.sub1, x.sub2, x.sub3), collapse = " + "), random, ", data = df),")}) -> sub_fit4
          c(sub_fit2, sub_fit3, sub_fit4) %>% unlist() -> sub_fit

        if(grepl('step', sem.f)){ifelse(grepl('.1', sem.f) == FALSE & grepl('.05', sem.f) == FALSE, 'aicc', sem.f) -> step; # aicc is default
          suppressMessages(suppressWarnings(x.sub2 %>% imap(~{paste0(.x, " ~ ", paste0(x.level1, collapse = " + "), random)}) %>%
            set_names(x.sub2) %>% imap(~{lmer(formula = .x, data = df) %>% y_fit_step(., step = step, output = "best_fit")}) -> sub_fit2))
          suppressMessages(suppressWarnings(x.sub3 %>% imap(~{paste0(.x, " ~ ", paste0(c(x.level1, x.sub2), collapse = " + "), random)}) %>%
            set_names(x.sub3) %>% imap(~{lmer(formula = .x, data = df) %>% y_fit_step(., step = step, output = "best_fit")}) -> sub_fit3))
          suppressMessages(suppressWarnings(x.sub4 %>% imap(~{paste0(.x, " ~ ", paste0(c(x.level1, x.sub2, x.sub3), collapse = " + "), random)}) %>%
            set_names(x.sub4) %>% imap(~{lmer(formula = .x, data = df) %>% y_fit_step(., step = step, output = "best_fit")}) -> sub_fit4))

          c(sub_fit2, sub_fit3, sub_fit4) %>% imap(~{getCall(.x) %>% deparse1() %>% paste0(",")}) %>% unlist() -> sub_fit
          c(sub_fit2, sub_fit3) %>% imap(~{all.vars(stats::formula(.x, fixed.only = T))[-1]}) %>% unlist() %>% unique() %>% .[!. %in% x.sub2] -> x.sub1.new
          c(x.sub1, x.sub1.new[!x.sub1.new %in% x.sub1]) %>% unique() -> x.sub1; # add the new x.sub1, generated in step

          if(grepl('full', sem.f)){ # add all links

            y_fit_update(fit, str = paste0('+', paste0(x.sub1, collapse = ' + '))) %>%
              getCall() %>% deparse1() %>%  str_replace(., "data = .*,", "data = df,") %>% paste0(",") -> main_fit

            x.sub2 %>% imap(~{paste0("lmer(", .x, " ~ ", paste0(x.sub1, collapse = " + "), random, ", data = df),")}) -> sub_fit2
            x.sub3 %>% imap(~{paste0("lmer(", .x, " ~ ", paste0(c(x.sub1, x.sub2), collapse = " + "), random, ", data = df),")}) -> sub_fit3
            x.sub4 %>% imap(~{paste0("lmer(", .x, " ~ ", paste0(c(x.sub1, x.sub2, x.sub3), collapse = " + "), random, ", data = df),")}) -> sub_fit4
            c(sub_fit2, sub_fit3, sub_fit4) %>% unlist() -> sub_fit}
          }

        if(corr == TRUE){
          if(length(x.sub1) > 1){combn(x.sub1, 2, FUN = function(pair){paste0(pair[1], " %~~% ", pair[2], ",")}) -> cov_fit1} else{NULL -> cov_fit1}
          if(length(x.sub2) > 1){combn(x.sub2, 2, FUN = function(pair){paste0(pair[1], " %~~% ", pair[2], ",")}) -> cov_fit2} else{NULL -> cov_fit2}
          if(length(x.sub3) > 1){combn(x.sub3, 2, FUN = function(pair){paste0(pair[1], " %~~% ", pair[2], ",")}) -> cov_fit3} else{NULL -> cov_fit3}
          if(length(x.sub4) > 1){combn(x.sub4, 2, FUN = function(pair){paste0(pair[1], " %~~% ", pair[2], ",")}) -> cov_fit4} else{NULL -> cov_fit4}
          c(cov_fit1, cov_fit2, cov_fit3, cov_fit4) -> cov_fit} else{NULL -> cov_fit}

      if(is.list(sub_fit.0)){sub_fit.0 %>% imap(~{getCall(.x) %>% deparse1() %>%  str_replace(., "data = .*,", "data = df,") %>% paste0(",") }) %>% unlist -> sub_fit}
      paste0("piecewiseSEM::psem(\n", paste0(c(main_fit, sub_fit, cov_fit), collapse = "\n") %>% sub(",[^,]*$", "", .), "\n)") -> sem.f; cat(sem.f, '\n\n')
      return(invisible(list("sem.formula" = sem.f, "sem.data" = df)))
    }



#' @title Update a piecewiseSEM formula
#' @description Updates selected paths in an existing piecewiseSEM formula.
#'
#' @param sem.formula A piecewiseSEM formula string or a list returned by
#' y_psem_formula() .
#' @param sem.data Data used to fit the SEM.
#' @param id Identifiers of the paths to remove.
#' @param str Formula update text, such as + x1 - x2 .
#' @return A list containing an updated formula and model data.
#'
#' @family piecewiseSEM
#' @export
#'

y_psem_formula_update <- function(sem.formula, sem.data, id = NULL, str = NULL){
      # id = c(3,2,6,17);  str = c("-x1+x2", '+', '', '&'); see y_fit_update

      if(is.list(sem.formula)){sem.formula$sem.data -> sem.data; sem.formula$sem.formula -> sem.formula}
      if(is.null(str)){rep('', time = length(id)) -> str}; sem.data -> df

      sem.formula %>% gsub("\\s*\\n\\s*", "\n", .) -> sem.formula.1
      strsplit(sem.formula.1, "\n")[[1]] %>% .[-length(.)] %>% .[-1] -> sem.f_id
      sem.f_id[!grepl('%~~%', sem.f_id)] -> sem.fit;

      id[id > length(sem.f_id)] -> id.up; id[id <= length(sem.fit)] -> id.fit; id[!id %in% c(id.up, id.fit)] -> id.corr

      for(i in id.corr){'' -> sem.f_id[i]}
      for(j in id.up){c(sem.f_id, str[match(j, id)]) -> sem.f_id}
      for(k in id.fit){if(str[match(k, id)] == ''){'' -> sem.f_id[k]} else{
        eval(parse(text = sem.f_id[k] %>% str_replace(',$', ''))) %>% y_fit_update(str = str[match(k, id)]) %>%
          getCall() %>% deparse1() %>% paste0(., ',') -> fit.f; ifelse(grepl('\\+', fit.f), fit.f, '') -> sem.f_id[k]}
      }

      # delete corr not in fit
      sem.f_id[grepl('%~~%', sem.f_id)] %>% strsplit(' %~~% ') %>% unlist() %>% str_replace(',', '') %>% unique() -> vars.corr
      vars.corr %>% imap(~{grepl(.x, sem.f_id[!grepl('%~~%', sem.f_id)] %>% paste0(collapse = '\n'))}) %>% unlist() -> vars.corr.id
      vars.corr[vars.corr.id == FALSE] -> A; if(length(A) > 0){A %>% c(paste0(., ' %~~%'), paste0('%~~% ', .)) %>%
        imap(~{sem.f_id[grepl(.x, sem.f_id)]}) %>% unlist() -> corr.no; sem.f_id[!sem.f_id %in% corr.no] -> sem.f_id}

      sem.f_id[sem.f_id != ''] %>% paste0(., ',') %>% str_replace_all(',,', ',') -> sem.f_id; cat(' id', '  str\n');
      paste0("piecewiseSEM::psem(\n",  paste0(sem.f_id, collapse = '\n') %>% sub(",[^,]*$", "", .), "\n)") -> sem.f; print(sem.f_id)
    return(invisible(list("sem.formula" = sem.f, "sem.data" = df)))
    }



#' @title Prune paths from a piecewiseSEM result
#' @description Removes paths by identifier, variable, p-value or step count.
#' Paths can be removed by identifier, variable, highest p-value, or a maximum
#' number of steps.
#'
#' @param psem.res A result returned by y_psem_boot() .
#' @param id Identifiers of the paths to remove.
#' @param vars Variable names whose paths should be removed.
#' @param step.n Maximum number of highest-p-value paths to remove.
#' @param step.p Prune paths until all p-values are below this threshold.
#' @param corr.p Prune correlations above this p-value threshold.
#' @param nrep Number of bootstrap replicates; zero fits the model without
#' bootstrapping.
#' @param show.p Significance threshold used for display.
#' @param show.est Effect size to plot: standardized or unstandardized.
#' @return The pruned piecewiseSEM result, invisibly.
#'
#' @family piecewiseSEM
#' @export
#'

y_psem_step <- function(psem.res, id = NULL, vars = NULL, step.n = NULL, step.p = NULL, corr.p = NULL,
                        nrep = 0, show.p = 0.1, show.est = "std"){
      # id = 2: delete the 2 in psem.coef (coefficients)
      # vars = c("x1"..): delete the all x1 and x.. in psem model
      # step.n = 3: delete the max P.value in psem.coef, 1 by 1, until 3
      # step.p = .1: delete the max P.value in psem.coef, 1 by 1, until P <= .1
      # corr.p = .1: delete corr that all P.value > .1 in psem.coef, 1 times but all

      yy_check_sem_dependencies()

      if(length(psem.res) == 2){capture.output(y_psem_boot(psem.res, nrep = nrep, show.p = show.p, show.est = show.est) -> psem.res)}

      psem_update.f.no <- function(psem.res, id.no = NULL){

        psem.res$sem.results$coefficients %>% data.frame() %>% rownames_to_column('id') -> psem.coef
        psem.coef$Response %>% unique() %>% .[!grepl('~~', .)] -> y;  psem.coef$Predictor %>% unique() %>% .[!grepl('~~', .)] -> x

        psem.coef %>% dplyr::filter(!grepl('~~', Response)) %>% mutate(id.fit = match(Response, y)) -> psem.fit
        psem.fit %>% dplyr::filter(id %in% id.no) %>% group_by(id.fit) %>% summarise(str = paste0('-', paste0(Predictor, collapse = '-'))) -> id.no.fit

        if(nrow(psem.coef) > nrow(psem.fit)){
          psem.coef %>% dplyr::filter(grepl('~~', Response)) %>% mutate(id.corr = c(1:length(.$Response)) + length(y)) -> psem.corr
          psem.corr %>% dplyr::filter(id %in% id.no) %>% .$id.corr -> id.no.corr} else{id.no.corr = NULL}

        suppressMessages(suppressWarnings(
          y_psem_formula_update(psem.res, id = c(id.no.fit$id.fit, id.no.corr), str = id.no.fit$str) %>%
          y_psem_boot(nrep = nrep, show.p = show.p, show.est = show.est) %>% return()))
      }

     if(is.null(id) == FALSE){psem_update.f.no(psem.res, id.no = id) -> psem.res;
       cat('\n0_o| done, delete id =', id, '(in psem coef..)\n---------------------------------------------------^\n')}

     if(is.null(vars) == FALSE){
       psem.res$sem.results$coefficients %>% data.frame() %>% rownames_to_column('id') -> psem.coef
       vars %>% imap(~{c(psem.coef$id[grepl(.x, psem.coef$Response)], psem.coef$id[grepl(.x, psem.coef$Predictor)])}) %>% unlist() -> id
       psem_update.f.no(psem.res, id.no = id) -> psem.res;
       cat('\n0_o| done, delete vars =', vars, '(in psem & all)\n---------------------------------------------------^\n')}

     if(is.null(corr.p) == FALSE){
       psem.res$sem.results$coefficients %>% data.frame() %>% rownames_to_column('id') -> psem.coef
       psem.coef %>% dplyr::filter(P.Value > corr.p & grepl('~~', Predictor)) %>% .$id %>% as.numeric() %>% psem_update.f.no(psem.res, id.no = .) -> psem.res
       cat('\n0_o| done, delete all corr that P.Value >', corr.p, '\n---------------------------------------------------^\n')}

     if(is.null(step.p) == FALSE){i = 0;
       while(psem.res$sem.results$coefficients$P.Value %>% max() > step.p){
         psem.res$sem.results$coefficients %>% data.frame() %>% rownames_to_column('id') -> psem.coef
         psem.coef %>% dplyr::filter(P.Value == max(P.Value)) %>% .$id %>% as.numeric() %>% psem_update.f.no(psem.res, id.no = .) -> psem.res
         i + 1 -> i; cat('\n0_o| done', i, 'times, delete P.Value = ', max(psem.coef$P.Value), '\n---------------------------------------------------^\n')}
       }

     if(is.null(step.n) == FALSE){if(step.n > 0){for(i in 1:step.n){
       psem.res$sem.results$coefficients %>% data.frame() %>% rownames_to_column('id') -> psem.coef
       psem.coef %>% dplyr::filter(P.Value == max(P.Value)) %>% .$id %>% as.numeric() %>% psem_update.f.no(psem.res, id.no = .) -> psem.res
       cat('\n0_o| done', i, 'times, delete P.Value = ', max(psem.coef$P.Value), '\n---------------------------------------------------^\n')}}
      }

    invisible(psem.res) %>% return()
  }



#' @title Fit and bootstrap a piecewiseSEM model
#' @description Fits a piecewiseSEM model with optional bootstrap confidence
#' intervals. Bootstrap failures are counted and reported. Setting nrep = 0
#' returns the direct piecewiseSEM fit without bootstrap confidence intervals.
#' Bootstrap failures are counted and reported. Setting nrep = 0 returns
#'
#' @param sem.formula A piecewiseSEM formula string or a list returned by
#' y_psem_formula() .
#' @param sem.data Data used to fit the SEM.
#' @param nrep Number of bootstrap replicates; zero fits the model without
#' bootstrapping.
#' @param show.p Significance threshold used for display.
#' @param show.est Effect size to plot: standardized or unstandardized.
#' @param output Output mode; see the function description.
#' @return A piecewiseSEM fit or a list containing fit results and plot data.
#'
#' @family piecewiseSEM
#' @export
#'

y_psem_boot <- function(sem.formula, sem.data, nrep = 0, show.p = 0.1, show.est = "std", output = ""){
      # nrep = 0, return the raw summary of psem with raw data
      # show.p = 0.1, alpha = 0.1, arrows with p < 0.1 in sem-plot are solid

      yy_check_sem_dependencies()

      if(is.list(sem.formula)){sem.formula$sem.data -> sem.data; sem.formula$sem.formula -> sem.formula}

      sem.data -> df; eval(parse(text = sem.formula)) -> sem.fit0
      sem.fit0 %>% summary(conserve = TRUE) -> sem.results0

      if(nrep == 0){

        sem.results0 -> sem.results
        } else{
            sem_bootstrap <- function(sem.data, sem.formula){
              sampleid <- sample(nrow(sem.data), replace = TRUE)
              df <- sem.data[sampleid, ]

              tryCatch({
                sem_model <- eval(parse(text = sem.formula))
                sem.summary <-  summary(sem_model, conserve = TRUE)
                return(list(sem.summary$Cstat, sem.summary$R2, sem.summary$coefficients))
              }, error = function(e) {return(NULL)})
            }

            # calculating median and 95% Ci
            m.ci <- function(data_col){
              median <- median(data_col, na.rm = TRUE) %>% round(3)
              low.ci <- stats::quantile(data_col, .025, na.rm = TRUE) %>% round(3)
              up.ci <- stats::quantile(data_col, .975, na.rm = TRUE) %>% round(3)
              cbind(median, low.ci, up.ci) %>% data.frame(row.names = "1") %>% return()
              }

        # run bootstrapping
          replicate(nrep, sem_bootstrap(sem.data, sem.formula)) -> boot_results

          boot_results[seq(1, length = nrep, by = 3)] %>% bind_rows() %>% apply(., 2, m.ci) %>% do.call(rbind, .) -> sem_Fisher.C_ci

          boot_results[seq(2, length = nrep, by = 3)] %>% bind_rows() %>%
            rename('M' = 'Marginal', 'C' = 'Conditional') %>% group_by(Response) %>% reframe(
              Marginal = m.ci(M)[[1]], Marginal.low = m.ci(M)[[2]], Marginal.up = m.ci(M)[[3]],
              Conditional = m.ci(C)[[1]], Conditional.low = m.ci(C)[[2]], Conditional.up = m.ci(C)[[3]]) -> sem_R.square_ci

          boot_results[seq(3, length = nrep, by = 3)] %>%  bind_rows() %>% dplyr::filter(!grepl("=", Predictor)) %>%
            rename('Est' = 'Estimate', 'Std.E' = 'Std.Estimate', 'P' = 'P.Value') %>% group_by(Response, Predictor) %>% reframe(
              Estimate = m.ci(as.numeric(Est))[[1]], Estimate.low = m.ci(as.numeric(Est))[[2]], Estimate.up = m.ci(as.numeric(Est))[[3]],
              Std.Estimate = m.ci(as.numeric(Std.E))[[1]], Std.Estimate.low = m.ci(as.numeric(Std.E))[[2]], Std.Estimate.up = m.ci(as.numeric(Std.E))[[3]],
              P.Value = m.ci(P)[[1]], P.Value.low = m.ci(P)[[2]], P.Value.up = m.ci(P)[[3]]) -> sem_coefficients

          # reorder
            sem_R.square_ci %>% merge(sem.results0$R2 %>% select(Response), ., by = c("Response"), sort = FALSE) -> sem_R.square_ci
            sem_coefficients %>% merge(sem.results0$coefficients %>% select(Response, Predictor), ., by = c("Response", "Predictor"), sort = FALSE)-> sem_coefficients
            list(Cstat = sem_Fisher.C_ci$median,  coefficients = sem_coefficients, R2 = sem_R.square_ci) -> sem.results
        }

          # std.Estimate if it is NA
          strsplit(sem.formula, "\n")[[1]] %>% .[-length(.)] %>% .[-1] -> sem.f_id
          sem.results$coefficients$Response %>% unique() %>% .[!grepl('~~', .)] -> y; sem.results$coefficients$Predictor %>% unique() %>% .[!grepl('~~', .)] -> x
          df %>% select(all_of(c(x, y))) %>%  apply(2, function(x){stats::sd(x)}) -> sd.xy; df %>% select(all_of(c(y))) %>%  apply(2, function(x){mean(x)}) -> mean.y

          sem.f_id[1:length(y)] %>% imap(~{ifelse(grepl('family', .x) == TRUE & grepl('log', .x) == TRUE, TRUE, FALSE)}) %>% set_names(y) %>% unlist() -> link.log

          sem.results$coefficients %>% data.frame()  %>%
            mutate(sd.x = sd.xy[match(Predictor, names(sd.xy))], sd.y = sd.xy[match(Response, names(sd.xy))],
                   mean.y = mean.y[match(Response, names(mean.y))], link.log = link.log[match(Response, names(link.log))]) %>%
             mutate(std.Est = ifelse(link.log == TRUE, Estimate * mean.y * sd.x/sd.y, Estimate * sd.x/sd.y)) %>%
             mutate(Std.Estimate = ifelse(Std.Estimate == '-' | is.na(Std.Estimate), round(std.Est, 3), Std.Estimate)) %>%
           select(!c(sd.x, sd.y, mean.y, link.log, std.Est)) -> sem.results$coefficients

          if("Estimate.low" %in% names(sem.results$coefficients)){sem.results$coefficients %>%
              mutate(Std.Estimate.low = ifelse(is.na(Std.Estimate.low), round(Std.Estimate/Estimate*Estimate.low, 3), Std.Estimate.low),
                     Std.Estimate.up = ifelse(is.na(Std.Estimate.up), round(Std.Estimate/Estimate*Estimate.up, 3), Std.Estimate.up)) -> sem.results$coefficients}

         # SEM plot
         if(show.est == "ustd"){"Estimate" -> effect.size} else{"Std.Estimate" -> effect.size}

            sem.results$Cstat %>% round(3) -> Fisher.C
            sem.results$R2 %>% select(nodes = Response, m.R = Marginal) -> RR
            sem.results$coefficients %>% select(nodes = Response, Est = all_of(effect.size), P.Value) %>%
              mutate(P.sig = cut(P.Value, breaks = c(-Inf, .01, .05, Inf), labels = c("**", "* ", "  "))) -> PP

          plot(sem.fit0, layout = "tree", alpha = show.p, show = "ustd", digits = 2, return = TRUE) -> sem_graph

          sem_graph$global_attrs <-
            rbind(sem_graph$global_attrs,
                  data.frame(attr = "label", value = paste0(c("Fisher C = ", "DF = ", "P-value = "), Fisher.C, collapse = ", "), attr_type = "graph"),
                  data.frame(attr = c("labelloc", "fontsize", "fontname"), value = c("b", "10", "Times-Italic"), attr_type = "graph"))

          sem_graph$nodes_df <- sem_graph$nodes_df %>%
            mutate(fontsize = 7, fontcolor = "black", width = .8, height = .42, color = rgb(0, 0, 0, .5), fillcolor = rgb(0, 0, 0, alpha = 0)) %>%
            mutate(m.R = ifelse(nodes %in% RR$nodes, RR$m.R[match(nodes, RR$nodes)], NA)) %>%
            mutate(label = ifelse(!is.na(m.R), paste0(label, "\nm.R^2 = ", round(m.R, 2)), label))

          sem_graph$edges_df <- sem_graph$edges_df %>%
            mutate(fontsize = 6, arrowsize = .5, label = round(as.numeric(PP$Est), 2), P.sig = PP$P.sig) %>%
            mutate(color = ifelse(as.numeric(label) >= 0, "steelblue2", "deeppink2"), fontcolor = "black",
                   style = ifelse(P.sig == "  ", style, "solid"), label = paste0(label, P.sig))

          DiagrammeR::render_graph(sem_graph, layout = "tree") -> sem.plot; print(sem.plot); print(sem.results$coefficients)

      # output
      if(output == "sem.fit"){sem.fit0 %>% return()} else{
        invisible(list('sem.formula' = sem.formula, 'sem.data' = df, "sem.results" = sem.results, "sem.plot" = sem.plot)) %>% return()}
     }
