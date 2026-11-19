
# 0_o .....
# -------------------------------------------------------------------------------------------------------------------- Read ME #
# **there is no env: unconstrained ordination
# PCA (if no too much 0 in community) else CA, DCA
# PCoA, NMDS (if based on Bray-Curtis and need to PERMANOVA)
#
# ** there is env: constrained ordination:
# db-RDA (if based on Bray-Curtis and need to PERMANOVA)
# RDA, CCA (if need to anova.cca & rdacca.hp; 'decorana' if DCA1 <= 4 RDA else CCA)
#
# 数据标准化（decostand & scale）; 环境因子全部建议scale消除量纲
# PCA 必须 优先 decostand(spe, "standardize") 或 scale = TRUE（两者相同）
# RDA 必须 优先 decostand(spe, "hellinger") 次 scale = TRUE
# NMDS 可 distance = 'bray','jaccard' 时不需要, 欧式距离(euclidean)时必须 scale(spe)
# PCoA, db-RDA 可 vegdist('bray') 时不需要，欧式距离(euclidean)时必须 scale(spe)
# CCA，DCA，CA 不可 必须原始物种矩阵

# ---------------------------------------------------------------------------------------------------------------------------- #
# # 群落β-dissimilarity，基于Bray-Curtis的相关检验      # spe需为raw 适用于NMDS，PCoA，db-RDA (Must based on Bray-Curtis)
# adonis2(spe ~ trt, method = "bray")                   # 基于距离的方差分解（看组间是否不同, 输出 R² 和 p 值）
# pairwise.adonis(spe, trt, sim.method = "bray")        # adonis2之后的组间比较（两两分组比较处理间的群落差异性）
# anosim(spe, trt, distance = "bray")                   # 基于秩的检验（组内是否比组间更相似, R > 0.7, 强分离; R < 0.25 弱重叠）
#
# meandist(vegdist(spe), trt, distance = "bray")        # 组间平均距离矩阵 (组内、组间相异度dissimilarity)
# simper(spe, trt, distance = "bray") -> sim            # contributions of each species to dissimilarity
# multipatt(spe, trt, duleg = TRUE)                     # 寻找指示物种

# ---------------------------------------------------------------------------------------------------------------------------- #
# # 判断环境因子是否显著解释群落差异                    # 适用于RDA，CCA，db-RDA
# decostand(spe, "hellinger") -> spe.hel                # hellinger transform to avoid "bow(弓) effect"
# rda(spe.hel ~., env, scale = scale) -> rda_model      # 检验环境因子总效应的R方和P值，通用于RDA，CCA and their partial patterns
# RsquareAdj(RDA); anova.cca(rda_model)                 # margin'检验各因子的净效应; 既多因子共存时，因子有其独特作用吗？受共线性影响
# anova.cca(RDA, by = 'margin')                         # 因子单独与群落结构有关吗？不控制其它因子，类似相关性
# envfit(RDA, env)                                      # envfit & anova.acc 同时显著 --> 核心因子(独立贡献显著)

# varpart(spe, env[1:3], env[4:6], transfo = 'hel')         # chisquare = FALSE is RDA; chisquare = T is CCA                                                 #
# rdacca.hp(spe.hel, env, method = "RDA", type = "adjR2")   # 给出各个环境因子的相对贡献，即使共线变量也可适用(不高于9变量相对可靠)
# rdacca.hp(spe.hel, list(A =  env[1:3], B = env[4:6]),...) # 为环境因子分组给出贡献度Proportion of explained variation (%)
# ---------------------------------------------------------------------------------------------------------------------------- #

# 0_0...test
# library(vegan)
# data("varespec"); data("varechem"); varespec -> spe; varechem[1:9] -> env;
# rep(c("CK", "A", "B", "C"), each = 6) %>% factor(level = c("CK", "A", "B", "C")) -> trt
# paste0(trt, rep(1:6)) -> row.names(spe) -> row.names(env)
# # test(spe, env, trt三者数据行名称必须对应样品编号)
#
# row.names(spe) -> A; paste0(A[grepl('CK', A)], ' vs ', A[!grepl('CK', A)]) -> pairs.ck
#
# y_price_part(spe[1:3,]) %>% y_price_part_plot()          # 价格方程
# ytips:::yy_adonis_summary(spe, show.p.pairs = '')        # β-dissimilarity
# y_beta_part(spe, pairs.id = pairs.ck)                    # 分解 turn-over and nest
# y_beta_simper(spe, pairs.id = pairs.ck)                  # 识别物种贡献
# y_beta_indispecies(spe)                                  # 识别指示物种
# y_pcoa_plot(spe, trt, env)                               # PCoA & adonis & envfit
# y_rda_plot(spe, trt, env, show.arrow = 'biplot')         # RDA & anova.cca
# y_rdacca_hp(spe, env, env.list = list(a = 1:3, b = 4:7, c = 8:9))  # 分解环境因子贡献
# ---------------------------------------------------------------------------------------------------------------------------- #



#' @title Convert a phyloseq object to community matrices
#' @description Converts phyloseq abundance, sample and taxonomy tables.
#'
#' @param ps A phyloseq object.
#' @param i.trt Column index or name defining treatment groups.
#' @param i.env Optional environmental-variable columns.
#' @param output Output mode; see the function description.
#' @return A community matrix or a list of converted data.
#'
#' @family ordination
#' @export
#'

y_ps_to_spe <- function(ps, i.trt = 1, i.env = NULL, output = 'spe'){

    yy_check_phyloseq()

    ps %>% phyloseq::otu_table() %>% as.data.frame() -> otu.df;
    ps %>% phyloseq::sample_data() %>% as.data.frame() -> sample.df;
    ps %>% phyloseq::tax_table() %>% as.data.frame() -> taxa.df

  otu.df %>% t() %>% as.data.frame() -> spe; sample.df[[i.trt]] -> trt;
  if(!is.null(i.env)){sample.df[i.env] -> env} else{NULL -> env}

  if(output == 'spe'){return(spe)} else{return(list('spe' = spe, 'trt' = trt, 'env' = env))}
}



#' @title Calculate alpha-diversity metrics
#' @description Calculates abundance, richness, Shannon, Simpson and Pielou
#' metrics.
#'
#' @param spe A community matrix with samples in rows.
#' @return A data frame of alpha-diversity metrics.
#'
#' @family ordination
#' @export
#'

y_alpha_diversity <- function(spe){
  abundance = rowSums(spe)
  richness = vegan::specnumber(spe, MARGIN = 1)                     # r = n/sum(n);
  shannon = vegan::diversity(spe, index = "shannon", MARGIN = 1)    # shannon = -sum(r*log(r))
  simpson = vegan::diversity(spe, index = "simpson", MARGIN = 1)    # simpson = 1 - sum((r)^2)
  pielou = shannon/log(richness)
  data.frame(richness, abundance, shannon, simpson, pielou) %>% round(4) %>% rownames_to_column("sample") %>%
  return()
}



#' @title Decompose community change with the Price equation
#' @description Calculates Price-equation components for sample pairs. Splits
#' community change into richness, composition, colonization, extinction and
#' resident-species effects.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param func.inrr Optional function expression used to transform a sample
#' pair.
#' @param pairs.id Optional sample-pair identifiers separated by ' vs ' .
#' @return A data frame of Price-equation components.
#'
#' @examples
#' spe <- matrix(c(1, 2, 0, 1, 3, 1), nrow = 3,
#'   dimnames = list(c("CK1", "A1", "B1"), c("s1", "s2")))
#' price <- y_price_part(spe)
#' y_price_part_plot(price)
#'
#' @family ordination
#' @export
#'

y_price_part <- function(spe, trt, func.inrr = NULL, pairs.id = NULL){
  # func.inrr 为适用于spe.2的function，计算出功能值后用于分解
  # func.inrr = 'y_beta_simper()' β不相似性中的每个物种贡献度，分解β-diversity
  # func.inrr = apply(., 2, function(x) log(x[2]/x[1])); 效应值（不可包含0值）
  # pairs.id = c('CK1 vs A1', ...); pairs.ck 提取block条件下处理与CK的配对

  if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}; data.frame(row.names(spe), trt) -> Trt
  if(is.null(pairs.id)){combn(row.names(spe), 2, simplify = FALSE) -> pairs.id} else{str_split(pairs.id, ' vs ') -> pairs.id}

  price_part_spe.2 <- function(spe.2, func.inrr = NULL){ spe.2 %>% t() %>% as.data.frame() %>%
      rename_with(~c('func.x', 'func.y')) %>% rownames_to_column('species') %>% dplyr::filter(func.x > 0 | func.y > 0) %>%
      mutate(wvec = ifelse(func.x > 0 & func.y > 0, 1, 0), xvec = ifelse(func.x > 0, 1, 0), yvec = ifelse(func.y > 0, 1, 0)) -> comm

    x.rich <- sum(comm$xvec);  y.rich <- sum(comm$yvec); c.rich <- sum(comm$wvec); x.func <- sum(comm$func.x); y.func <- sum(comm$func.y)

    if(is.null(func.inrr)){
      SRE.L <- (c.rich - x.rich)/x.rich * x.func
      SRE.G <- (y.rich - c.rich)/y.rich * y.func
      SIE.L <- -sum(comm$func.x[comm$xvec == 1 & comm$yvec == 0]) - SRE.L
      SIE.G <-  sum(comm$func.y[comm$xvec == 0 & comm$yvec == 1]) - SRE.G
      CDE <- sum(comm$func.y[comm$wvec == 1] - comm$func.x[comm$wvec == 1])

    } else{
      eval(parse(text = paste0('spe.2 %>% ', func.inrr,' -> func.inrr')))
      func.inrr[,1:2] %>% rename_with(~c('species', 'func.z')) %>% merge(comm, by = 'species', sort = FALSE) -> comm
      z.func <- sum(comm$func.z); x.func <- 0; y.func <- z.func

      SRE.L <- (x.rich - c.rich)/nrow(comm) * z.func
      SRE.G <- (y.rich - c.rich)/nrow(comm) * z.func
      SIE.L <-  sum(comm$func.z[comm$xvec == 1 & comm$yvec == 0]) - SRE.L
      SIE.G <-  sum(comm$func.z[comm$xvec == 0 & comm$yvec == 1]) - SRE.G
      CDE <- sum(comm$func.z[comm$wvec == 1]) }

    SL <- SRE.L + SIE.L; SG <- SRE.G + SIE.G; SR <- SRE.L + SRE.G; CE <- SIE.L + SIE.G + CDE
    pp <- c(SRE.L, SRE.G, SIE.L, SIE.G, CDE, SL, SG, SR, CE, x.func, y.func, x.rich, y.rich, c.rich)
    names(pp) <- c("SRE.L", "SRE.G", "SIE.L", "SIE.G", "CDE", "SL", "SG", "SR", "CE", "x.func", "y.func", "x.rich", "y.rich", "c.rich")
    bind_rows(pp) %>% mutate(pairs.0 = paste0(row.names(spe.2), collapse = ' vs '), .before = 1) %>% return() }

   map(pairs.id, ~spe[.x, ]) %>% map(~price_part_spe.2(.x, func.inrr = func.inrr) ) %>%
     bind_rows() %>% separate(pairs.0, into = c('id', 'id.'), sep = ' vs ', remove = FALSE) %>%
     mutate(pairs = paste0(Trt[match(id, Trt[,1]), 2], ' vs ', Trt[match(id., Trt[,1]), 2]), .after = 1) %>%
     select(!c(id, id.)) %>% return()
 }



#' @title Plot Price-equation components
#' @description Plots two-, three- or five-component Price-equation paths.
#'
#' @param price.df Data returned by y_price_part() .
#' @param price Number of Price-equation components: 2, 3, or 5.
#' @param cols Optional color values.
#' @param xlim Optional x-axis limits.
#' @param show Whether to show all paths or only the mean path.
#' @return A ggplot object.
#'
#' @family ordination
#' @export
#'

y_price_part_plot <- function(price.df, price = 5, cols = NULL, xlim = NULL, show = c('all', 'mean')){

  if(is.null(cols)){c('#008002', '#56F256', '#DB0102', '#FF6500', '#001ECC', '#48C5E9', '#A01AE5') -> cols}
  if(price == 2){cols[1:2] -> cols} else if(price == 3){cols[c(3,5,7)] -> cols} else {cols[3:7] -> cols}
  if(is.null(xlim)){c(min(price.df$c.rich), max(c(price.df$x.rich, price.df$y.rich))) -> xlim}

  price_part_line <- function(cc, price, xlim, ...){

    with(cc, c(y.rich, x.func + SR)) ->  end.SR;
    with(cc, c(c.rich, x.func + SL)) ->  end.SL; with(cc, c(y.rich, end.SL[2] + SG)) ->  end.SG
    with(cc, c(c.rich, x.func + SRE.L)) -> end.SRE.L; with(cc, c(c.rich, end.SRE.L[2] + SIE.L)) ->  end.SIE.L
    with(cc, c(y.rich, end.SIE.L[2] + SRE.G)) ->  end.SRE.G; with(cc, c(y.rich, end.SRE.G[2] + SIE.G)) ->  end.SIE.G

    with(cc, c(x.rich, x.func)) -> start; with(cc, c(y.rich, y.func)) -> end;
    if(price == 2){rbind(start, end.SR, end.SR, end) %>% data.frame() %>% mutate(txt = rep(c('SR', 'TCE'), each = 2)) -> path.df}
    if(price == 3){rbind(start, end.SL, end.SL, end.SG, end.SG, end) %>% data.frame() %>% mutate(txt = rep(c('SL', 'SG', 'CDE'), each = 2)) -> path.df}
    if(price == 5){rbind(start, end.SRE.L, end.SRE.L, end.SIE.L, end.SIE.L, end.SRE.G, end.SRE.G, end.SIE.G, end.SIE.G, end) %>%
                  data.frame() %>% mutate(txt = rep(c('SRE.L', 'SIE.L', 'SRE.G', 'SIE.G', 'CDE'), each = 2)) -> path.df
                  if(cc$SIE.G * cc$CDE < 0){path.df[6:8, 1] <- cc$y.rich -.01 * (xlim[2] - xlim[1])} }

    path.df %>% mutate(txt = factor(txt, levels = unique(.$txt))) -> path.df
    geom_path(data = path.df, aes(x = X1, y = X2, colour = txt, group = 1), arrow = arrow(length = unit(.2,'cm')), ...) }

  price.df -> c.0; c.0 %>% group_by(pairs) %>% reframe(across(where(is.numeric), mean)) -> c.mean

  if(show[1] == 'mean'){ggplot() -> p} else{
    ggplot() + geom_point(data = c.0, aes(x = y.rich, y = y.func), shape = '_', size = 3, alpha = .3) +
    geom_point(data = c.0, aes(x = x.rich, y = x.func), shape = '*', size = 5, alpha = .3)  -> p
    for(j in 1:nrow(c.0)){p <- p + price_part_line(c.0[j, ], price = price, xlim = xlim, alpha = .2, size = .2)} }

  for(i in 1:nrow(c.mean)){p <- p + price_part_line(c.mean[i, ], price = price, xlim = xlim, alpha = 1, size = .5)}

  p + scale_color_manual("Component", values = c(cols)) +
    geom_point(data = c.mean, aes(x = y.rich, y = y.func), shape = '_', size = 5) +
    geom_point(data = c.mean, aes(x = x.rich, y = x.func), shape = '*', size = 8) +
    labs(x = 'Species richness', y = 'Ecosystem function') + theme_classic() + xlim(xlim) %>% return()
}



#' @title Partition beta diversity
#' @description Partitions Sorensen and Bray-Curtis beta diversity. Sorensen
#' dissimilarity is split into turnover and nestedness. Bray-Curtis
#' dissimilarity is split into balanced variation and gradient.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param pairs.id Optional sample-pair identifiers separated by ' vs ' .
#' @param output Output mode; see the function description.
#' @return A ggplot object or a list containing results and plots.
#'
#' @examples
#' \dontrun{
#'  beta <- y_beta_part(spe, trt)
#' }
#'
#' @family ordination
#' @export
#'

# Sorensen indices(βsor): 0, 1数据集（分解为turnover + nestedness = βsor）
# Bray-Curtis indices(dBC): 丰富数据集（分解为balanced variation + abundance gradient = dBC）
# 相同的β-diversity值，e.g., Bray-Curtis indices内部的四个机制可能不同，分为两大类：
# subset patterns(nestedness, gradient) & substitution patterns(turnover, balanced)
# row.names(spe) -> A; paste0(A[grepl('CK', A)], ' vs ', A[!grepl('CK', A)]) -> pairs.ck
y_beta_part <- function(spe, trt, pairs.id = NULL, output = ''){

  yy_check_pkg("betapart")

  if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}; data.frame(row.names(spe), trt) -> Trt
  if(is.null(pairs.id)){combn(row.names(spe), 2, simplify = FALSE) %>% map(~paste0(.x, collapse = ' vs ')) %>% unlist() -> pairs.id}

  spe %>% mutate(across(everything(), ~ifelse(. > 0, 1, 0))) -> spe.01  # spe 必须为丰度原始值，不可r.abun
  betapart::beta.pair(spe.01, index.family = "sorensen") -> soren       # Sorensen 分解
  betapart::beta.pair.abund(spe, index.family = "bray") -> brayc        # Bray-Curtis 分解（mean = meandist）

  c(soren, brayc) %>% map(~{.x %>% as.matrix() %>% as.data.frame() %>% rownames_to_column('id') %>%
     pivot_longer(cols = !id, names_to = 'id.') %>% mutate(pairs.0 = paste0(id, ' vs ', id.), .before = 1) %>%
     mutate(pairs = paste0(Trt[match(id, Trt[,1]), 2], ' vs ', Trt[match(id., Trt[,1]), 2]), .after = 1) %>%
     dplyr::filter(!id == id.) %>% select(!c(id, id.)) }) %>% bind_rows(.id = 'index') %>%

    mutate(index = recode(index, 'beta.sim' = 'Turnover', 'beta.sne' = 'Nested', 'beta.sor' = 'Sorensen',
      'beta.bray.bal' = 'Balance', 'beta.bray.gra' = 'Gradient', 'beta.bray' = "Bray.Curtis")) %>%
      mutate(pairs = factor(pairs, levels = unique(pairs)), index = factor(index, levels = unique(index))) %>%
      dplyr::filter(., pairs.0 %in% pairs.id) -> part.res

  list("Sorensen indices" = c('Turnover', 'Nested'), "Bray-Curtis indices" = c('Balance', 'Gradient')) %>%
    imap(~{part.res %>% dplyr::filter(index %in% .x) %>% group_by(index, pairs) %>% reframe(value = mean(value)) %>%
    ggplot() + geom_bar(aes(x = pairs, y = value, fill = index), stat = "identity", width = .8) + scale_fill_manual(values = c("red", "blue")) +
      theme_classic() + theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = c(.8, .95), legend.background = element_blank()) +
      theme(legend.key.size = unit(.16, 'in')) + ylim(0,1) + labs(fill = '', x = '', y = .y) }) -> p

  if(output == 'all'){return(list('beta_part.res' = part.res, 'plot' = p))} else{p %>% return()}
}



#' @title Estimate species contributions to beta diversity
#' @description Calculates species contributions to Bray-Curtis dissimilarity.
#' Species contributions are aggregated by sample pair and the remaining
#' contributions are grouped as Others.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param pairs.id Optional sample-pair identifiers separated by ' vs ' .
#' @param cols Optional color values.
#' @param ylim Upper y-axis limit.
#' @param show.contri Cumulative contribution retained before grouping the
#' remainder as Others .
#' @param show.adonis Whether to add PERMANOVA results to the plot.
#' @param output Output mode; see the function description.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object or a list containing contribution results.
#'
#' @family ordination
#' @export
#'

y_beta_simper <- function(spe, trt, pairs.id = NULL, cols = NULL, ylim = 1,
                          show.contri = 0.7, show.adonis = TRUE, output = '', ... ){

  if(is.null(cols)){RColorBrewer::brewer.pal(n = 9, "Set1") %>% scales::alpha(., .8) -> cols}
  if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt; factor(trt, levels = unique(trt)) -> trt}; data.frame(row.names(spe), trt) -> Trt
  if(is.null(pairs.id)){combn(row.names(spe), 2, simplify = FALSE) -> pairs.id} else{str_split(pairs.id, ' vs ') -> pairs.id}

   pairs.id %>% map(~spe[.x, ]) %>% map(~apply(.x, 2, function(x) abs(x[2] - x[1])/ sum(rowSums(.x))) %>%
     data.frame(sim = .) %>% rownames_to_column('species') %>% mutate(pairs.0 = paste0(row.names(.x), collapse = ' vs ')) %>%
     mutate(pairs = paste0(Trt[match(row.names(.x), Trt[,1]), 2], collapse = ' vs ')) ) %>% bind_rows() %>% as_tibble() -> simper.res

   # 按pairs求平均: vegan:::simper(spe, group = trt, distance = "bray") 可以给出p值
   simper.res %>% mutate(pairs  = factor(pairs, levels = unique(pairs))) %>% split(f = .$pairs) %>%
     map(~.x %>% group_by(species) %>% reframe(sim.mean = mean(sim), sim.sd = stats::sd(sim)) %>% arrange(desc(sim.mean)) %>%
     mutate(cumsum.p = cumsum(sim.mean)/sum(sim.mean)) %>% mutate(species = ifelse(cumsum.p <= show.contri, species, 'Others')) %>%
     group_by(species) %>% reframe(sim.mean = sum(sim.mean), sim.sd = mean(sim.sd)) ) %>% bind_rows(.id = 'pairs') -> simper.res.mean

   simper.res.mean %>% group_by(species) %>% reframe(i = sum(sim.mean)) %>% arrange(-i) %>% dplyr::filter(species != 'Others') -> rank.species
   simper.res.mean %>% mutate(pairs = factor(pairs, levels = unique(pairs)), species = factor(species, levels = c(rank.species$species, 'Others'))) %>%
     ggplot() + geom_bar(aes(x = pairs, y = sim.mean, fill = species), stat = "identity", width = .8) + scale_fill_manual(values = cols) + theme_classic() +
     theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.key.size = unit(.16, 'in')) + ylim(0, ylim) + labs(fill = 'Species', x = '', y = 'Bray-Curtis indices') -> p

  if(show.adonis){yy_adonis_summary(spe, trt, xylim = .5*ylim, ...) -> adonis.txt;
    p + annotate("text", x = .6, y = Inf, label = adonis.txt[[1]], parse = TRUE, size = 2.8, hjust = 'left', vjust = 1.3) +
        geom_text(data = adonis.txt[[2]], aes(x = .6, y = .48*ylim + y, label = txt), size = 2.5, vjust = "bottom", hjust = 'left') -> p}

  if(output == 'all'){return(list('simper.res' = simper.res, 'simper.res.mean' = simper.res.mean, 'plot' = p))} else{p %>% return()}
}



#' @title Identify indicator species
#' @description Tests indicator species with multipatt.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param output Output mode; see the function description.
#' @return A ggplot object or a list containing indicator-species results.
#'
#' @family ordination
#' @export
#'

y_beta_indispecies <- function(spe, trt, output = ''){
  # $str是IndVal值, $A是平均丰度，$B是频率：str = sqrt(A * B)

  if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt; factor(trt, levels = unique(trt)) -> trt}

  indicspecies::multipatt(spe, trt, control = permute::how(nperm = 999), duleg = TRUE) -> indval
  list(indval$str, indval$sign) %>% map(~{.x %>% as.data.frame() %>% rownames_to_column("species")}) -> indval.res
  indval.res[[1]] %>% pivot_longer(cols = !species, names_to = "trt", values_to = 'Indval') ->indval.str
  indval.res[[2]] %>% pivot_longer(cols = contains('s.'), names_to = "trt") %>% dplyr::filter(value == 1) %>%
    mutate(trt = gsub('s.', '', trt)) %>% select(species, trt, p.value) %>% merge(indval.str, by = c('species', 'trt'), all = TRUE) %>%
    mutate(p.sig = ifelse((p.value > 0.05 | is.na(p.value)), '', '*')) %>% mutate(trt = factor(trt, levels = levels(indval$cluster))) -> indval.df

  ggplot(indval.df, aes(x = trt, y = species)) + geom_tile(color = 0, fill = NA, size = .05) + geom_point(aes(size = Indval, colour = p.sig), pch = 16) +
    scale_color_manual('p-value', values = c("grey", "red"), labels = c("Not sig.", "p < 0.05")) + theme_minimal() + labs(x = '', y = '') +
    theme(axis.text = element_text(colour = 1), axis.text.x = element_text(angle = 45, hjust = 1), legend.key.size = unit(.16, "in")) -> p

  if(output == 'all'){return(list('indval.res' = indval.res[[2]], 'plot' = p))} else{p %>% return()}
}



#' @title Run PERMANOVA and format the result
#' @description Runs PERMANOVA and pairwise post-hoc comparisons. Returns
#' formatted labels for total and pairwise PERMANOVA results.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param show.p Significance threshold used for display.
#' @param show.p.pairs Text filter for pairwise PERMANOVA results.
#' @param xylim Symmetric ordination-axis limits.
#' @param rr.inc Numeric controls for label spacing.
#' @param ... Additional arguments passed to underlying methods.
#' @return A list containing formatted PERMANOVA labels.
#'
#' @family ordination
#' @keywords internal
#'

yy_adonis_summary <- function(spe, trt, show.p = 0.05, show.p.pairs = c('CK', ''), xylim = 1, rr.inc = .05, ...){

  if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}; as.matrix(trt) -> trt
  spe.dist.bray <- vegan::vegdist(spe, method = "bray") # Bray-Curtis本身为相对距离不需要scale

    # PERMANOVA: Permutation multivariate analysis of variance
    vegan::adonis2(spe ~ trt, method = "bray", permutations = 999) -> Adonis   # 检验群落组间差异显著性(R & p)
    vegan::betadisper(spe.dist.bray, group = trt) -> dispersion;                # plot(dispersion, hull = F, ellipse = F)
    vegan::permutest(dispersion); Adonis$R2[1] -> RR; Adonis$`Pr(>F)`[1] -> P   # 检验组内齐性（PERMANOVA需满足组内齐性）
    c(RR, P) %>% map(~ifelse(round(.x, 2) < 0.01, "~`<`~0.01", paste0("~`=`~", round(.x, 2)))) -> RR
    paste0("R^2", RR[1], "~~italic(p)", RR[2]) -> rr.txt

   # pairwise Adonis: adonis2事后检验，两两组间差异
   pairwiseAdonis::pairwise.adonis(x = spe, factors = trt, sim.function = "vegdist", sim.method = "bray", p.adjust.m = "fdr") -> pair.adonis

   pair.adonis %>% select(pairs, p = p.value) %>% dplyr::filter(p < show.p & grepl(show.p.pairs[1], pairs)) -> pair.adonis
   if(nrow(pair.adonis) > 0){pair.adonis %>%
       mutate(p.sig = ifelse(p < 0.01, " **", ifelse(p < 0.05, " * ", "   "))) %>% mutate(txt = paste0(pairs, p.sig)) %>%
       mutate(x = -xylim, y = .98*xylim - rr.inc * c(1:nrow(.))*2*xylim) -> pair.txt} else{data.frame(x = 0, y = 0, txt = ' ') -> pair.txt}

  return(list('rr.txt' = rr.txt, 'pair.adonis.txt' = pair.txt))
}



#' @title Format constrained-ordination tests
#' @description Formats global and marginal tests for a constrained ordination
#' model.
#'
#' @param model.ord An ordination model returned by vegan.
#' @return A list containing formatted test labels and arrow colors.
#'
#' @family ordination
#' @keywords internal
#'

yy_anova_cca_summary <- function(model.ord){

  # Permutation test
  vegan::RsquareAdj(model.ord)$r.squared -> RR; vegan::anova.cca(model.ord)$`Pr(>F)`[1] -> P
  c(RR, P) %>% map(~ifelse(round(.x, 2) < 0.01, "~`<`~0.01", paste0("~`=`~", round(.x, 2)))) -> cca.RRP
  paste0("R^2", cca.RRP[1], "~~italic(p)", cca.RRP[2]) -> rr.txt

  # unique contribution of each variable after accounting for other variables
  vegan::anova.cca(model.ord, by = 'margin') -> cca.env
  cca.env$`Pr(>F)` %>% .[-length(.)] %>% data.frame(color = .) %>% mutate(color = ifelse(color < .05, "#ef1828", "black")) -> arrow.p

  return(list('rr.txt' = rr.txt, 'arrow.p' = arrow.p))
}



#' @title Calculate ordination environment arrows
#' @description Projects environmental variables into an ordination.
#'
#' @param model.ord An ordination model returned by vegan.
#' @param env Environmental variables.
#' @return A data frame of projected environmental arrows.
#'
#' @family ordination
#' @keywords internal
#'

yy_envfit_arrow <- function(model.ord, env){
  vegan::envfit(model.ord, env) -> envfit
  envfit$vectors$pvals %>% data.frame(color = .) %>% mutate(color = ifelse(color < .05, "#ef1828", "black")) -> arrow.p
  (envfit$vectors$arrows[, 1:2] * envfit$vectors$r*2) %>% data.frame() %>% rownames_to_column('index') %>% cbind(arrow.p) -> arrow.df
  return(arrow.df)
}



#' @title Calculate individual environmental contributions
#' @description Calculates individual and grouped environmental contributions.
#' Calculates individual environmental contributions and, optionally, the
#' contribution of predefined variable groups.
#'
#' @param spe A community matrix with samples in rows.
#' @param env Environmental variables.
#' @param env.list Optional list defining groups of environmental variables.
#' @param method Ordination or grouping method.
#' @param type Explained-variation metric.
#' @param cols Optional color values.
#' @return A list containing individual and grouped contributions.
#'
#' @examples
#' \dontrun{
#'  y_rdacca_hp(spe, env, env.list = list(A = 1:3, B = 4:6))
#' }
#'
#' @family ordination
#' @export
#'

y_rdacca_hp <- function(spe, env, env.list = NULL, method =  c('RDA', 'dbRDA', 'CCA'), type = c('R2', 'adjR2'), cols = NULL){
  # env.list = list(A = 1:3, B = 4:7, C = 8:9) group of env metrics

  yy_check_pkg("rdacca.hp")

  if(is.null(cols)){RColorBrewer::brewer.pal(9, "Set3") -> cols}; data.frame(names(env), cols) -> e.cols
  if(!is.null(env.list)){cols[1:length(env.list)] %>% rep(times = purrr::map_int(env.list, length)) -> cols.grp; data.frame(names(env), cols.grp) -> e.cols}

  rdacca.hp::rdacca.hp(spe, env, method = method[1], type = type[1], scale = FALSE) -> VPA.ind; VPA.ind
  VPA.ind$Hier.part[, "I.perc(%)"] %>% sort(decreasing = TRUE) -> p.ind
  barplot(p.ind, col = e.cols[match(names(p.ind), e.cols[,1]), 2], las = 2)

  if(is.null(env.list)){ NULL -> VPA.gro} else{
    env.list %>% map(~{env[.x]}) -> env.grp; rdacca.hp::rdacca.hp(spe, env.grp, method = method, type = type, scale = F) -> VPA.gro
    VPA.gro$Hier.part[, "I.perc(%)"] %>% cbind() -> p.gro
    barplot(p.gro, col = unique(e.cols[[2]], sort = FALSE), ylab = "Proportion of explained variation (%)",  ylim = c(0,100))}

  return(list(VPA.ind = VPA.ind, VPA.grp = VPA.gro))
}



#' @title Plot points and environmental arrows
#' @description Draws ordination points, hulls, arrows and fit labels.
#'
#' @param point.df Point coordinates with columns for the two axes and
#' treatment.
#' @param arrow.df Arrow coordinates and colors.
#' @param eig Optional ordination eigenvalues.
#' @param rr.txt Optional formatted fit text.
#' @param n.mark Hull or ellipse marker configuration.
#' @param fig.name Optional data frame: input labels in column 1 and display
#' labels in column 2.
#' @param arrow.length Multiplier for environmental-arrow length.
#' @param xylim Symmetric ordination-axis limits.
#' @param show.point Whether to draw sample points instead of group summaries.
#' @param cols Optional color values.
#' @param shape Point shapes.
#' @param lab.xy Ordination-axis prefix.
#' @param legend.title Legend title.
#' @param output Output mode; see the function description.
#' @param show.density Whether to add marginal density plots.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object or a list containing plotting data.
#'
#' @family ordination
#' @keywords internal
#'

yy_ord_plot <- function(point.df, arrow.df = NULL, eig = NULL, rr.txt = NULL, n.mark = 0, fig.name = NULL, arrow.length = 1, xylim = 1,
                       show.point = TRUE, cols = NULL, shape = NULL, lab.xy = "Dim", legend.title = 'Group', output = 'p', show.density = FALSE, ...){
  # point.df should include Dim1, Dim2, trt.
  # arrow.df should include index, xend, yend, color.
  # eig is the feature root; Dim is the feature value.
  # n.mark = 3 is triangles; The larger the n.mark, the rounder the shadow.

  if(isTRUE(show.density)){yy_check_pkg(c("ggExtra", "ggplotify"),'install.packages(c("ggExtra", "ggplotify"))')}

  if(is.null(cols)){RColorBrewer::brewer.pal(8, "Set2") -> cols}; if(is.null(shape)){rep(15, 9) -> shape}
  point.df %>% rename_with(~c("Dim1", "Dim2", "trt")) %>% mutate(trt = factor(trt, levels = unique(trt))) -> point.df

  if(!show.point){ point.df %>% group_by(trt) %>%
    reframe(n = n(), se.1 = stats::sd(Dim1)/sqrt(n), se.2 = stats::sd(Dim2)/sqrt(n), Dim1 = mean(Dim1), Dim2 = mean(Dim2)) %>%
    ggplot(aes(x = Dim1, y = Dim2)) + geom_errorbar(aes(ymin = Dim2 - se.2, ymax = Dim2 + se.2, color = trt), width = 0) +
    geom_errorbarh(aes(xmin = Dim1 - se.1, xmax = Dim1 + se.1, color = trt), width = 0) + geom_point(aes(shape = trt, color = trt), size = 2) -> p
    } else{
    point.df %>% ggplot(aes(x = Dim1, y = Dim2)) + geom_point(aes(shape = trt, color = trt), size = 1.5) -> p
    if(n.mark == 0){p + ggforce::geom_mark_hull(aes(fill = trt, color = trt), expand = unit(2, "mm"), alpha = .2, size = .2) -> p
    } else{p + ggforce::geom_mark_ellipse(aes(fill = trt, color = trt), n = n.mark, expand = .05, tol = .01, alpha = .2, linetype = 0) -> p} }

  p + geom_hline(yintercept = 0, linetype = 2, size = .2) + geom_vline(xintercept = 0 ,linetype = 2, size = .2) +
    scale_colour_manual(legend.title, values = cols) + scale_shape_manual(legend.title, values = shape) + scale_fill_manual(legend.title, values = cols) +
    theme_classic() + theme(legend.background = element_blank(), legend.key.size = unit(.16, 'in')) + ylim(-xylim, xylim) + xlim(-xylim, xylim) -> p

  if(!is.null(eig)){paste0(' (', format(100*c(eig[1],eig[2])/sum(eig), digits = 3), '%)') -> eig};
  p + labs(x = paste0(lab.xy, ' 1', eig[1]), y = paste0(lab.xy, ' 2', eig[2])) -> p

  if(!is.null(arrow.df)){
    arrow.df %>% rename_with(~c("index", "xend", "yend", 'color')) %>% mutate(xend = xend * arrow.length, yend = yend * arrow.length) -> arrow.df
    if(!is.null(fig.name)){arrow.df %>% mutate(index = fig.name[match(index, fig.name[,1]), 2]) -> arrow.df}

    p + geom_segment(data = arrow.df, aes( x = 0, y = 0, xend = xend, yend = yend), size = .2, colour = arrow.df$color, arrow = arrow(length = unit(.1, 'cm'))) +
        geom_text(data = arrow.df, aes(x = xend, y = yend, label = index), size = 2.3, colour = arrow.df$color, hjust = (1 - sign(arrow.df$xend))/2, angle = (180/pi)*atan(arrow.df$yend/arrow.df$xend)) -> p
    }

  if(!is.null(rr.txt)){p + annotate("text", x = -xylim, y = xylim, label = rr.txt, parse = TRUE, size = 2.8, hjust = 'left', vjust = .38) -> p}
  if(show.density){p + theme(legend.position = c(.9,.8)) -> p; ggExtra::ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE) %>% ggplotify::as.ggplot() -> p}
  if(output != 'p'){return(list('point.df' = point.df, 'plot' = p))} else{return(p)}
}



#' @title Draw a principal-components ordination
#' @description Runs PCA and delegates plotting to yy_ord_plot.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param show.arrow Whether or how environmental arrows are displayed.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @examples
#' if (requireNamespace("vegan", quietly = TRUE)) {
#'   data(varespec, package = "vegan")
#'   trt <- factor(rep(c("A", "B"), each = 12))
#'   y_pca_plot(varespec[1:24, 1:20], trt)
#' }
#'
#' @family ordination
#' @export
#'

y_pca_plot <- function(spe, trt, show.arrow = TRUE, ...){

    if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}
    PCA <- vegan::rda(vegan::decostand(spe, "standardize"), scale = FALSE)
      PCA$CA$u[, 1:2] %>% cbind.data.frame(., trt) -> point.df
      PCA$CA$v[, 1:2] %>% data.frame() %>% rownames_to_column("index") %>% mutate(color = "black") -> arrow.df

    if(!show.arrow){NULL -> arrow.df}
    yy_ord_plot(point.df, arrow.df, PCA$CA$eig, lab.xy = 'PCA', ...) %>%
    return()
    }



#' @title Draw a principal-coordinates ordination
#' @description Runs PCoA and delegates plotting to yy_ord_plot.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param env Environmental variables.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'  y_pcoa_plot(spe, trt, env)
#' }
#'
#' @family ordination
#' @export
#'

y_pcoa_plot <- function(spe, trt, env = NULL, ...){

    if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}
    PCoA <- vegan::vegdist(spe, method = "bray") %>% stats::cmdscale(., k = 3, eig = TRUE)
      PCoA$points[, 1:2] %>% cbind.data.frame(., trt) -> point.df
      yy_adonis_summary(spe = spe, trt = trt, ...) -> adonis.txt

    if(!is.null(env)){yy_envfit_arrow(PCoA, vegan::decostand(env, "standardize")) -> arrow.df} else{NULL -> arrow.df}

    yy_ord_plot(point.df, arrow.df, PCoA$eig, rr.txt = adonis.txt$rr.txt, lab.xy = 'PCoA', ...) +
      geom_text(data = adonis.txt[[2]], aes(x = x, y = y, label = txt), size = 2.5, vjust = "bottom", hjust = "left") %>%
      return()
    }



#' @title Draw a non-metric multidimensional scaling ordination
#' @description Runs NMDS and delegates plotting to yy_ord_plot.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param env Environmental variables.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @family ordination
#' @export
#'

y_nmds_plot <- function(spe, trt, env = NULL, ...){

    if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}
    NMDS <- vegan::metaMDS(spe, distance = 'bray', k = 3, trymax = 100)
      NMDS$points[, 1:2] %>% cbind.data.frame(., trt) -> point.df
      yy_adonis_summary(spe = spe, trt = trt, ...) -> adonis.txt

    if(is.null(env) == FALSE){yy_envfit_arrow(NMDS, vegan::decostand(env, "standardize")) -> arrow.df} else{NULL -> arrow.df}

    yy_ord_plot(point.df, arrow.df, eig = NULL, rr.txt = adonis.txt$rr.txt, lab.xy = 'NMDS', ...) +
      geom_text(data = adonis.txt[[2]], aes(x = x, y = y, label = txt), size = 2.5, vjust = "bottom", hjust = "left") +
      annotate('text', x = Inf, y = -Inf, label = paste0('Stress = ', round(NMDS$stress, 2)), size = 3, vjust = -1, hjust = "right") %>%
    return()
    }



#' @title Draw a redundancy-analysis ordination
#' @description Runs RDA and delegates plotting to yy_ord_plot.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param env Environmental variables.
#' @param env.condition Optional conditioning variables for partial ordination.
#' @param show.arrow Whether or how environmental arrows are displayed.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#'  y_rda_plot(spe, trt, env)
#' }
#'
#' @family ordination
#' @export
#'

y_rda_plot <- function(spe, trt, env, env.condition = NULL, show.arrow = c('biplot', 'envfit'), ...){
    # env are variables of constraining/focus; env.condition = env[3:5, ] is the condition/covariates(partial-RDA)
    # envfit is ok, but biplot are precise (more strict in p-value). where env.condition should not appear as arrows.
    if(vegan::decorana(spe)$rproj[1] > 4){warning('0_o| DCA1 > 4 you should use CCA instead of RDA')}

    if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}
    vegan::decostand(spe, "hellinger") -> spe; vegan::decostand(env, "standardize") -> env

    if(is.null(env.condition)){vegan::rda(spe ~., data = env, scale = FALSE) -> RDA
    } else{
      cbind(env, env.condition) -> env.all; names(env.condition) -> cov.names
      paste("spe ~ . -", paste0(cov.names, collapse = " - "), "+ Condition(", paste(cov.names, collapse = " + "), ")") -> partial.model
      vegan::rda(as.formula(partial.model), data = env.all, scale = FALSE) -> RDA}

      RDA$CCA$u[, 1:2] %>% cbind.data.frame(., trt) -> point.df
      RDA$CCA$biplot[, 1:2] %>% data.frame() %>% rownames_to_column('index') -> arrow.df
      yy_anova_cca_summary(RDA) -> anova.cca.res

    if(show.arrow[1] == 'envfit'){yy_envfit_arrow(RDA, env) -> arrow.df} else {anova.cca.res$arrow.p %>% cbind(arrow.df, .) -> arrow.df}

    yy_ord_plot(point.df, arrow.df, RDA$CCA$eig, rr.txt = anova.cca.res$rr.txt, lab.xy = 'RDA', ...) %>%
    return()
    }



#' @title Draw a canonical-correspondence ordination
#' @description Runs CCA and delegates plotting to yy_ord_plot.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param env Environmental variables.
#' @param env.condition Optional conditioning variables for partial ordination.
#' @param show.arrow Whether or how environmental arrows are displayed.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @family ordination
#' @export
#'

y_cca_plot <- function(spe, trt, env, env.condition = NULL, show.arrow = c('biplot', 'envfit'), ...){

    if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}
    vegan::decostand(env, "standardize") -> env; # here, raw spe is Must!

    if(is.null(env.condition)){vegan::cca(spe ~., data = env) -> CCA
    } else{
      cbind(env, env.condition) -> env.all; names(env.condition) -> cov.names
      paste("spe ~ . -", paste0(cov.names, collapse = " - "), "+ Condition(", paste(cov.names, collapse = " + "), ")") -> partial.model
      vegan::cca(as.formula(partial.model), data = env.all) -> CCA}

      CCA$CCA$u[, 1:2] %>% cbind.data.frame(., trt) -> point.df
      CCA$CCA$biplot[, 1:2] %>% data.frame() %>% rownames_to_column('index') -> arrow.df
      yy_anova_cca_summary(CCA) -> anova.cca.res

    if(show.arrow[1] == 'envfit'){yy_envfit_arrow(CCA, env) -> arrow.df} else {anova.cca.res$arrow.p %>% cbind(arrow.df, .) -> arrow.df}

    yy_ord_plot(point.df, arrow.df, CCA$CCA$eig, rr.txt = anova.cca.res$rr.txt, lab.xy = 'CCA', ...) %>%
    return()
    }



#' @title Draw a distance-based redundancy ordination
#' @description Runs db-RDA and delegates plotting to yy_ord_plot.
#'
#' @param spe A community matrix with samples in rows.
#' @param trt Treatment grouping vector.
#' @param env Environmental variables.
#' @param env.condition Optional conditioning variables for partial ordination.
#' @param show.arrow Whether or how environmental arrows are displayed.
#' @param ... Additional arguments passed to underlying methods.
#' @return A ggplot object.
#'
#' @family ordination
#' @export
#'

y_db_rda_plot <- function(spe, trt, env, env.condition = NULL, show.arrow = c('biplot', 'envfit'), ...){

    if(missing(trt)){str_remove_all(row.names(spe), "\\d+") -> trt}
    vegan::vegdist(spe, "bray") -> dist; vegan::decostand(env, "standardize") -> env

    if(is.null(env.condition)){vegan::dbrda(dist ~ ., data = env, scale = FALSE) -> dbRDA
    } else{
      cbind(env, env.condition) -> env.all; names(env.condition) -> cov.names
      paste("spe ~ . -", paste0(cov.names, collapse = " - "), "+ Condition(", paste(cov.names, collapse = " + "), ")") -> partial.model
      vegan::dbrda(as.formula(partial.model), data = env.all, scale = FALSE) -> dbRDA}

    dbRDA$CCA$u[, 1:2] %>% cbind.data.frame(., trt) -> point.df
    dbRDA$CCA$biplot[, 1:2] %>% data.frame() %>% rownames_to_column('index') -> arrow.df
    yy_anova_cca_summary(dbRDA) -> anova.cca.res

    if(show.arrow[1] == 'envfit'){yy_envfit_arrow(dbRDA, env) -> arrow.df} else {anova.cca.res$arrow.p %>% cbind(arrow.df, .) -> arrow.df}

    yy_ord_plot(point.df, arrow.df, dbRDA$CCA$eig, rr.txt = anova.cca.res$rr.txt, lab.xy = 'db-RDA', ...) %>%
      return()
  }
