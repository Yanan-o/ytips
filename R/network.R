
# 0_o...test
# ---------------------------------------------------------------------------------------------------------------------------- #
# library(phyloseq)
# read.csv(onedrive('netpp/com_plant_2022_w.csv')) %>% mutate(ID = 'Plant') -> otu.df
# read.csv(onedrive('netpp/otu_fungi_2022_w.csv')) %>% mutate(ID = 'Fungi')-> otu.df.f
# data.frame(trt = rep(c('CK', 'K', 'P', 'R', 'G', 'F', 'L'), each = 7), block = 1:7) %>%
#   mutate(sample = paste0(trt, block), .before = trt) -> sample.df
#
# y_ps(otu.df, sample.df, rarefy = F) -> ps; ps
# ps %>% y_ps_bar(method = 'rich');
# ps %>% y_ps_bar(group = 'FG', method = 'r.abun')
# ps %>% subset_samples(trt == 'CK') %>% y_ps_network(v.group = 'FG', v.rank = T)
#
# y_ps(otu.df.f, sample.df) %>% y_ps_filter(filter.r.abun = 0.01) -> ps.f; ps.f
# merge_phyloseq(ps, ps.f) %>% subset_samples(trt == 'CK') -> ps.dd; ps.dd
#
# ps.dd %>% y_ps_igraph() -> igraph; igraph::vertex_attr(igraph)
# ps.dd %>% y_ps_network(p.adjust = 'fdr', c.group = 'ID', v.group = 'FG', show.module = 4)
# ps.dd %>% y_ps_network(spiec.easi = T, c.group = 'ID', v.group = 'FG', show.module = 3)
# ps.dd %>% y_ps_network(taxa.shuffle = T, p.adjust = 'none', c.group = 'ID', v.group = 'FG')
# ---------------------------------------------------------------------------------------------------------------------------- #



#' @title Build a phyloseq object from tabular data
#' @description Creates a phyloseq object from abundance, sample and taxonomy
#' tables.
#'
#' @param otu.df Abundance table whose first column contains taxon identifiers.
#' @param sample.df Sample-data table whose first column contains sample
#' identifiers.
#' @param taxa.df Optional taxonomy table.
#' @param rarefy Whether to rarefy samples to an even depth.
#' @return A phyloseq object.
#'
#' @examples
#' \dontrun{
#'  ps <- y_ps(otu.df, sample.df, taxa.df)
#' }
#'
#' @family network
#' @export
#'

y_ps <- function(otu.df, sample.df, taxa.df = NULL, rarefy = FALSE){
  # first column of sample.df are sample.id; otu.df & FG.df are species; all.df must have no row.names

  yy_check_phyloseq()

  if(!is.null(taxa.df)){merge(otu.df, taxa.df, by.x = names(otu.df)[1], by.y = names(taxa.df)[1], all.x = TRUE, sort = FALSE) -> otu.df}
  list(sample.df, otu.df) %>% map(~as_tibble(.x) %>% column_to_rownames(., var = names(.x)[1])) -> df_list

  df_list[[1]] %>% row.names() -> sample.ID;  df_list[[2]] %>% select(all_of(sample.ID)) -> otu;
  df_list[[2]] %>% select(!all_of(sample.ID)) %>% mutate(Species = row.names(.), .before = 1) -> taxa

  phyloseq::phyloseq(phyloseq::sample_data(df_list[[1]]), phyloseq::otu_table(as.matrix(otu), taxa_are_rows = TRUE), phyloseq::tax_table(as.matrix(taxa))) -> ps

  if(rarefy){phyloseq::rarefy_even_depth(ps) -> ps} # 稀释到相同测序深度
  return(ps)
  }



#' @title Filter and normalise a phyloseq object
#' @description Removes empty or constant taxa and applies abundance filters.
#' Removes empty samples and taxa before applying the abundance and prevalence
#' filters. prune = FALSE retains constant taxa. prevalence filters. prune =
#' FALSE retains constant taxa.
#'
#' @param ps A phyloseq object.
#' @param filter.n Keep the top N OTUs by total abundance when nonzero.
#' @param filter.r.abun Minimum mean relative abundance.
#' @param filter.sample Minimum sample prevalence and minimum count.
#' @param prune Whether to remove zero-variance taxa.
#' @return A filtered phyloseq object.
#'
#' @family network
#' @export
#'

y_ps_filter <- function(ps, filter.n = 0, filter.r.abun = 0, filter.sample = c(0, 3), prune = TRUE){
  # filter.n != 0, # 提取所有样本加总后，丰度前N个的otu
  # filter.r.abun = 1e-3, 提取所有样本平均后，相对丰度大于千份之一的OTU
  # filter.sample = (0.2, 3), 提取至少20％样本中，出现过3次及以上的OTU
  # prune = TRUE, ps.pruned 删除各样本中丰度均为0或相等的otu(值相等无法计算相关性)

  yy_check_phyloseq()

  ps %>% phyloseq::prune_samples(phyloseq::sample_sums(.)> 0, .) -> ps_0_sample
  message("0_0| Samples without any otu:~ ", paste0(setdiff(colnames(phyloseq::otu_table(ps)), colnames(phyloseq::otu_table(ps_0_sample))), collapse = ', '))

  ps_0_sample %>% phyloseq::prune_taxa(phyloseq::taxa_sums(.) > 0, .) -> ps_0
  message("0_0| otu = 0:~ ", dim(phyloseq::otu_table(ps))[1] - dim(phyloseq::otu_table(ps_0))[1], " otus were pruned.")

  if(prune){ps_0 %>% phyloseq::filter_taxa(., function(x) stats::sd(x) > 0, TRUE) -> ps.pruned} else{ps_0 -> ps.pruned}
  message("0_0| otu = equal value:~ ", dim(phyloseq::otu_table(ps_0))[1] - dim(phyloseq::otu_table(ps.pruned))[1], " otus were pruned.")

  # 按丰度筛选
  ps.pruned -> ps; phyloseq::transform_sample_counts(ps, function(x) x/sum(x)) -> ps.r # 转化OTU count数为相对丰度

  if(filter.n != 0){phyloseq::otu_table(ps) %>% as.data.frame() %>% arrange(desc(rowSums(.))) %>% head(filter.n) %>% row.names() -> otu.names.keep
    phyloseq::prune_taxa(otu.names.keep, ps) -> ps.filter} else{
      ps.r %>% phyloseq::filter_taxa(., function(x) mean(x) >= filter.r.abun, TRUE) -> ps.1
      ps %>% phyloseq::filter_taxa(., function(x) sum(x >= filter.sample[2]) >= (filter.sample[1] * length(x)), TRUE) -> ps.2
      intersect(phyloseq::otu_table(ps.1) %>% row.names(), phyloseq::otu_table(ps.2) %>% row.names()) -> otu.names.keep
    phyloseq::prune_taxa(otu.names.keep, ps) -> ps.filter}

  message("0_0| After filter:~ ", dim(phyloseq::otu_table(ps))[1] - dim(phyloseq::otu_table(ps.filter))[1], " otus were pruned. ",
          round(phyloseq::taxa_sums(ps.filter) %>% sum()/phyloseq::taxa_sums(ps) %>% sum(), 6) * 100, "% otu.abundance were kept(compare to ps.pruned): ")

  round(phyloseq::otu_table(ps.filter) %>% apply(., 2, function(x) sum(x))/phyloseq::otu_table(ps) %>% apply(., 2, function(x) sum(x)), 3) -> sum.abun_p
  if(!length(unique(sum.abun_p)) == 1){sum.abun_p %>% print()}  # 查看filter之后每个sample删除的物种丰度

 return(ps.filter)
 }



#' @title Aggregate phyloseq taxa into groups
#' @description Groups taxa and optionally collapses small groups into Others.
#'
#' @param ps A phyloseq object.
#' @param group Grouping variable, column name, or grouping table.
#' @param method Ordination or grouping method.
#' @param max.group Maximum number of retained groups.
#' @param rank Whether to rank taxa or groups.
#' @param show.abun Whether to include abundance in generated labels.
#' @param output Output mode; see the function description.
#' @return A data frame describing grouped taxa.
#'
#' @family network
#' @export
#'

y_ps_group <- function(ps, group = "", method = 'r.abun', max.group = 8, rank = FALSE, show.abun = TRUE, output = ""){
  # group = 'FG', 按FG进行分组，""时不分组
  # method = c('r.abun', 'abun', 'rich'), 计算relative abun./abundance/richness
  # max.group = 8, 按总丰度最多保留8个分组，其余的汇总至Others
  # show.abun = TRUE, 分组信息中展示分组水平的relative abun./abundance/richness (rowMeans)

  suppressMessages(y_ps_filter(ps)) -> ps; phyloseq::tax_table(ps) %>% as.data.frame() -> tax.0
  if(!group %in% names(tax.0)){message('0_0| group should be: ', paste0(names(tax.0), collapse = ', '));
     mutate(tax.0, group = 'na') -> tax.0}; tax.0 %>% select(group) %>% rename_with(~"group") -> tax

  if(method[1] == 'abun'){phyloseq::otu_table(ps) -> g.otu}
  if(method[1] == 'r.abun'){phyloseq::otu_table(ps) %>% apply(2, function(x) x/sum(x)) -> g.otu}
  if(method[1] == 'rich'){phyloseq::otu_table(ps) %>% apply(2, function(x) ifelse(x > 0, 1, 0)) -> g.otu}

  g.otu %>% as.data.frame() %>% mutate(row.mean = rowMeans(.), num.otu = ifelse(row.mean > 0, 1, 0)) %>%
    merge(., tax, by = "row.names", sort = FALSE) %>% group_by(group) %>% reframe(across(where(is.numeric), sum)) %>%
    ungroup() %>% mutate(rank = 1 + max(dense_rank(row.mean)) - dense_rank(row.mean)) -> g.otu

  g.otu %>% arrange(rank) %>% slice(1:max.group) -> g.otu.rank; g.otu %>% dplyr::filter(rank <= max(g.otu.rank$rank)) -> rank.2
  if(nrow(g.otu.rank) != nrow(rank.2)){warning('0_o| Share the same rank: ', paste(sort(rank.2$rank), collapse = ', '))}

  g.otu %>% dplyr::filter(!group %in% g.otu.rank$group) %>%
    reframe(across(where(is.numeric), sum)) %>% mutate(group = 'Others', rank = max.group) %>% rbind(g.otu.rank, .) -> g.otu.ok

  tax %>% rownames_to_column("otu") %>% mutate(group = ifelse(group %in% g.otu.ok[["group"]], group, "Others")) %>% merge(.,
    g.otu.ok, by = "group") %>% column_to_rownames("otu") %>% mutate(group = as.character(factor(group, levels = c(g.otu.ok[["group"]]))) ) -> tax.out

    if(method[1] != 'r.abun'){tax.out %>% mutate(group.r = paste0(group, " (", num.otu, ", ", formatC(row.mean, format = "f", digits = 1), ")")) -> tax.out}
    if(method[1] == 'r.abun'){tax.out %>% mutate(group.r = paste0(group, " (", num.otu, ", ", formatC(row.mean*100, format = "f", digits = 1), "%)")) -> tax.out}

  if(rank){tax.out %>% arrange(rank) -> tax.out}; if(show.abun){select(tax.out, group.r) -> A} else{select(tax.out, group) -> A};
  if(output == "otu"){A %>% rename_with(~group) %>% return()} else{tax.out %>% unique() %>% remove_rownames() %>% return()}
  }



#' @title Plot grouped taxon abundance
#' @description Plots grouped abundance with optional model-based significance
#' labels.
#'
#' @param ps A phyloseq object.
#' @param trt Treatment grouping vector.
#' @param group Grouping variable, column name, or grouping table.
#' @param method Ordination or grouping method.
#' @param cols Optional color values.
#' @param max.group Maximum number of retained groups.
#' @param rank.by Reference treatment used to rank groups.
#' @param rank Whether to rank taxa or groups.
#' @param pairwise Whether to perform all pairwise comparisons.
#' @param p.adjust P-value adjustment method.
#' @param fit.m Model expression used to test group differences. The default
#' uses a random intercept for a column named `block`; override it for other
#' random-effect or fixed-effect structures.
#' @return A ggplot object.
#'
#' @family network
#' @export
#'

y_ps_bar <- function(ps, trt = "", group = "", method = c('abun', 'r.abun', 'rich'), cols = NULL, max.group = 8, rank.by = 'CK',
                     rank = FALSE, pairwise = FALSE, p.adjust = 'default', fit.m = 'lmer(y ~ x + (1|block), data = .x)'){

  # rank.by = 'CK' 按CK排序提取最多8个组
  # fit: **模型检验y ~ x & data = .x 形式不变, 其它可写(变量必须存在于sample_data)
  if(isTRUE(pairwise)){yy_check_pkg("rcompanion")}

  if(is.null(cols)){RColorBrewer::brewer.pal(n = 9, "Set1") %>% scales::alpha(., .8) -> cols}
  if(p.adjust == 'default'){ifelse(pairwise, 'tukey', 'dunnett') -> p.adjust}

  if(group == '' && method[1] == 'r.abun'){stop('0_o| r.abun make no sense!')}
  switch(method[1], r.abun = 'Relative abundance', abun  = 'Abundance', rich = 'Richness') -> label.y
  ifelse(trt == '', names(phyloseq::sample_data(ps))[1], trt) -> x; phyloseq::sample_data(ps) %>% data.frame() %>% rename('x' = x) -> trt.df

  if(rank.by %in% trt.df$x){eval(parse(text = paste0("phyloseq::subset_samples(ps, ", x, "=='",  rank.by, "')"))) -> ps.rank} else{ps -> ps.rank}
  y_ps_group(ps.rank, group = group, method = method[1], max.group = max.group, rank = rank) %>% .$group %>% .[.!= "Others"] -> g.name

  y_ps_group(ps, group = group, method = method[1], max.group = 99, rank = FALSE) %>%
    mutate(group = ifelse(group %in% g.name, group, 'Others')) %>% group_by(group) %>%
    reframe(across(where(is.numeric), sum)) %>% pivot_longer(cols = !group, values_to = 'y') %>%
    merge(trt.df, by.x = 'name', by.y = 'row.names', sort = FALSE) -> gg.after.group

  gg.after.group %>% split(f = factor(.$group)) %>% imap(~{
     eval(parse(text = fit.m)) -> fit; emmeans::emmeans(fit, ~x) -> emm
     if(pairwise){
       emmeans::contrast(emm,  method = "pairwise", adjust = p.adjust) %>% as.data.frame() %>%
       rcompanion::cldList(p.value ~ contrast, reversed = FALSE, data = .) %>% mutate(p.sig = Letter, x = Group) -> est.p;
       if(length(unique(est.p$p.sig)) == 1){' ' -> est.p$p.sig} } else{
      emmeans::contrast(emm,  method = "trt.vs.ctrl", adjust = p.adjust) %>% as.data.frame() %>%
      mutate(p.sig = ifelse(p.value < .01, "**", ifelse(p.value < .05, '*', " ")), x = gsub(' -.*', '',contrast)) -> est.p }

    merge(data.frame(emm), est.p[c('x', 'p.sig')], by = 'x', all.x = T) }) %>% bind_rows(.id ='group') %>%
    mutate(x = gsub('x', '', x), p.sig = ifelse(is.na(p.sig), ' ',p.sig)) %>% select(group, x, p.sig) -> gg.p.sig

  gg.after.group %>% group_by(x, group) %>% reframe(y = mean(y)) %>% merge(gg.p.sig, by = c('group', 'x'), all.x = TRUE, sort = FALSE) %>%
    dplyr::filter(y != 0) %>% mutate(group = factor(group, levels = c(g.name, "Others")), x = factor(x, levels = unique(trt.df$x))) -> df.p

  ggplot(df.p, aes(x = x, y = y, fill = group)) + geom_bar(stat = "identity", width = .8) +
    geom_text(aes(label = p.sig), position = position_stack(vjust = ifelse(length(unique(df.p$group)) == 1, 1.1, .3)), stat = "identity", size = 3.3) +
    scale_fill_manual(values = cols) + theme_bw() + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1), axis.text = element_text(colour = 1)) +
    theme(legend.key.size = unit(.16, 'in') , plot.title = element_text(hjust = .5, face = "bold", size = 9)) + labs(fill = group, x = '', y = label.y) %>%
    return()
}



#' @title Build and display a co-occurrence network
#' @description Builds a network and returns graph, module and Zi-Pi results.
#' Builds the igraph object, calculates network properties, draws the network
#' and returns module and Zi-Pi information. network, and returns module and
#' Zi-Pi information.
#'
#' @param ps A phyloseq object.
#' @param r.threshold Absolute correlation threshold.
#' @param p.threshold Correlation p-value threshold.
#' @param r.method Correlation method.
#' @param p.adjust P-value adjustment method.
#' @param taxa.shuffle Whether to generate p-values by taxon shuffling.
#' @param spiec.easi Whether to use SpiecEasi instead of pairwise correlations.
#' @param c.group Optional cross-kingdom or higher-level grouping table.
#' @param v.group Optional node-grouping table.
#' @param v.method Abundance metric used for node groups.
#' @param v.max.group Maximum number of node groups.
#' @param v.rank Whether to rank node groups.
#' @param show.abun Whether to include abundance in generated labels.
#' @param v.size Optional node sizes.
#' @param v.cols Node-group colors.
#' @param m.cols Module colors.
#' @param show.module Minimum module size to display.
#' @param m.gamma Leiden resolution parameter.
#' @param isolate.node Whether to retain isolated nodes.
#' @param layout igraph layout expression.
#' @param layout.size Layout scaling values.
#' @return A list containing the network, node and edge data, properties, and
#' Zi-Pi results.
#'
#' @examples
#' \dontrun{
#'  network <- y_ps_network(ps, c.group = "Kingdom", v.group = "FG")
#' }
#'
#' @family network
#' @export
#'

# 输入ps文件，输出cor, igraph, otu.g, property，并展示共线网络图和主要模块的组成信息
# otu.g是由otu，taxa，module，degree等信息构成的df文件，适合后续分析与作图
y_ps_network <- function(ps, r.threshold = 0.7, p.threshold = 0.05, r.method = "spearman", p.adjust = "fdr", taxa.shuffle = FALSE, spiec.easi = FALSE,
                       c.group = NULL, v.group = NULL, v.method = 'r.abun', v.max.group = 8, v.rank = FALSE, show.abun = TRUE, v.size = NULL, v.cols = NULL,
                       m.cols = NULL, show.module = 8, m.gamma = 0, isolate.node = TRUE, layout = "igraph::layout_nicely()", layout.size = c(.6, .8, 1)){

  # r., p.: 控制相关系数的阈值、计算方法
  # v., m.: 定义节点和模块分组、属性、计算方法
  # spiec.easi = TRUE 自动稀疏模型；taxa.shuffle = TRUE 但采用分类洗牌计算p值
  # layout: 控制布局, see igraph, layout_in_circle()*.8, 可控制大小, 同layout.size
  # c.group: 控制多级网络布局, 跨界网络分组, e.g., plant, nematode, microbe...
  # v.group: 控制节点分组（与otu匹配的数据框或提供分类名由y_ps_group生成, 同c.group）

  yy_check_network_dependencies(spiec.easi, taxa.shuffle)

  if(is.character(c.group)){c.group -> c.gro; y_ps_group(ps, group = c.gro, max.group = 99, rank = FALSE, show.abun = FALSE, output = 'otu') -> c.group}

  if(is.character(v.group)){v.group -> v.gro;
    if(is.null(c.group)){y_ps_group(ps, group = v.gro, method = v.method, max.group = v.max.group, rank = v.rank, show.abun = show.abun, output = 'otu') -> v.group
    } else{
      data.frame() -> v.group
      for(i in unique(c.group)[,1]){
        eval(parse(text = paste0("phyloseq::subset_taxa(ps, ", c.gro, "=='",  i, "')"))) -> sub.ps
        y_ps_group(sub.ps, group = v.gro, method = v.method, max.group = v.max.group, rank = v.rank, show.abun = show.abun, output = 'otu') -> v.group.i
        rbind.data.frame(v.group, v.group.i) -> v.group}}
  }

  # igraph
    y_ps_igraph(ps, spiec.easi = spiec.easi, r.threshold = r.threshold, p.threshold = p.threshold,
                r.method = r.method, taxa.shuffle = taxa.shuffle, p.adjust = p.adjust, isolate.node = isolate.node) -> igraph

    yy_network_properties(igraph) -> net_property;

    yy_network_plot(igraph, c.group = c.group, v.group = v.group, v.size = v.size, v.cols = v.cols,
                 m.cols = m.cols, show.module = show.module, m.gamma = m.gamma, layout = layout, layout.size = layout.size) -> df_igraph

  # Zi-Pi (Zi > 2.5, 模块内核心物种; Pi > 0.62, 模块间桥梁物种; Zi > 2.5 & Pi > 0.62, 全局关键物种; Others, 边缘物种)
  df_igraph$V.df[, c('name', 'color')] -> cols.v
  if(is.null(c.group)){NULL -> pchs.c} else{df_igraph$V.df[, c('name', 'c.group')] -> pchs.c}
  yy_network_zipi(igraph, cols.df = cols.v, pchs.df = pchs.c) -> zipi; # zipi$p
  zipi$res.df %>% mutate(cols.0 = cols.v[match(Node, cols.v$name), 3]) -> zipi$res.df

  list("ps" = ps, "igraph" = igraph, 'V.df' = df_igraph$V.df, 'E.df' = df_igraph$E.df, "property" = net_property, 'zipi.res' = zipi$res.df) -> netpp
  return(invisible(netpp))
}



#' @title Build an igraph network from phyloseq data
#' @description Builds a correlation or SpiecEasi network as an igraph object.
#' Uses either pairwise correlations or SpiecEasi to build the network.
#' P-values can be adjusted or estimated by taxon shuffling.
#'
#' @param ps A phyloseq object.
#' @param r.threshold Absolute correlation threshold.
#' @param p.threshold Correlation p-value threshold.
#' @param r.method Correlation method.
#' @param p.adjust P-value adjustment method.
#' @param taxa.shuffle Whether to generate p-values by taxon shuffling.
#' @param spiec.easi Whether to use SpiecEasi instead of pairwise correlations.
#' @param isolate.node Whether to retain isolated nodes.
#' @return An igraph object.
#'
#' @family network
#' @export
#'

# 分类洗牌计算p值（taxa shuffle，Column-wise independent shuffle，将otu丰度打乱计算cor）
y_ps_igraph <- function(ps, r.threshold = 0.7, p.threshold = 0.05, r.method = 'spearman',
                        p.adjust = 'fdr', taxa.shuffle = FALSE, spiec.easi = FALSE, isolate.node = TRUE){

  # spiec.easi = TRUE 忽视r,p参数，自动选择稀疏模型
  # taxa.shuffle = TRUE 基于相关性cor，但采用分类洗牌计算p值
  # 基于相关性是转化成r.abun后计算网络是推荐的，但SpiecEasi则必须使用OTU数目
  # p.adjust = "fdr"，能较好地平衡网络稀疏性与假阳性率，被广泛接受（'none'不设置）
  # p & r.threshold，用来过滤边，假阳性已双重控制
    yy_check_network_dependencies(spiec.easi, taxa.shuffle)

    y_ps_filter(ps) -> ps

    # 基于SPIEC-EASI (采用稀疏逆协方差估计 + StARS稳定性选择,自动选择稀疏模型，无需控制r或p值)
    if(spiec.easi){ warning('0_o| Using Spiec.Easi to built "igraph", the r.& p. parameters are invalid')
      SpiecEasi::spiec.easi(ps, method = 'mb', lambda.min.ratio = 1e-2, nlambda = 20, pulsar.params = list(rep.num = 50)) -> Spiec.E
      SpiecEasi::getRefit(Spiec.E) -> adj; SpiecEasi::getOptBeta(Spiec.E) -> sig
      phyloseq::taxa_names(ps) -> rownames(adj); phyloseq::taxa_names(ps) -> colnames(adj);
      adj[sig < 0] <- -1; adj[sig > 0] <- 1; as.matrix(adj) -> adj  } else{

    # 基于相关性矩阵('pearson', 'spearman')
    warning('0_o| Using Cor to built "igraph", the r. & p. parameters are passed in')
    phyloseq::transform_sample_counts(ps, function(x) x/sum(x)) -> ps.r
    phyloseq::otu_table(ps.r) %>% as.data.frame() %>% t() -> otu_mat
    otu_mat %>% cor(., method = r.method) -> cor

    #分类洗牌计算p值（taxa shuffle，Column-wise independent shuffle，将otu丰度打乱计算cor）
    if(taxa.shuffle){
        n <- ncol(otu_mat); null_cor <- array(0, c(n, n, 999))  # 存放零分布cor

        for (k in 1:999) { # 每列（otu）内独立打乱行顺序，计算零分布cor
          apply(otu_mat, 2, sample) -> shuffled;
          cor(shuffled, method = r.method) -> null_cor[,,k] }

        p_val <- matrix(1, n, n) # 计算经验 p 值（双尾）
        for (i in 1:(n-1)) { for (j in (i+1):n) { # p既零分布中有多少比例r的绝对值 ≥ 观测绝对
            p_val[i,j] <- mean(abs(null_cor[i,j,]) >= abs(cor[i,j]))
            p_val[j,i] <- p_val[i,j]} }

        } else{ # 基于正态分布，通常测序数据不具备该条件
          otu_mat %>% Hmisc::rcorr(., type = r.method) -> cor.rp; cor.rp$r -> cor; cor.rp$P -> p_val
         }

    # 校正p，筛选r, p值，并构建igraph
    p_val[lower.tri(p_val)] %>% stats::p.adjust(., method = p.adjust) -> p_adj
    matrix(1, nrow(cor), ncol(cor)) -> p; p[lower.tri(p)] <- p_adj; p[upper.tri(p)] <- t(p)[upper.tri(p)]

    cor[p > p.threshold | abs(cor) < r.threshold] = 0; diag(cor) <- 0; as.matrix(cor) -> adj } # 对角线为0的cor
    igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE) -> igraph
    if(!isolate.node){igraph::delete.vertices(igraph, names(igraph::degree(igraph)[igraph::degree(igraph) == 0])) -> igraph}
  return(igraph)
  }



#' @title Draw a grouped co-occurrence network
#' @description Draws and annotates a module-aware network. Supports one or
#' more c.group clusters arranged around the graph. Node colors and sizes can
#' be supplied through v.group and v.size. Supports one or more c.group
#' clusters arranged around the graph. Node colors and sizes can be supplied
#' through v.group and v.size.
#'
#' @param igraph An igraph network object.
#' @param c.group Optional cross-kingdom or higher-level grouping table.
#' @param v.group Optional node-grouping table.
#' @param v.size Optional node sizes.
#' @param v.cols Node-group colors.
#' @param m.cols Module colors.
#' @param show.module Minimum module size to display.
#' @param m.gamma Leiden resolution parameter.
#' @param layout igraph layout expression.
#' @param layout.size Layout scaling values.
#' @return A list containing node and edge tables, invisibly.
#'
#' @family network
#' @keywords internal
#'

yy_network_plot <- function(igraph, c.group = NULL, v.group = NULL, v.size = NULL, v.cols = NULL,
                         m.cols = NULL, show.module = 8, m.gamma = 0,
                         layout = "igraph::layout_nicely()", layout.size = c(.6, .8, 1)){

    # r., p.: 控制相关系数的阈值、计算方法
    # v., m.: 定义节点和模块分组、属性、计算方法
    # show.module = 8: 只展示node个数 >= 8的模块
    # layout: 控制布局, see igraph, layout_in_circle()*.8, 可控制大小, 同layout.size
    # c.group: 控制多级网络布局, 跨界网络分组, e.g., plant, nematode, microbe...
    # v.group: 控制节点分组（与otu匹配的数据框或提供分类名由y_ps_group生成, 同c.group）

    if(is.null(m.cols)){RColorBrewer::brewer.pal(n = 8, "Set2") -> m.cols}
    if(is.null(v.cols)){RColorBrewer::brewer.pal(n = 9, "Set1") -> v.cols}

    # 调整V(igraph)和E(igraph)属性
    module.group <- yy_network_modularity(igraph, m.gamma = m.gamma)
    igraph::V(igraph)$module <- igraph::membership(module.group); igraph::V(igraph)$degree <- igraph::degree(igraph)

    which(igraph::sizes(module.group) >= show.module) -> show.module.id
    m.colbar <- rep("gray90", 999); m.colbar[show.module.id] <- m.cols[1:length(show.module.id)]
    module.group[show.module.id] -> mark.groups; scales::alpha(m.colbar[show.module.id], .3) -> mark.col

    c = igraph::V(igraph)$degree; if(diff(range(c)) == 0){rep(2, igraph::vcount(igraph)) -> igraph::V(igraph)$size} else{((c-min(c))/(max(c)-min(c)))*6 + 2 -> igraph::V(igraph)$size}
    igraph::V(igraph)$color <- m.colbar[igraph::membership(module.group)]
    igraph::E(igraph)$color <- ifelse(igraph::E(igraph)$weight > 0, "skyblue", "pink")

    if(is.null(v.group)){warning("0_o| v.color is drawed by module")} else{
      v.group %>% as.data.frame() %>% rename_with(~"v.group") -> v.group
      v.group[["v.group"]][match(igraph::V(igraph)$name, row.names(v.group))] -> igraph::V(igraph)$v.group
      data.frame(unique(v.group), color = v.cols[1:nrow(unique(v.group))]) -> v.colbar
      v.colbar$color[match(igraph::V(igraph)$v.group, v.colbar$v.group)]  -> igraph::V(igraph)$color}

    if(is.null(v.size)){warning("0_o| v.size is drawed by degree")} else{
      if(is.numeric(v.size)){igraph::V(igraph)$size <- v.size} else {
        v.size %>% as.data.frame() %>% rename_with(~"c") %>%
          mutate(size = ((c-min(c))/(max(c)-min(c)))*6 + 2) -> v.size
        v.size[["size"]][match(igraph::V(igraph)$name, row.names(v.size))] -> igraph::V(igraph)$size}}


    if(is.null(c.group)){
      eval(parse(text = paste0('igraph %>%', layout, '*', layout.size[1],' -> layout.c')));
      } else{

      # 为每个分组生成独立圆形布局
      c.group %>% as.data.frame() %>% rename_with(~"c.group") -> c.group
      c.group[["c.group"]][match(igraph::V(igraph)$name, row.names(c.group))] -> igraph::V(igraph)$c.group

      row.names(c.group)

        igraph::vertex_attr(igraph, 'c.group') -> groups; unique(groups) %>% length() -> n
        layout.size <- rep(layout.size, length.out = n)
        matrix(0, nrow = igraph::vcount(igraph), ncol = 2) -> layout.c

        for(i in 1:n){
          igraph::induced_subgraph(igraph, groups == unique(groups)[[i]]) -> sub.igraph
          eval(parse(text = paste0('sub.igraph %>% ', layout, '*', layout.size[i],' -> sub_layout')))

          if(n == 1){x_center <- 0; y_center <- 0}

          if(n == 2){
            x_center <- (i - 1) * 8  # 水平间隔
            y_center <- (i - 1) * -8} # 垂直间隔

          if(n == 3){
            angles <- c(90, 210, 330) * pi / 180  # 三角形分布（弧度）
            radius <- 8
            x_center <- radius * cos(angles[i]) # 计算顶点坐标
            y_center <- radius * sin(angles[i])} # 计算顶点坐标

          if(n > 3){
            angles <- 2 * pi * c(0:(n - 1))/n
            radius <- 4 + 2 * n
            x_center <- radius * cos(angles[i]); y_center <- radius * sin(angles[i])}

          layout.c[groups == unique(groups)[[i]], ] <- t(t(sub_layout) + c(x_center, y_center))} # 偏移坐标

       # 为每个分组构建模块图（网络分组中同色为相同模块）
        list() -> mark.groups
        c.group %>% split(f = c.group) %>% imap(~{row.names(.x)}) -> c.group.list

        for(j in 1:n){
          module.group[show.module.id] %>% imap(~{.x[.x %in% c.group.list[[j]]]}) -> sub.mark.groups
          c(mark.groups, sub.mark.groups) -> mark.groups}
       }

      # plot
      if(is.null(igraph::E(igraph)$weight)){1 -> edge.width}else{abs(igraph::E(igraph)$weight) -> edge.width}
      igraph::plot.igraph(igraph, layout = layout.c,
        vertex.label = NA, vertex.shape = "circle", edge.curved = FALSE,
        vertex.size = igraph::V(igraph)$size, vertex.frame.color = "black", edge.width = edge.width,
        mark.groups = mark.groups, mark.col = mark.col, mark.border = NA)


      igraph::vertex_attr(igraph) %>% as.data.frame() -> V.df; igraph::as_data_frame(igraph) -> E.df
      if(nrow(E.df) == 0){mutate(E.df, weight = numeric()) -> E.df}
      E.df %>% mutate(edge.w = ifelse(weight > 0, ' + ', ' - ')) -> E.df; E.df %>% count(edge.w, name = "n", sort = F) -> edge.w

      # 计数节点
       if(!is.null(c.group)){V.df %>% count(c.group, name = "n.node", sort = F) %>% print(); cat('\n')} else{NA -> c.group}
       if(!is.null(v.group)){V.df %>% count(v.group, name = "n.node", sort = F) %>% print(); cat('\n')} else{NA -> v.group}

      # 计数边
        merge(v.group, c.group, by = 'row.names', all = T, sort = F) %>% rename_with(~c('name', 'v.group', 'c.group')) -> vc.group
        E.df %>% merge(., vc.group, by.x = 'to', by.y = 'name', all.x = T, sort = F) %>%
          merge(., vc.group, by.x = 'from', by.y = 'name', all.x = T, sort = F, suffixes = c('.to','.from')) %>%
          mutate(edge.v = (paste0(edge.w, v.group.from, "<-->", v.group.to)), edge.c = (paste0(edge.w, c.group.from, "<-->", c.group.to))) %>%
          mutate(edge.c.in = paste0(edge.w, ifelse(c.group.from == c.group.to, 'Intra', 'Inter'))) %>%  select(where(~ !all(is.na(.)) | length(.) == 0)) -> E.df

        if(!is.null(nrow(c.group))){print(count(E.df, edge.c.in, name = "n", sort = F)); print(count(E.df, edge.c, name = "n", sort = F)); cat('\n')}
        if(!is.null(nrow(v.group))){count(E.df, edge.v, name = "n", sort = F) %>% print(); cat('\n')}

      if(!is.null(nrow(v.group))){
        legend("topright", title = ifelse(grepl("%", v.colbar[['v.group']][1]), "Group (N, r.abun)", "Group"),
               legend = v.colbar[['v.group']], pt.cex = 1.5, pch = 16, col = v.colbar[['color']], text.col = 1, border = 1, bty = "n", xjust = 0, title.adj = .2)}

        legend("right", inset = c(0, -.2), title = "Correlation (N)", col = c("skyblue", "pink"),
               legend = paste0(edge.w$edge.w, "(", edge.w$n, ")"), lty = 1,  bty = "n", xjust = 0, title.adj = .5)

        legend("bottomleft", title = "Module (N)", fill = scales::alpha(m.colbar[show.module.id], .8), border = 1,
           legend = paste0("#", show.module.id, " (", igraph::sizes(module.group)[show.module.id], ")"),  bty = "n", xjust = 0, title.adj = .5)

      list(V.df = V.df, E.df = E.df) -> result; return(invisible(result))
    }



#' @title Detect network modules
#' @description Runs Leiden community detection.
#'
#' @param igraph An igraph network object.
#' @param m.gamma Leiden resolution parameter.
#' @return An igraph community-detection result.
#'
#' @family network
#' @keywords internal
#'

# Leiden + CPM 支持正负边的权重计算，算法会将正边集中到社区内部，将负边放到社区之间
# 理论上：划分应使模块内部以正边为主，模块之间以负边（或弱边）为主
# m.gamma = 0 表示 weight < 0的边视为排斥，> 0的边视为吸引，与spearman相关性符合
yy_network_modularity <- function(igraph, m.gamma = 0){
    igraph::cluster_leiden(igraph, objective_function = "CPM", weights = igraph::E(igraph)$weight, resolution_parameter = m.gamma, n_iterations = 99)
    }



#' @title Calculate Zi-Pi network roles
#' @description Calculates within- and among-module connectivity.
#'
#' @param igraph An igraph network object.
#' @param cols.df Optional node color table.
#' @param pchs.df Optional node shape table.
#' @param m.gamma Leiden resolution parameter.
#' @return A list containing a Zi-Pi plot and role table.
#'
#' @family network
#' @keywords internal
#'

yy_network_zipi <- function(igraph, cols.df = NULL, pchs.df = NULL, m.gamma = 0){
  # Zi > 2.5, 模块内核心物种;Pi > 0.62, 模块间桥梁物种; Zi > 2.5 & Pi > 0.62, 全局关键物种; Others, 边缘物种

    yy_network_modularity(igraph, m.gamma = m.gamma) %>% igraph::membership() -> modules

    calculate_zi <- function(igraph, modules){
      igraph::as_adjacency_matrix(igraph, attr = "weight", sparse = FALSE) -> adj
      sapply(1:igraph::vcount(igraph), function(i){same_module <- which(modules == modules[i]); sum(adj[i, same_module])}) -> k_in
      tapply(k_in, modules, function(k_module){s <- stats::sd(k_module)
      if (is.na(s)||s == 0) return(rep(0, length(k_module))); (k_module - mean(k_module))/s }) -> zi; unlist(zi)[order(modules)]}

    calculate_pi <- function(igraph, modules){
      unique(modules) -> am; setNames(seq_along(am), am)[as.character(modules)] -> all_modules
      sapply(1:igraph::vcount(igraph), function(i){ if(igraph::degree(igraph)[i] == 0) return(0)
        tabulate(all_modules[igraph::neighbors(igraph, i)], nbins = length(all_modules)) -> mod_connections
        mod_connections / igraph::degree(igraph)[i] -> proportions; 1 - sum(proportions^2)}) -> pi_values; pi_values}

    zi_vals <- calculate_zi(igraph, modules)
    pi_vals <- calculate_pi(igraph, modules)

    data.frame(Node = igraph::V(igraph)$name, zi = zi_vals, pi = pi_vals) %>%
      mutate(roles = ifelse(zi_vals > 2.5 & pi_vals > 0.62, "Network Hub",
                            ifelse(zi_vals > 2.5, "Module Hub", ifelse(pi_vals > 0.62, "Connector", "Peripheral")))) -> res.df.0

   res.df.0 -> df.p
   if(is.null(cols.df)){df.p$roles -> df.p$cols} else{df.p %>% mutate(cols = cols.df[match(Node, cols.df[,1]), 2]) -> df.p}
   if(is.null(pchs.df)){16 -> df.p$pchs} else{df.p %>% mutate(pchs.0 = pchs.df[match(Node, pchs.df[,1]), 2],) %>%
       mutate(pchs = c(15, 17, 16, 18)[match(pchs.0, unique(pchs.0))])-> df.p}

    message('0_0| Network Hub: OTUs in both ~ zi & pi ~ ........ \n')
    message('Module Hub(zi > 2.5): '); df.p[df.p$zi > 2.5, 'Node'] %>% print()
    message('Connector(pi > 0.62): '); df.p[df.p$pi > 0.62, 'Node'] %>% print()

    ggplot(df.p, aes(x = pi, y = zi)) +
     geom_hline(yintercept = 2.5, lty = 2, linewidth = .3) + geom_vline(xintercept = .62 ,lty = 2, linewidth = .3) +
     geom_point(color = df.p$cols, shape = df.p$pchs) + theme_bw() + theme(panel.grid = element_blank()) +
     labs(x = 'Among-module connectivity (Pi)', y = 'Within-module connectivity (Zi)') -> p
   return(invisible(list(plot = p, res.df = df.p)))
  }



#' @title Calculate network properties
#' @description Returns common global network properties.
#'
#' @param igraph An igraph network object.
#' @param m.gamma Leiden resolution parameter.
#' @param n.hub Whether to return hub information.
#' @return A data frame of global network properties.
#'
#' @family network
#' @keywords internal
#'

yy_network_properties <- function(igraph, m.gamma = 0, n.hub = FALSE){

    num.edges <- length(igraph::E(igraph)) #network connectivity (L): 边数
    num.vertices <- length(igraph::V(igraph)) #network size (n): 节点数

    #connectance: 连接度，网络中实际相互作用（连接数之和）占潜在相互作用（潜在连接数）的比例，反应网络复杂程度
    connectance <- igraph::edge_density(igraph, loops = FALSE)
    average.degree <- mean(igraph::degree(igraph)) #average degree: 平均度(2*num.edges/num.vertices)

   if(!is.null(igraph::E(igraph)$weight)){igraph.weight <- igraph::E(igraph)$weight; igraph::E(igraph)$weight = abs(igraph::E(igraph)$weight)}
    average.path.length <- igraph::average.path.length(igraph, directed = FALSE, unconnected = TRUE, weights = NULL) #average path length: 平均路径长度
    diameter <- igraph::diameter(igraph, directed = FALSE, unconnected = TRUE, weights = NULL) #diameter: 直径

   if(!is.null(igraph::E(igraph)$weight)){igraph::E(igraph)$weight = igraph.weight}
    edge.connectivity <- igraph::edge_connectivity(igraph)  #edge connectivity or group adhesion: 群连通度

    #average clustering coefficient (average CC): 聚集系数
    #反映网络中节点紧密关系，也称传递性。分为整体和局域聚集系数。整体聚集系数C表征整个网络的平均的“成簇性质”
    clustering.coefficient <- igraph::transitivity(igraph, type = "average")
    no.clusters <- igraph::no.clusters(igraph) #  no.clusters

    centralization.degree <- igraph::centralization.degree(igraph)$centralization #degree centralization: 度中心
    centralization.betweenness <- igraph::centralization.betweenness(igraph)$centralization #between centralization: 介数中心性
    centralization.closeness <- igraph::centralization.closeness(igraph)$centralization

    if(!is.null(igraph::E(igraph)$weight)){
      num.pos.edges <- sum(igraph.weight > 0); num.neg.edges <- sum(igraph.weight < 0)}else{num.pos.edges <- 0; num.neg.edges <- 0}


    igraph.network.pro <- rbind(num.edges, num.pos.edges, num.neg.edges,
                                num.vertices, connectance, average.degree, average.path.length,
                                diameter, edge.connectivity, clustering.coefficient,
                                no.clusters, centralization.degree, centralization.betweenness,
                                centralization.closeness)

    rownames(igraph.network.pro) <- c("num.edges(L)", "num.pos.edges",
                                      "num.neg.edges", "num.vertices(n)", "Connectance(edge_density)",
                                      "average.igraph::degree(Average K)", "average.path.length", "diameter",
                                      "edge.connectivity", "mean.clustering.coefficient(Average.CC)",
                                      "no.clusters", "centralization.degree", "centralization.betweenness",
                                      "centralization.closeness")

    return(data.frame(value = round(igraph.network.pro, 4)))
  }



#' @title Plot module composition
#' @description Plots module composition with radial stacked bars.
#'
#' @param otu.g Network node table returned by y_ps_network() .
#' @param group Grouping variable, column name, or grouping table.
#' @param output Output mode; see the function description.
#' @param v.cols Node-group colors.
#' @param module.n Maximum number of modules to display.
#' @param module.id Optional module identifiers.
#' @return A ggplot object.
#'
#' @family network
#' @keywords internal
#'

yy_network_module_plot <- function(otu.g, group = "v.group", output = "r.abun",
                                   v.cols = NULL, module.n = 7, module.id = NULL){

  # module.n = 7 展示7个最大的module
  # module.id = c(2, 3) 展示模块2和3
  # output = "r.abun" 展示每个模块中不同分类otu的相对丰度
  # output = "num.otu" 展示每个模块中不同分类otu的多度
  # output = "r.abun" 展示每个模块中不同分类otu的重要性（num.otu * degree）

  # otu.g: the output from y_ps_network (tiff, 6:4)
  if(is.null(v.cols)){RColorBrewer::brewer.pal(n = 9, "Set1") -> v.cols}

  otu.g %>% mutate("group" = factor(.[, group], levels = unique(.[, group]))) %>% group_by(module, group) %>%
    reframe(r.abun = sum(r.abun), num.otu = length(group), mean.degree = mean(degree)) %>% ungroup()  %>% group_by(module) %>%
    mutate(mo.r.abun = sum(r.abun), mo.num.otu = sum(num.otu), mo.mean.degree = sum(num.otu*mean.degree)/mo.num.otu,
           "num.otu*degree" = num.otu*mean.degree, "mo.num.otu*degree" = mo.num.otu*mo.mean.degree) %>% ungroup() -> otu.mo


  if(output == "r.abun"){y = "r.abun"; mo.y = "mo.r.abun"; mo.txt = paste0(formatC(100*otu.mo$mo.r.abun, format = "f", digits = 1), "%")}
  if(output == "num.otu"){y = "num.otu"; mo.y = "mo.num.otu"; mo.txt = otu.mo$mo.num.otu}
  if(output == "num.otu*degree"){y = "num.otu*degree"; mo.y = "mo.num.otu*degree"; mo.txt = paste0(otu.mo$mo.num.otu, "; ", formatC(otu.mo$mo.mean.degree, format = "f", digits = 0))}


  otu.mo %>% select(module, group, y = y, mo.y = mo.y) %>% mutate(mo.txt = mo.txt) %>%
    mutate(rank = 1+max(dense_rank(mo.y))-dense_rank(mo.y)) %>% dplyr::filter(rank <= module.n) -> df_p1

  df_p1 %>% select(module, rank) -> id
  while(length(unique(id$module)) > module.n){distinct(id) %>% arrange(rank) %>% slice_head(n = -1) -> id}
  arrange(id, rank) -> id; setdiff(unique(df_p1$module), unique(id$module)) -> mo.id.loss
  if(length(mo.id.loss > 0)){warning("0_0| module # ", paste0(mo.id.loss, collapse = ", "), " are lost due to duplication")}

  if(is.null(module.id)){unique(id$module) -> mo.id} else{module.id -> mo.id}
  mo.id %>% data.frame(module = .) %>% mutate(rank.id = c(1:length(module))) -> df_mo.id


  df_p1 %>% dplyr::filter(module %in%  mo.id) %>% merge(., df_mo.id, by = "module") %>%
    mutate(y1 = y/max(mo.y), mo.y1 = mo.y/max(mo.y)) %>%
    group_by(group) %>% mutate(y.sum.group = sum(y)) %>% ungroup() %>%
    mutate(group.p = paste0(gsub("\\(.*?\\)", "", group), "(", formatC(y.sum.group/sum(y)*100, format = "f", digits = 1), "%)")) %>%
    arrange(group) %>% mutate(group.p = factor(group.p, levels = unique(group.p))) -> df_p

  data.frame(g = unique(otu.g$v.group)) %>% mutate(col = v.cols[1:length(g)]) %>% dplyr::filter(g %in% df_p$group) %>% pull(col) -> v.cols

  ggplot(df_p, aes(x = as.factor(rank.id), y = y1, fill = group.p)) +
    geom_bar(stat = "identity", position = "stack") + scale_fill_manual(values = v.cols) +
    coord_radial(start = -.1 * pi, end = 1.7/7 * max(df_p$rank.id) * pi, inner.radius = 0) +
    annotate(geom = "text", x = df_p$rank.id, y = df_p$mo.y1 + .2, size = 2.8,
             label = paste0("# ", df_p$module, "\n(", df_p$mo.txt, ")"), angle = -(360*.87/7*(df_p$rank.id -.8))) +
    ylim(-.5, NA) + theme_void(10) + guides(fill = guide_legend(nrow = 2)) + labs(fill = '')
    theme(legend.position = "bottom", legend.direction = "horizontal", legend.key.size = unit(.18, "in")) %>%
    return()
}
