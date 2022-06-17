# loading libraries
library(tcR)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(ggbeeswarm)
library(rstatix)
library(igraph)

give.n <- function(y){ # adds number of points in a group the plot
  return(data.frame(
    y = 0,
    label = length(y))) 
}

# Loading samples ----
# NB: the names of samples MUST have the following structure: diagnosis_donorID_source_fraction_timepoint
# Allowed values:
# diagnosis: 'as' for ankylosing spondylitis; 'psa' for psoriatic arthritis, 'h' for healthy donors
# source: 'PB' for peripheral blood, 'SF' for synovial fluid; 
# fraction: - 'F' for total T cells; '4' for CD4+ T cells; '8' for CD8+ T cells
# timepoint: 'p1' or 'p2' (only for donors with recurrent synovitis)
#
# Examples: as_Bel_PB_CD8, psa_Stu_SF_F_p1
#
#
datalist = parse.folder(.folderpath = 'PATH-TO-MIXCR-EXPORTED-CLONESETS', .format = 'mixcr')

pb_f_inframes = get.inframes(datalist[grep('_PB_F', names(datalist))])
pb_4 = get.inframes(datalist[grep('_PB_4', names(datalist))])
pb_8 = get.inframes(datalist[grep('_PB_8', names(datalist))])
sf_f_single = get.inframes(datalist[grep('_SF_F', names(datalist))])
sf_f_single = sf_f_single[-grep('_p2$', names(sf_f_single))] # removes repertoires of the second timepoint for patients with recurring synovitis
spa_sf_4_clear = get.inframes(datalist[grep('_SF_4', names(datalist))])
spa_sf_8_clear = get.inframes(datalist[grep('_SF_8', names(datalist))])

# 1. Diversity ----
pb_has_34k_umi = which(sapply(c(pb_f_inframes, pb_4, pb_8), function(x) sum(x[,'Read.count']))>=34000)
pb_34k = clonotype.sampler.new(c(pb_f_inframes, pb_4, pb_8)[pb_has_34k_umi], 'Read.count', .number = 34000)
pb_34k = lapply(pb_34k, function(x) {x$Read.proportion= x$Read.count/sum(x$Read.count);x})

div_pb_34k = data.frame(clones_strict = sapply(pb_34k, nrow)) %>%
  mutate(gini = sapply(pb_34k, function(x) gini(x$Read.proportion))) %>%
  mutate(sample = names(pb_34k)) %>% 
  separate(sample, c('status', 'donor','source', 'fraction'), '_') %>%
  gather('metric', 'value', 1:3) %>%
  mutate(metric = factor(metric, levels = c('clones_strict', 'gini')))

sf_has_34k_umi = which(sapply(c(sf_f_single, spa_sf_4_clear, spa_sf_8_clear), function(x) sum(x[,'Read.count']))>=34000)
sf_34k = clonotype.sampler.new(c(sf_f_single, spa_sf_4_clear, spa_sf_8_clear)[sf_has_34k_umi], 'Read.count', .number = 34000)
sf_34k = lapply(sf_34k, function(x) {x$Read.proportion= x$Read.count/sum(x$Read.count);x})

div_sf_34k = data.frame(clones_strict = sapply(sf_34k, nrow)) %>%
  mutate(gini = sapply(sf_34k, function(x) gini(x$Read.proportion))) %>%
  mutate(sample = names(sf_34k)) %>% 
  separate(sample, c('status', 'donor','source', 'fraction'), '_') %>%
  gather('metric', 'value', 1:3) %>%
  mutate(metric = factor(metric, levels = c('clones_strict', 'gini')))


div_34k = rbind(div_pb_34k, div_sf_34k) %>% 
  filter(status %in% c('as', 'psa')) %>% 
  mutate(fraction = factor(fraction, levels=c('4','8','F'), labels = c('CD4','CD8','F')))

## Fig 1A ----

div_paired = div_34k %>% 
  filter(metric=='clones_strict') %>% 
  group_by(fraction) %>% 
  count(donor) %>% 
  filter(n==2)

plot_div_clones_strict = ggboxplot(div_34k %>% inner_join(div_paired) %>% filter(metric=='clones_strict'), 
                                   x='source', y='value', color='source',
                                   xlab = F, ylab='# clonotypes', width = 0.5,
                                   add='point', add.params = list(alpha=0.3, size=1, position=position_quasirandom(0.2)),
                                   palette = 'npg', outlier.shape = NA) + 
  facet_wrap(~fraction, ncol = 1, labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  stat_compare_means(method = 'wilcox.test', paired = T, label = "p.signif", tip.length = 0, hide.ns = T, step.increase = 0, 
                     label.x.npc = 'middle') + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        legend.title = element_blank(),
        panel.spacing = unit(2, "lines")) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = 0))


## Fig 1B ----
plot_div_gini = ggboxplot(div_34 %>% inner_join(div_paired) %>% filter( metric=='gini'), 
                          x='source', y='value', color='source',
                          xlab = F, ylab='Index Gini', width = 0.5,
                          # title = 'gini',
                          add='point', add.params = list(alpha=0.3, size=1, position=position_quasirandom(0.2)),
                          palette = 'npg', outlier.shape = NA) +
  facet_wrap(~fraction, ncol = 1, labeller = label_wrap_gen(multi_line=FALSE)) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  stat_compare_means(method = 'wilcox.test', paired = T, label = "p.signif", tip.length = 0, hide.ns = T, step.increase = 0, 
                     label.x.npc = 'middle') + 
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(2, "lines")) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = 0))

# 2. Cross ----
hla.crs = read.table(file = 'hla.crs.txt', header = T, sep = '\t', stringsAsFactors = F)

crs_norm_sf_8 = intersectClonesets(sf_34k[grep('_8', names(sf_34k))], .type = 'ave', .norm = T)
crs_norm_sf_8[lower.tri(crs_norm_sf_8, diag=T)] = NA
crs_norm_sf_8 = crs_norm_sf_8 %>% 
  as.data.frame() %>%
  mutate(sample1=names(.)) %>%
  gather('sample2', 'value', 1:(length(.)-1)) %>%
  separate(sample1, c('status1', 'donor1', 'source1', 'fraction1'), '_') %>%
  separate(sample2, c('status2', 'donor2', 'source2', 'fraction2'), '_') %>% 
  mutate(group=paste(status1, status2, sep=' x ')) %>% 
  left_join(hla.crs, by = c('donor1'='s1', 'donor2'='s2'))

crs_norm_sf_4 = intersectClonesets(sf_34k[grep('_4', names(sf_34k))], .type = 'ave', .norm = T)
crs_norm_sf_4[lower.tri(crs_norm_sf_4, diag=T)] = NA
crs_norm_sf_4 = crs_norm_sf_4 %>% 
  as.data.frame() %>%
  mutate(sample1=names(.)) %>%
  gather('sample2', 'value', 1:(length(.)-1)) %>%
  separate(sample1, c('status1', 'donor1', 'source1', 'fraction1'), '_') %>%
  separate(sample2, c('status2', 'donor2', 'source2', 'fraction2'), '_') %>% 
  mutate(group=paste(status1, status2, sep=' x ')) %>% 
  left_join(hla.crs, by = c('donor1'='s1', 'donor2'='s2'))

crs_norm_sf_f = intersectClonesets(sf_34k[grep('_F', names(sf_34k))], .type = 'ave', .norm = T)
crs_norm_sf_f[lower.tri(crs_norm_sf_f, diag=T)] = NA
crs_norm_sf_f = crs_norm_sf_f %>% 
  as.data.frame() %>%
  mutate(sample1=names(.)) %>%
  gather('sample2', 'value', 1:(length(.)-1)) %>%
  separate(sample1, c('status1', 'donor1', 'source1', 'fraction1'), '_') %>%
  separate(sample2, c('status2', 'donor2', 'source2', 'fraction2'), '_') %>% 
  mutate(group=paste(status1, status2, sep=' x ')) %>% 
  left_join(hla.crs, by = c('donor1'='s1', 'donor2'='s2'))


crs_norm_pb_sf = rbind(crs_norm_pb_f,crs_norm_pb_4,crs_norm_pb_8,
                       crs_norm_sf_f,crs_norm_sf_4,crs_norm_sf_8) %>%
  filter(status1 !='h', status2 !='h') %>%
  mutate(fraction = factor(fraction, levels=c('F', '8', '4'), labels = c('F', 'CD8', 'CD4')))

crs_norm_pb_sf_median = crs_norm_pb_sf %>%
  group_by(source, fraction, donor1) %>% 
  summarise(med_d1 = median(value)) 

## Fig 1C ----
plot_crs_norm_pb_sf_median = ggboxplot(crs_norm_pb_sf_median, 
                                       x="source", y="med_d1", color="source",
                                       ylab='Normalized number of shared clonotypes', width = 0.5,
                                       add = 'point', add.params = list(alpha=0.2, size=1, position=position_quasirandom(0.2)),
                                       palette = 'npg', outlier.shape = NA) +
  stat_compare_means(label = "p.signif", tip.length = 0, hide.ns = T, step.increase = 0, 
                     label.x.npc = 'middle') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + 
  facet_wrap(~fraction, ncol = 1) + 
  rremove('legend') + rremove('xlab') +
  theme(strip.background = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(2, "lines")) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = 0))

## Fig 1D ----
crs_norm = intersectClonesets(c(pb_34k[grep('^as_|^psa_', names(pb_34k))], sf_34k), .type = 'ave', .norm = T) %>%
  as.data.frame() %>%
  mutate(sample1=names(.)) %>%
  gather('sample2', 'value', 1:(length(.)-1)) %>%
  na.omit() %>%
  separate(sample1, c('status1', 'donor1', 'source1', 'fraction1'), '_') %>%
  separate(sample2, c('status2', 'donor2', 'source2', 'fraction2'), '_') %>% 
  mutate(group1=paste(source1, source2, sep=' x '),
         group2=paste(fraction1, fraction2, sep=' x ')) %>%
  left_join(hla.crs, by = c('donor1'='s1', 'donor2'='s2')) %>%
  filter(group1 %in% c('PB x PB', 'SF x SF'), 
         group2 %in% c('F x F', '4 x 4', '8 x 8'))



crs_norm_f = crs_norm %>% 
  filter(fraction1=='F') %>% 
  group_by(hla.total, source1, donor1) %>% 
  summarise(med_d1 = median(value)) 

plot_crs_norm_f = ggboxplot(crs_norm_f %>%
                              filter(hla.total<=6), x="source1", y="med_d1", color="source1",
                            xlab = '', ylab='', title = '# matched HLA alleles',
                            outlier.shape = NA, add = 'point', add.params = list(alpha=0.3, size=1, position=position_quasirandom(0.2)),
                            palette = 'npg') +
  stat_compare_means(label = "p.signif", tip.length = 0, hide.ns = T, step.increase = 0, 
                     label.x.npc = 'middle') + 
  theme(title = element_text(size = 10)) +
  scale_y_log10(expand = expansion(mult = c(0.1, 0.2)), n.breaks=4) + 
  rremove('legend') + rremove('xylab') +
  facet_wrap(.~hla.total, ncol=7) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = log(5e-08,10))) +
  font('x.text', color='black', size=10)


crs_norm_8 = crs_norm %>% 
  filter(fraction1=='8') %>% 
  group_by(hla1, source1, donor1) %>% 
  summarise(med_d1 = median(value)) 

plot_crs_norm_8 = ggboxplot(crs_norm_8 %>%
                              filter(hla1<=3), x="source1", y="med_d1", color="source1",
                            xlab = '', ylab='', title = '# matched HLA-I alleles',
                            outlier.shape = NA, add = 'point', add.params = list(alpha=0.3, size=1, position=position_quasirandom(0.2)),
                            palette = 'npg') +
  stat_compare_means(label = "p.signif", tip.length = 0, hide.ns = T, step.increase = 0, 
                     label.x.npc = 'middle') + 
  theme(title = element_text(size = 10)) +
  scale_y_log10(expand = expansion(mult = c(0.1, 0.2)), n.breaks = 4) + 
  rremove('x.axis') + rremove('x.text') + rremove('x.ticks') + rremove('legend') + rremove('xylab') +
  facet_wrap(.~hla1, ncol=7) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = log(3e-08,10)))

crs_norm_4 = crs_norm %>% 
  filter(fraction1=='4') %>% 
  group_by(hla2, source1, donor1) %>% 
  summarise(med_d1 = median(value)) 

plot_crs_norm_4 = ggboxplot(crs_norm_4 %>%
                              filter(hla2<=3), x="source1", y="med_d1", color="source1",
                            xlab = '', ylab='', title = '# matched HLA-II alleles',
                            outlier.shape = NA, add = 'point', add.params = list(alpha=0.3, size=1, position=position_quasirandom(0.2)),
                            palette = 'npg') +
  stat_compare_means(label = "p.signif", tip.length = 0, hide.ns = T, step.increase = 0, 
                     label.x.npc = 'middle') + 
  theme(title = element_text(size = 10)) +
  scale_y_log10(expand = expansion(mult = c(0.1, 0.2)), n.breaks = 4) + 
  rremove('legend') + rremove('xylab') +
  facet_wrap(.~hla2, ncol=7) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = log(5e-08, 10))) 


plot_crs_norm = ggarrange(plot_crs_norm_f,
                          plot_crs_norm_8, 
                          plot_crs_norm_4, 
                          nrow = 3, ncol = 1, legend = 'none', align = 'h') 

plot_crs_norm = annotate_figure(plot_crs_norm, left = 'Normalized number of shared clonotypes')


# 3. ALICE for SF CD8+ TCRbeta ----
sf8beta = get.inframes(spa_sf_8_clear)
sf8beta = lapply(sf8beta, function(x) {x = x %>% 
  filter(Read.count>1) %>%
  rename(bestVGene = V.gene, bestJGene = J.gene) %>%
  as.data.table
}
)

sf8beta_alice<-ALICE_pipeline(DTlist=sf8beta,folder="sf_cd8_b",cores=7,iter=25,nrec=5e6, 
                              Read_count_filter = 0)
sf8beta_alice = lapply(sf8beta_alice, function(x) x= as.data.frame(x))

q = replicate(length(table(sh_sf8beta_alice$bestVGene)), paste(sample(LETTERS, 5, replace=TRUE), collapse=""))
min(stringdist::stringdistmatrix(q, method = 'hamming'))
vgenes = as.data.frame(table(sh_sf8beta_alice$bestVGene), stringsAsFactors = F)
vgenes$chifr = q
for (i in 1:nrow(vgenes)){
  sh_sf8beta_alice[grep(paste(vgenes[i,1],'$', sep = ''), sh_sf8beta_alice$bestVGene),1] = paste(sh_sf8beta_alice[grep(paste(vgenes[i,1],'$', sep = ''), sh_sf8beta_alice$bestVGene),1],vgenes[i,3], sep = '')
}
gr_sf8beta_alice = mutation.network(sh_sf8beta_alice, .method = 'hamm', .max.errors = 1, 
                                    .label.col = 'CDR3.amino.acid.sequence', .seg.col = 'bestVGene')
gr_sf8beta_alice = set.vertex.attribute(gr_sf8beta_alice, name = 'label2', index = V(gr_sf8beta_alice), 
                                        value = substr(sh_sf8beta_alice$CDR3.amino.acid.sequence, start = 0, stop = nchar(sh_sf8beta_alice$CDR3.amino.acid.sequence)-5))
sh_sf8beta_alice = sh_sf8beta_alice %>% 
  mutate(CDR3.amino.acid.sequence = substr(CDR3.amino.acid.sequence, start = 0, stop = nchar(CDR3.amino.acid.sequence)-5))

# Clustering 
cluster = clusters(gr_sf8beta_alice) 
sh_sf8beta_alice$cluster = cluster$membership 

# Selected clusters graph
sf8beta_clusters = igraph::groups(cluster)

a = list()
for (i in 1:length(igraph::groups(cluster))) {
  clst = V(gr_sf8beta_alice)$people[igraph::groups(cluster)[[i]]]
  a[[i]] = rowSums(sapply(clst, function(x) as.numeric(unlist(str_split(x, pattern = '')))))
}
ppl = sapply(a, function(x) sum(x>0)) # number of donors that share a luster

sf8beta_clusters = lapply(sf8beta_clusters, function(x) x=V(gr_sf8beta_alice)$label2[x])
sf8beta_clusters = sf8beta_clusters[ppl>2] # select clusters shared by at least 3 donors

gr_sf8beta_alice = set.vertex.attribute(gr_sf8beta_alice, name = 'cluster_id', index = V(gr_sf8beta_alice), 
                                        value = cluster$membership)
gr_gut = subgraph(gr_sf8beta_alice, v = which(V(gr_sf8beta_alice)$cluster_id %in% as.numeric(names(sf8beta_clusters))))

selected = sh_sf8beta_alice %>% 
  filter(cluster %in% V(gr_gut)$cluster_id) %>%
  relocate(cluster, .after=People)

# 4. Annotation of clusters with VDJdb----
vdjdb = read.delim("vdjdb.slim_2021-02-02.txt", header=T)
vdjdb_in_selected = selected[,1:41] %>% inner_join(vdjdb, 
                                                   by=c('CDR3.amino.acid.sequence'='cdr3', 'bestVGene'='v.segm'))

vdjdb_in_selected.good <- vdjdb_in_selected %>%
  mutate(mhc.a = str_split_fixed(mhc.a, "[:,]", 2)[,1]) %>%
  select(CDR3.amino.acid.sequence, antigen.epitope, mhc.a, vdjdb.score, reference.id) %>% 
  unique %>%
  #--- filtering of VDJdb ---
  group_by(CDR3.amino.acid.sequence) %>%
  mutate(vdjdb.score.max = max(vdjdb.score)) %>%
  filter(vdjdb.score == vdjdb.score.max) %>%
  group_by(CDR3.amino.acid.sequence) %>%
  mutate(num.pub = str_count(reference.id, ","),
         num.pub.max = max(num.pub)) %>%
  filter(num.pub == num.pub.max) %>%
  # Remove all remaining ambigous cases
  group_by(CDR3.amino.acid.sequence) %>%
  mutate(n.epitopes = length(unique(antigen.epitope))) %>% 
  filter(n.epitopes == 1) %>%
  ungroup

vdjdb_in_selected <- vdjdb_in_selected %>% 
  merge(vdjdb_in_selected.good, by=c('CDR3.amino.acid.sequence', 'antigen.epitope'))


# 5. B27 vs B38 enrichment----
b2705 = c("Abr","Ali","Ash-110","Bal","Bel","Bost","Dv","Evst","Gar","GE","Gonch","Kos",
"Kzh","Leb","LK","LN","Mart","Mikh","Mor","Non","Nos","OB","Pach","Pon","Sass","Sem","Shep",
"Shor","Luk","Stu","SZ","TumE","TumV","TwHM","Uv","Vats","Vol","Zakh","Zar","Zbot","ZK","Zlt","ZN")

b3801 = c("Azh","YZh","Eluk","EL","Kal","Lem","MStr","Pet","Rod","Sbkv","Luk","Uma","Vor","Zar")

selected_clusters_summary = selected %>% 
  select(CDR3.amino.acid.sequence, bestVGene, cluster_id) %>%
  rename(V.gene = bestVGene, cluster = cluster_id) %>%
  left_join(sh_sf_8_beta, by = c("CDR3.amino.acid.sequence", "V.gene" )) %>%
  group_by(cluster) %>% 
  summarise(across(c(as_Bel_SF_8:psa_Zlt_SF_8_p1), ~sum(.x, na.rm = T))) %>%
  select(-contains('Zar')) %>% ## IMPORTANT, removes HLA-B*27+B*38+ patient
  rowwise() %>% 
  mutate(people = rowSums(across(as_Bel_SF_8:psa_Zlt_SF_8_p1)>0)) %>%
  ungroup %>%
  mutate(b27pos = rowSums(.[,grep(paste(b2705, collapse = '|'), names(.))]>0),
         b38pos = rowSums(.[,grep(paste(b3801, collapse = '|'), names(.))]>0)) %>%
  mutate(total_b27pos = length(grep(paste(b2705, collapse = '|'), names(.))),
         total_b38pos = length(grep(paste(b3801, collapse = '|'), names(.)))) %>%
  mutate(rest_b27 = total_b27pos-b27pos, 
         rest_b38=total_b38pos-b38pos) %>% 
  # mutate(n.clonotypes = selected %>% count(cluster_id) %>% select(n)) %>%
  relocate(people:total_b38pos, .after=cluster) %>%
  mutate(cluster = factor(cluster, levels=c(1,10,12,404,4,6,15,16,701)))

fisher.p = c()  
for (n in 1:nrow(selected_clusters_summary)){
  print(n)
  cln = matrix(unlist(selected_clusters_summary[n,c('b27pos','b38pos','rest_b27','rest_b38')]), ncol=2, 
               dimnames = list(c('b28pos','b38pos'),c('pos','neg')))
  fisher.p[n] = pairwise_fisher_test(cln)$p
}

selected_clusters_summary = selected_clusters_summary %>% 
  mutate(fisher.p = fisher.p) %>%
  mutate(fisher.p.adj = p.adjust(fisher.p, method = 'fdr'))

## Fig 2B ----
plot_selected_clusters_b27_vs_b38 = ggplot(selected_clusters_summary, aes(x=b27pos/total_b27pos, y=b38pos/total_b38pos, 
                                                                          fill=cluster)) + 
  geom_abline(slope = 1, intercept = 0, color='gray80') +
  geom_point(position=position_jitter(0.02,0.02, seed = 9), size=3, shape=21, alpha=0.9) + 
  theme_classic() + 
  scale_fill_d3(guide=F, na.value = 'gray90', drop=F) + 
  coord_fixed(ratio = 1, ylim = c(0,1), xlim = c(0,1)) +
  xlab('n=15 \n Occurrence in HLA-B*27+ donors, \n proportion') +
  ylab('Occurrence in HLA-B*38+ donors, \n proportion \n n=7') +
  font('xy.text', color='black', size=12) +
  font('xylab', color='black', size=12)

# 6. Enrichment of clusters in CD8+ SF vs CD8+ PB ----

# selected_clusters_summary_enriched = selected_clusters_summary %>% 
#   filter(fisher.p.adj<0.05)

# selected_clusters_summary_enriched %>% 
#   gather(as_Bel_SF_8:psa_Zlt_SF_8_p1, key = 'sample', value = 'proportion', na.rm = F) %>% 
#   separate(sample, c('status', 'donor','source', 'fraction'), '_') %>% 
#   filter(proportion>0) %>% 
#   mutate(b27 = if_else(donor %in% b2705, 'b27', 'NA'),
#          b38 = if_else(donor %in% b3801, 'b38', 'NA')) %>%
#   group_by(cluster) %>% 
#   summarise(b27_as = n_distinct(donor[status=='as' & b27=='b27']),
#             b27_psa = n_distinct(donor[status=='psa' & b27=='b27']),
#             b38_as = n_distinct(donor[status=='as' & b38=='b38']),
#             b38_psa = n_distinct(donor[status=='psa' & b38=='b38']),
#             other_psa = n_distinct(donor[status=='psa' & b27!='b27' & b38!='b38'])) %>%
#   mutate(b27_total = b27_as + b27_psa,
#          b38_total = b38_as + b38_psa)

prop_selected_clusters_summary_enriched_in_pb_8 = sh_pb_8_beta %>% 
  right_join(enriched_clusters_cumulative, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene')) %>%
  gather(`as_Ash-110_PB_8`:psa_Zlt_PB_8, key = 'sample', value = 'proportion', na.rm = F)

prop_selected_clusters_summary_enriched_in_sf_8 = sh_sf_8_beta %>% 
  right_join(enriched_clusters_cumulative, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene')) %>%
  gather(as_Bel_SF_8:psa_Zlt_SF_8_p1, key = 'sample', value = 'proportion', na.rm = F)

prop_selected_clusters_summary_enriched_pb8_vs_sf8 = rbind(prop_selected_clusters_summary_enriched_in_pb_8, prop_selected_clusters_summary_enriched_in_sf_8) %>%
  separate(sample, c('status', 'donor','source', 'fraction'), '_') %>%
  mutate(cluster = factor(cluster, levels = c(1,10,12,404,4,6,15,16,701)))

stat_selected_clusters_summary_enriched_in_pb8_vs_sf8 = prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>% 
  group_by(cluster,source, donor) %>% 
  summarise(cum_prop = sum(proportion, na.rm = T)) %>%
  # mutate(bts.pval = if_else(cluster %in% selected_clusters_summary_enriched_b27, 'b27', 'b38')) %>%
  mutate(cluster = factor(cluster, levels = c(1,10,12,404,4,6,15,16,701)))

stat.test_pb_sf_b27 = stat_selected_clusters_summary_enriched_in_pb8_vs_sf8 %>%
  filter(donor %in% b2705) %>% 
  filter(cluster %in% c(1,10,12,404))  %>%
  group_by(cluster) %>% 
  wilcox_test(cum_prop ~ source) %>% 
  mutate(p.adj = p.adjust(p, 'fdr')) %>% 
  add_significance('p.adj') %>%
  mutate(y.position = 0.01)

## Fig 2C ----
plot_selected_clusters_summary_enriched_pb8_vs_sf8_b27 = ggboxplot(stat_selected_clusters_summary_enriched_in_pb8_vs_sf8 %>%
                                                                     filter(donor %in% b2705) %>% 
                                                                     filter(cluster %in% c(1,10,12,404)) %>%
                                                                     mutate(cum_prop = cum_prop+1e-8), 
                                                                   x='source', y='cum_prop',
                                                                   color='cluster', 
                                                                   # ylab='Cumulative proportion',
                                                                   # title ='B27+ donors' ,
                                                                   add = 'point', outlier.shape = NA,
                                                                   add.params = list(alpha=0.5, size=2,
                                                                                     position = position_quasirandom(0.2))) +
  facet_grid(~cluster) +
  scale_color_d3(drop=F,  guide=F) +
  scale_y_log10(expand = expansion(mult = c(0.01, 0.2)), 
                labels = trans_format(trans = 'log10',format = math_format())) + 
  stat_pvalue_manual(stat.test_pb_sf_b27, label = "p.adj.signif",
                     tip.length = 0, hide.ns = T, step.increase = 0) + 
  rremove('xylab') +
  font("xy.text", size = 10)


stat.test_pb_sf_b38 = stat_selected_clusters_summary_enriched_in_pb8_vs_sf8 %>%
  filter(donor %in% b3801) %>% 
  filter(cluster %in% c(4,6,15,16,701))  %>%
  group_by(cluster) %>% 
  wilcox_test(cum_prop ~ source) %>% 
  mutate(p.adj = p.adjust(p, 'fdr')) %>% 
  add_significance('p.adj') %>%
  mutate(y.position = 0.01)

plot_selected_clusters_summary_enriched_pb8_vs_sf8_b38 = ggboxplot(stat_selected_clusters_summary_enriched_in_pb8_vs_sf8 %>%
                                                                     filter(donor %in% b3801) %>% 
                                                                     filter(cluster %in% c(4,6,15,16,701)) %>%
                                                                     mutate(cum_prop = cum_prop+1e-8), 
                                                                   x='source', y='cum_prop',
                                                                   color='cluster', 
                                                                   # ylab='Cumulative proportion',
                                                                   # title ='B38+ donors' , 
                                                                   add = 'point', outlier.shape = NA,
                                                                   add.params = list(alpha=0.3, size=2,
                                                                                     position = position_quasirandom(0.2))) +
  facet_grid(~cluster) +
  scale_color_d3(drop=F,  guide=F) +
  scale_y_log10(expand = expansion(mult = c(0.01, 0.2)), 
                labels = trans_format(trans = 'log10',format = math_format())) + 
  stat_pvalue_manual(stat.test_pb_sf_b38, label = "p.adj.signif",
                     tip.length = 0, hide.ns = T, step.increase = 0) + 
  rremove('xylab') +
  font("xy.text", size = 10) 


plot_selected_clusters_summary_enriched_pb8_vs_sf8 = ggarrange(plot_selected_clusters_summary_enriched_pb8_vs_sf8_b27, 
                                                               plot_selected_clusters_summary_enriched_pb8_vs_sf8_b38,
                                                               ncol = 1, nrow = 2, align = 'h')

plot_selected_clusters_summary_enriched_pb8_vs_sf8 = annotate_figure(plot_selected_clusters_summary_enriched_pb8_vs_sf8, 
                                                                     left = 'Cumulative proportion', bottom = '')

# 7. Occurrence in total PB repertoires of SpA patients and healthy donors ----
# 'emerson' is a dataset from Emerson et al. 2017 Nat Gen, available from https://clients.adaptivebiotech.com/pub/Emerson-2017-NatGen.
emerson = lapply(emerson, function(x) { x = x[,1:11]; names(x) = names(healthy_f[[1]])[c(1:5,7:8,11:14)]; x})
emerson_b27 = emerson_metadata[grep('HLA-B*27', emerson_metadata$hla, fixed = T),]$sample_id
emerson_b38 = emerson_metadata[grep('HLA-B*38', emerson_metadata$hla, fixed = T),]$sample_id

healthy_f = lapply(healthy_f, function(x) {x = x %>% inner_join(sh_sf8beta_alice[,c(1:3,5)], 
                                                                by = c('CDR3.amino.acid.sequence', 
                                                                       'V.gene' = 'bestVGene'))})
enriched_clusters_cumulative = sh_sf8beta_alice %>% 
  filter(cluster %in% enriched_clusters_cumulative$cluster) %>%
  select(CDR3.amino.acid.sequence, bestVGene, bestJGene, cluster) %>%
  filter(!(cluster==6 & bestJGene=='TRBJ2-5'))

selected_clusters_summary_enriched_in_pb_hd = shared.repertoire(c(emerson, healthy_f), .type = 'avrp') %>%
  right_join(enriched_clusters_cumulative, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene')) %>%
  relocate(cluster, .after=V.gene) %>%
  group_by(cluster) %>% 
  summarise(across(HIP00110.tsv.mixcr.gz:h_ZN_PB_F, ~sum(.x, na.rm=T))) %>%
  mutate(b27pos = rowSums(.[,grep(paste(c(emerson_b27, b2705), collapse = '|'), names(.))]>0),
         b38pos = rowSums(.[,grep(paste(c(emerson_b38, b3801),collapse = '|'), names(.))]>0),
         total_b27pos = length(grep(paste(c(emerson_b27, b2705), collapse = '|'), names(.))),
         total_b38pos = length(grep(paste(c(emerson_b38, b3801), collapse = '|'), names(.)))) %>%
  relocate(b27pos:total_b38pos, .after=cluster)

selected_clusters_summary_enriched_in_pb_spa = sh_spa_pb_f %>%
  right_join(enriched_clusters_cumulative, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene')) %>%
  relocate(cluster, .after=V.gene) %>%
  group_by(cluster) %>% 
  summarise(across(as_Abd_PB_F:psa_Uma_PB_F, ~sum(.x, na.rm=T))) %>%
  mutate(b27pos = rowSums(.[,grep(paste(c(emerson_b27, b2705), collapse = '|'), names(.))]>0),
         b38pos = rowSums(.[,grep(paste(c(emerson_b38, b3801),collapse = '|'), names(.))]>0),
         total_b27pos = length(grep(paste(c(emerson_b27, b2705), collapse = '|'), names(.))),
         total_b38pos = length(grep(paste(c(emerson_b38, b3801), collapse = '|'), names(.)))) %>%
  relocate(b27pos:total_b38pos, .after=cluster)

selected_clusters_summary_enriched_in_pb = selected_clusters_summary_enriched_in_pb_hd[,1:5] %>%
  full_join(selected_clusters_summary_enriched_in_pb_spa[,1:5], suffix = c('.h', '.spa'), by = 'cluster') %>%
  mutate(rest_b27pos.spa = total_b27pos.spa-b27pos.spa,
         rest_b38pos.spa = total_b38pos.spa-b38pos.spa,
         rest_b27pos.h = total_b27pos.h-b27pos.h,
         rest_b38pos.h = total_b38pos.h-b38pos.h) %>%
  mutate(cluster = factor(cluster, levels = c(1,10,12,404,4,6,15,16,701)))

fisher.p.b27 = c()  
for (n in 1:nrow(selected_clusters_summary_enriched_in_pb)){
  cln = matrix(unlist(selected_clusters_summary_enriched_in_pb[n,c('b27pos.spa','b27pos.h','rest_b27pos.spa','rest_b27pos.h')]), ncol=2, 
               dimnames = list(c('spa','h'),c('rest.spa','rest.h')))
  fisher.p.b27[n] = pairwise_fisher_test(cln)$p
}


fisher.p.b38 = c()  
for (n in 1:nrow(selected_clusters_summary_enriched_in_pb)){
  cln = matrix(unlist(selected_clusters_summary_enriched_in_pb[n,c('b38pos.spa','b38pos.h','rest_b38pos.spa','rest_b38pos.h')]), ncol=2, 
               dimnames = list(c('spa','h'),c('rest.spa','rest.h')))
  fisher.p.b38[n] = pairwise_fisher_test(cln)$p
}


selected_clusters_summary_enriched_in_pb = selected_clusters_summary_enriched_in_pb %>%
  mutate(fisher.p.b27 = fisher.p.b27,
         fisher.p.b38 = fisher.p.b38) %>%
  mutate(fisher.p.adj.b27 = p.adjust(fisher.p.b27, method = 'fdr'),
         fisher.p.adj.b38 = p.adjust(fisher.p.b38, method = 'fdr'))

## Fig 2D ----
plot_summary_clusters_b27 = ggplot(selected_clusters_summary_enriched_in_pb %>% filter(cluster %in% c(1,10,12,404)), 
                                   aes(x = b27pos.spa/total_b27pos.spa, 
                                       y=b27pos.h/total_b27pos.h, 
                                       fill=cluster,
                                       stroke = if_else(fisher.p.adj.b27<0.05, 1.5, 0.3))) + 
  geom_abline(slope = 1, intercept = 0, color='gray80') +
  geom_point(position=position_jitter(0.01,0.01), size=3, shape=21, alpha=0.8) + 
  theme_classic() + 
  scale_fill_d3(guide=F, drop=F) + 
  coord_fixed(ratio = 1, ylim = c(0,1), xlim = c(0,1)) +
  # ggtitle('B27+ donors') +
  xlab('n=29') + ylab('n=51') +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  font("xy.text", size = 12, color='black')

plot_summary_clusters_b38 = ggplot(selected_clusters_summary_enriched_in_pb %>% filter(cluster %in% c(4,6,15,16,701)), 
                                   aes(x = b38pos.spa/total_b38pos.spa, 
                                       y=b38pos.h/total_b38pos.h, 
                                       fill=cluster,
                                       stroke = if_else(fisher.p.adj.b38<0.05, 1.5, 0.3))) + 
  geom_abline(slope = 1, intercept = 0, color='gray80') +
  geom_point(position=position_jitter(0.01,0.01), size=3, shape=21, alpha=0.8) + 
  theme_classic() + 
  scale_fill_d3(guide=F, drop=F) + 
  coord_fixed(ratio = 1, ylim = c(0,1), xlim = c(0,1)) +
  # ggtitle('B38+ donors') +
  xlab('n=14') + ylab('n=37') +
  theme(plot.margin = unit(c(0,0,0,0), "lines")) +
  font("xy.text", size = 12, color='black')


plot_summary_clusters = ggarrange(plot_summary_clusters_b27, plot_summary_clusters_b38, ncol = 1, nrow = 2, align = 'hv')
plot_summary_clusters = annotate_figure(plot_summary_clusters, 
                                        left = 'Occurrence in healthy donors, proportion', 
                                        bottom = 'Occurrence in SpA patients, \n proportion')

# 8. Proportion of clusrers in total PB repertoires of SpA patients and healthy donors----
prop_selected_clusters_summary_enriched_in_pb_f_b27 = selected_clusters_summary_enriched_in_pb_spa[,-c(2:5)] %>% 
  full_join(selected_clusters_summary_enriched_in_pb_hd[,-c(2:5)], by = 'cluster') %>%
  gather(key = 'sample', value = 'cum_prop', as_Abd_PB_F:b27_ZN_PB_F) %>%
  filter(grepl(paste(c(b2705, emerson_b27), collapse = '|'), sample)) %>%
  mutate(status = if_else(sample %in% grep('^as_|^psa_', sample, value = T),'SpA', 'HD')) %>%
  mutate(cluster = factor(cluster, levels = c(1,10,12,404,4,6,15,16,701)))

prop_selected_clusters_summary_enriched_in_pb_f_b38 = selected_clusters_summary_enriched_in_pb_spa[,-c(2:5)] %>% 
  full_join(selected_clusters_summary_enriched_in_pb_hd[,-c(2:5)], by = 'cluster') %>%
  gather(key = 'sample', value = 'cum_prop', as_Abd_PB_F:b27_ZN_PB_F) %>%
  filter(grepl(paste(c(b3801, emerson_b38), collapse = '|'), sample)) %>%
  mutate(status = if_else(sample %in% grep('^as_|^psa_', sample, value = T),'SpA', 'HD')) %>%
  mutate(cluster = factor(cluster, levels = c(1,10,12,404,4,6,15,16,701)))

## Fig 2E ----
stat.test_pb_f_b27 = prop_selected_clusters_summary_enriched_in_pb_f_b27 %>% 
  filter(cluster %in% c(1,10,12,404)) %>%
  group_by(cluster) %>% 
  wilcox_test(cum_prop ~ status) %>% 
  mutate(p.adj = p.adjust(p, 'fdr')) %>% 
  add_significance('p.adj') %>%
  mutate(y.position = 0.0001)

plot_prop_selected_clusters_summary_enriched_in_pb_f_b27 = ggboxplot(prop_selected_clusters_summary_enriched_in_pb_f_b27 %>% 
                                                                       filter(cluster %in% c(1,10,12,404)) %>%
                                                                       mutate(cum_prop = cum_prop+1e-8),  
                                                                     x='status', y='cum_prop',
                                                                     color='cluster', 
                                                                     # ylab='Cumulative proportion in PB',
                                                                     # title ='B27+ donors' ,
                                                                     add = 'point', outlier.shape = NA,
                                                                     add.params = list(alpha=0.3, size=2,
                                                                                       position = position_quasirandom(0.2))) +
  scale_color_d3(drop=F, guide=F) +
  facet_grid(~cluster) +
  scale_y_log10(expand = expansion(mult = c(0.01, 0.2)), 
                labels = trans_format(trans = 'log10',format = math_format())) + 
  stat_pvalue_manual(stat.test_pb_f_b27, label = "p.adj.signif", 
                     tip.length = 0, hide.ns = T, step.increase = 0, label.x.npc = 'center') +
  rremove('xylab') + 
  font("xy.text", size = 10)


stat.test_pb_f_b38 = prop_selected_clusters_summary_enriched_in_pb_f_b38 %>% 
  filter(cluster %in% c(4,6,15,16,701)) %>%
  group_by(cluster) %>% 
  wilcox_test(cum_prop ~ status) %>% 
  mutate(p.adj = p.adjust(p, 'fdr')) %>% 
  add_significance('p.adj') %>%
  mutate(y.position = 0.01)

plot_prop_selected_clusters_summary_enriched_in_pb_f_b38 = ggboxplot(prop_selected_clusters_summary_enriched_in_pb_f_b38 %>% 
                                                                       filter(cluster %in% c(4,6,15,16,701)) %>%
                                                                       mutate(cum_prop = cum_prop+1e-8), 
                                                                     x='status', y='cum_prop',
                                                                     color='cluster', 
                                                                     # ylab='Cumulative proportion in PB',
                                                                     # title ='B38+ donors' ,
                                                                     add = 'point', outlier.shape = NA,
                                                                     add.params = list(alpha=0.3, size=2,
                                                                                       position = position_quasirandom(0.2))) +
  scale_color_d3(drop=F, guide=F) +
  facet_grid(~cluster) +
  scale_y_log10(expand = expansion(mult = c(0.01, 0.2)), 
                labels = trans_format(trans = 'log10',format = math_format())) + 
  stat_pvalue_manual(stat.test_pb_f_b38, label = "p.adj.signif", 
                     tip.length = 0, hide.ns = T, step.increase = 0, label.x.npc = 'center', ) +
  rremove('xylab') +
  font("xy.text", size = 10)


plot_prop_selected_clusters_summary_enriched_in_pb_f = ggarrange(plot_prop_selected_clusters_summary_enriched_in_pb_f_b27, 
                                                                 plot_prop_selected_clusters_summary_enriched_in_pb_f_b38,
                                                                 ncol = 1, nrow = 2, align = 'h')

plot_prop_selected_clusters_summary_enriched_in_pb_f = annotate_figure(plot_prop_selected_clusters_summary_enriched_in_pb_f, 
                                                                       left = 'Cumulative proportion in total PB T cells', bottom = '')


# 9. Presence of clusters in PD1+ and CD137+ T-cells ----
enriched_clusters_in_pd1pos = shared.repertoire(pd1[grepl('pos', names(pd1))]) %>%
  inner_join(enriched_clusters_cumulative, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene'))
enriched_clusters_in_cd137pos = shared.repertoire(cd137) %>%
  inner_join(enriched_clusters_cumulative, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene'))

selected_in_pd1pos = shared.repertoire(pd1[grepl('pos', names(pd1))]) %>%
  inner_join(selected, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene')) %>%
  relocate(cluster, .after=V.gene)

selected_in_cd137pos = shared.repertoire(cd137) %>%
  inner_join(selected, by = c('CDR3.amino.acid.sequence', 'V.gene'='bestVGene')) %>%
  relocate(cluster, .after=V.gene)



# Fig 3A SeqLogo ----
seqlogo_clst701 = ggplot() + 
  ggseqlogo::geom_logo(sh_sf8beta_alice %>% 
                         dplyr::filter(cluster==701) %>% 
                         pull(CDR3.amino.acid.sequence), 
                       method = 'probability') + theme_classic()

ggarrange(ggarrange(seqlogo_clst1, seqlogo_clst10, seqlogo_clst12, seqlogo_clst404, ncol = 1, nrow = 4, labels=c(1,10,12,404), vjust = -1, common.legend = T),
          ggarrange(seqlogo_clst4, seqlogo_clst6, seqlogo_clst15, seqlogo_clst16,seqlogo_clst701, ncol=1, nrow=5, labels = c(4,6,15,16,701), vjust = -1,legend = 'none'), 
          ncol=2, 
          common.legend = T)


# Supp Fig 1 Proportion of PD-1+ and CD137+ cells from PB and SF T lymphocytes ----
plot_cd137 = ggpaired(facs %>% filter(fraction=='CD137'), id = 'donor', x='source', y='proportion', color='source', 
                      palette='npg', line.color = 'gray80') +
  stat_summary(fun.data = give.n, geom = "text", size=3, position = position_nudge(y=-0.01)) +
  stat_compare_means(hide.ns = T, label = 'p.signif', label.x.npc = 'middle') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + theme(legend.title = element_blank()) +
  rremove('xlab') + ylab('Proportion of CD137+ cells \n from T lymphocytes') 

plot_pd1 = ggpaired(facs %>% filter(fraction=='PD-1'), id = 'donor', x='source', y='proportion', color='source', 
                    palette='npg', line.color = 'gray80') +
  stat_summary(fun.data = give.n, geom = "text", size=3, position = position_nudge(y=-0.01)) +
  stat_compare_means(hide.ns = T, label = 'p.signif', label.x.npc = 'middle') +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) + theme(legend.title = element_blank()) + 
  rremove('xlab') + ylab('Proportion of PD-1+ cells \n from T lymphocytes')

ggarrange(plot_pd1, plot_cd137, nrow = 1, align = 'v', common.legend = T, legend = 'right')

#Supp Fig 2 Proportion of EBV-specific clonotyes in PB and SF of SpA patients ----

plot_vdjdb_f = ggpaired(summary_vdjdb_in_data_total_hla %>% filter(antigen.species=='EBV'), 
                        x='source', y='cumulative_prop', color='source',
                        id = 'donor', line.color = 'gray70',
                        outlier.shape = NA, add='point', add.params = list(alpha=0.3, position=position_quasirandom()), 
                        palette='npg') + #title = 'Total T-cell repertoire'
  stat_pvalue_manual(stat.test_vdjdb_hla, tip.length = 0, hide.ns = T, 
                     label = "p.adj.signif", label.x.npc = 'middle') +
  scale_y_log10(expand = expansion(mult = c(0.1, 0.2)),
                labels = trans_format(trans = 'log10',format = math_format())) +
  rremove('xlab') + ylab('Cumulative proportion of EBV-specific\n clonotypes from total T cells') +
  theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank()) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = log(1e-07,10)))

plot_vdjdb_cd8 = ggpaired(summary_vdjdb_in_cd8_hla %>% filter(antigen.species=='EBV'), 
                          x='source', y='cumulative_prop', color='source',
                          id = 'donor', line.color = 'gray70',
                          outlier.shape = NA, add='point', add.params = list(alpha=0.3, position=position_quasirandom()), 
                          palette='npg') + #title = 'CD8+ T-cell repertoire'
  stat_pvalue_manual(stat.test_vdjdb_cd8_hla, tip.length = 0, hide.ns = T, 
                     label = "p.signif") +
  scale_y_log10(expand = expansion(mult = c(0.1, 0.2)),
                labels = trans_format(trans = 'log10',format = math_format())) +
  rremove('xlab') + ylab('Cumulative proportion of EBV-specific\n clonotypes from CD8+ T cells') +
  theme(plot.title = element_text(hjust = 0.5)) +
  stat_summary(fun.data = give.n, geom = "text",size=3, position = position_nudge(y = log(1e-07,10)))


ggarrange(plot_vdjdb_f, plot_vdjdb_cd8, nrow = 1, labels = 'AUTO', common.legend = T, legend = 'right')

# Supp Fig 3 Enrichment of clusters in CD8+ SF vs CD8+ PB: Paired test  ----
prop_selected_clusters_summary_enriched_pb8_vs_sf8 = prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>% 
  mutate(cdr3_donor = paste(CDR3.amino.acid.sequence, donor)) %>% 
  mutate(proportion = if_else(is.na(proportion), 0, proportion))

stat.test_paored_b27_clonotypes = prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>%
  filter(donor %in% pb8_vs_sf8_paired_b27$donor) %>%
  filter(cluster %in% c(1,10,12,404)) %>%
  group_by(cluster) %>%
  wilcox_test(proportion~source, paired = T) %>%
  mutate(p.adj = p.adjust(p, 'fdr')) %>% 
  add_significance('p.adj') %>%
  mutate(y.position = 0.01)

zero_b27 = prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>%
  filter(donor %in% pb8_vs_sf8_paired_b27$donor) %>%
  filter(cluster %in% c(1,10,12,404)) %>%
  filter(proportion==0) %>% 
  count(cdr3_donor) %>% 
  filter(n>1) %>% 
  pull(cdr3_donor)

plot_paired_b27_clonotypes = ggpaired(prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>%
                                        filter(donor %in% pb8_vs_sf8_paired_b27$donor) %>%
                                        filter(cluster %in% c(1,10,12,404)) %>%
                                        filter(!(cdr3_donor %in% zero_b27)) %>%
                                        mutate(proportion = proportion + 1e-8),id = 'cdr3_donor', 
                                      x='source', y='proportion',
                                      line.size = 0.3, line.color = 'gray70',
                                      color='cluster', 
                                      # ylab='Cumulative proportion',
                                      # title ='B38+ donors' , 
                                      add = 'point', outlier.shape = NA,
                                      add.params = list(alpha=0.3, size=1,position = position_jitter(0.2))) + 
  facet_grid(~cluster) +
  scale_color_d3(drop=F,  guide=F) +
  scale_y_log10(expand = expansion(mult = c(0.01, 0.2)), 
                labels = trans_format(trans = 'log10',format = math_format())) + 
  rremove('xylab') +
  stat_pvalue_manual(stat.test_paored_b27_clonotypes, label = 'p.adj.signif', tip.length = 0, hide.ns = T, step.increase = 0)

stat.test_paored_b38_clonotypes = prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>%
  filter(donor %in% pb8_vs_sf8_paired_b38$donor) %>%
  filter(cluster %in% c(4,6,15,16,701)) %>%
  group_by(cluster) %>%
  wilcox_test(proportion~source, paired = T) %>%
  mutate(p.adj = p.adjust(p, 'fdr')) %>% 
  add_significance('p.adj') %>%
  mutate(y.position = 0.01)

zero_b38 = prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>%
  filter(donor %in% pb8_vs_sf8_paired_b38$donor) %>%
  filter(cluster %in% c(4,6,15,16,701)) %>%
  filter(proportion==0) %>% 
  count(cdr3_donor) %>% 
  filter(n>1) %>% 
  pull(cdr3_donor)

plot_paired_b38_clonotypes = ggpaired(prop_selected_clusters_summary_enriched_pb8_vs_sf8 %>%
                                        filter(donor %in% pb8_vs_sf8_paired_b38$donor) %>%
                                        filter(cluster %in% c(4,6,15,16,701)) %>%
                                        filter(!(cdr3_donor %in% zero_b38)) %>%
                                        mutate(proportion = proportion + 1e-8),
                                      id = 'cdr3_donor', 
                                      x='source', y='proportion',
                                      line.size = 0.3, line.color = 'gray70',
                                      color='cluster', 
                                      # ylab='Cumulative proportion',
                                      # title ='B38+ donors' , 
                                      add = 'point', outlier.shape = NA,
                                      add.params = list(alpha=0.3, size=1,position = position_jitter(0.2))) + 
  facet_grid(~cluster) +
  scale_color_d3(drop=F,  guide=F) +
  scale_y_log10(expand = expansion(mult = c(0.01, 0.2)), 
                labels = trans_format(trans = 'log10',format = math_format())) + 
  rremove('xylab') +
  stat_pvalue_manual(stat.test_paored_b38_clonotypes, label = 'p.adj.signif', tip.length = 0, hide.ns = T, step.increase = 0)

annotate_figure(ggarrange(plot_paired_b27_clonotypes,
                          plot_paired_b38_clonotypes,
                          ncol = 1, nrow = 2, align = 'h'), left = 'Proportion from\n CD8+ T cells of SpA patients')


#Supp Fig 4 ALICE on TCRalpha ----
sf8alpha = get.inframes(spa_sf_8_a)
sf8alpha = lapply(sf8alpha, function(x) {x=x[x$Read.count>1,];x})
sf8alpha = lapply(sf8alpha, function(x) {x=as.data.table(x);x})
sf8alpha = lapply(sf8alpha, function(x) {names(x)[grep('V.gene$', names(x))]='bestVGene';x})
sf8alpha = lapply(sf8alpha, function(x) {names(x)[grep('J.gene$', names(x))]='bestJGene';x})

sf8alpha_alice<-ALICE_pipeline_alpha(DTlist=sf8alpha,folder="sf_cd8_a/",cores=7,iter=25,nrec=5e6, 
                                            Read_count_filter = 0)
sf8alpha_alice = lapply(sf8alpha_alice, function(x) x= as.data.frame(x))
sh_sf8alpha_alice = shared.repertoire(sf8alpha_alice, .by.col = c('CDR3.amino.acid.sequence', 'bestVGene', 'bestJGene'), .sum.col = 'Read.count')


q = replicate(length(table(sh_sf8alpha_alice$bestVGene)), paste(sample(LETTERS, 5, replace=TRUE), collapse=""))
min(stringdist::stringdistmatrix(q, method = 'hamming'))
vgenes = as.data.frame(table(sh_sf8alpha_alice$bestVGene), stringsAsFactors = F)
vgenes$chifr = q
for (i in 1:nrow(vgenes)){
  sh_sf8alpha_alice[grep(paste(vgenes[i,1],'$', sep = ''), sh_sf8alpha_alice$bestVGene),1] = paste(sh_sf8alpha_alice[grep(paste(vgenes[i,1],'$', sep = ''), sh_sf8alpha_alice$bestVGene),1],vgenes[i,3], sep = '')
}

gr_sf8alpha_alice = mutation.network(sh_sf8alpha_alice, .method = 'hamm', .max.errors = 1, 
                                     .label.col = 'CDR3.amino.acid.sequence', .seg.col = 'bestVGene')
gr_sf8alpha_alice = set.vertex.attribute(gr_sf8alpha_alice, name = 'label2', index = V(gr_sf8alpha_alice), 
                                         value = substr(sh_sf8alpha_alice$CDR3.amino.acid.sequence, start = 0, stop = nchar(sh_sf8alpha_alice$CDR3.amino.acid.sequence)-5))
