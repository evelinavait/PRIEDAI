setwd("C:/Users/win7/Desktop/PRIEDAI")
library(tidyverse)
library(DESeq2)
library(isomiRs)
library(outliers)
library(plotly)

# metadata
meta.se <- openxlsx::read.xlsx("gc_ap_hc_metafile.xlsx")

# data.frame
rownames(meta.se) <- meta.se$Sample_ID

# counts
ids.se <- openxlsx::read.xlsx("mirtop_rawData.xlsx")

# read count generation from mirtop output
ids_se <- IsomirDataSeqFromMirtop(ids.se, meta.se)

# saveRDS (.rda for saving multiple objects into a single file)
saveRDS(ids_se, file = "ids_se.rda")
ids_se <- readRDS(file = "ids_se.rda")

# summary of traits
colData(ids_se) %>% as.data.frame() %>% select(Condition) %>% group_by(Condition) %>% summarise(n=n())

# check traits on MDS
mds <- isoNorm(ids_se) %>% counts(., norm=TRUE) %>% t() %>% dist() %>% cmdscale() %>% as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% dplyr::rename("MDS1" = "V1", "MDS2" = "V2") %>%
  left_join(., colData(ids_se) %>% as.data.frame(), by = "Sample_ID")
p1 <- ggplot(mds,aes(MDS1,MDS2,text=paste("Sample_ID:", Sample_ID))) +
  geom_point(aes(color = Condition), data = mds) + stat_ellipse(aes(MDS1, MDS2, color = Condition), data = mds, inherit.aes = FALSE)
ggplotly(p1)
p2 <- p1 + scale_color_brewer(palette = "Dark2",name = "Tiriamuju grupe", label = c("AG", "SV-nav", "SV-?alia", "K"))
ggsave(p2, filename = "p2.jpg", width = 6, height = 4, dpi = 600)

# sequencing depth
lib.size <- counts(ids_se) %>% colSums() %>% as.data.frame() %>% dplyr::rename("mapped_to_mirnas" = ".") %>%
  rownames_to_column("Sample_ID") %>% left_join(colData(ids_se) %>% as.data.frame() %>%
                                                  select(Sample_ID, Condition), ., by = "Sample_ID")
# number of detected miRNAs
mirna.nr <- counts(ids_se) %>% as.data.frame() %>% mutate_all(~if_else(.x>0, 1, 0)) %>% colSums() %>%
  as.data.frame() %>% dplyr::rename("unique_mirnas" = ".") %>% rownames_to_column("Sample_ID")

# quality control tag
qc.tb <- left_join(lib.size, mirna.nr, by = "Sample_ID") %>%
  mutate(lib_outlier = if_else(mapped_to_mirnas < 1e6, "FAILED", "PASSED")) %>%
  mutate(mir_outlier = scores(log2(unique_mirnas), type = "iqr", lim = 1.5)) %>%
  mutate(mir_outlier = if_else(mir_outlier == TRUE, "FAILED", "PASSED")) %>%
  mutate(QC = if_else(lib_outlier == "PASSED" & mir_outlier == "PASSED", "PASSED", "FAILED"))

# check outliers on MDS
mds <- isoNorm(ids_se) %>% counts(., norm=TRUE) %>% t() %>% dist() %>% cmdscale() %>% as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% dplyr::rename("MDS1" = "V1", "MDS2" = "V2") %>%
  left_join(., qc.tb, by = "Sample_ID")

p3 <- ggplot(mds, aes(MDS1, MDS2,text=paste("Sample_ID:", Sample_ID, "\n", "Condition:", Condition))) +
  geom_point(aes(color = QC), data = mds) + stat_ellipse(aes(MDS1, MDS2, color = QC), data = mds, inherit.aes = FALSE)
ggplotly(p2)
p4 <- p3 + scale_color_brewer(palette = "Dark2", direction = -1, name = "Kokyb?s kontrol?", label = c("Nuo vidurkio nutolusi reik?m?", "Nuo vidurkio nenutolusi reik?m?"))
ggsave(p4, filename = "p4.jpg", width = 6, height = 3, dpi = 600)

# remove outliers (qced - qc editor)
ids_se_qced <- ids_se[, qc.tb %>% filter(QC %in% "PASSED") %>% .$Sample_ID]

# generate DESeq2 object
dds_se <- DESeqDataSetFromMatrix(countData = counts(ids_se_qced), colData = colData(ids_se_qced), design = ~Condition)
dds_se <- DESeq(dds_se)

# save RDS
saveRDS(dds_se, file = "dds_se.rda")
dds_se <- readRDS(file = "dds_se.rda")

# calculate variance and mean for every miRNA in raw and vst counts
vars.tb <- counts(dds_se, norm = FALSE) %>% as.data.frame() %>% rownames_to_column("mir") %>% gather(Sample_ID, value, -mir) %>% group_by(mir) %>%
  summarise(var_raw = var(value), mean_raw = mean(value), coef_var_raw = sd(value)/mean(value)) %>%
  left_join(., getVarianceStabilizedData(dds_se) %>% as.data.frame() %>% rownames_to_column("mir") %>%
              gather(Sample_ID, value, -mir) %>% group_by(mir) %>%
              summarise(var_vst = var(value), mean_vst = mean(value), coef_var_vst = sd(value)/mean(value)), by = "mir")

# mean-variance realationship raw data
p5 <- ggplot(vars.tb, aes(mean_raw, var_raw, text=paste(mir))) +
  geom_point() + geom_smooth(data = vars.tb, aes(mean_raw, var_raw), inherit.aes = FALSE, method = "lm") +
  scale_x_log10() + scale_y_log10() +
  labs(x="Vidurkis", y = "Variacija")
ggplotly(p5)
ggsave(p5, filename = "p5.jpg", width = 5, height = 3, dpi = 600)

# mean-variance realationship vst data
p6 <- ggplot(vars.tb, aes(mean_vst, var_vst, text=paste(mir))) +
  geom_point() + geom_smooth(data = vars.tb, aes(mean_vst, var_vst), inherit.aes = FALSE, method = "lm") +
  scale_x_log10() +
  labs(x="Vidurkis", y = "Variacija")
ggplotly(p6)
ggsave(p6, filename = "p6.jpg", width = 5, height = 3, dpi = 600)

# mean vs coef of var realationship raw data
p7 <- ggplot(vars.tb, aes(mean_raw, coef_var_raw, text=paste(mir))) +
  geom_point() + geom_smooth(data = vars.tb, aes(mean_raw, coef_var_raw), inherit.aes = FALSE) +
  geom_vline(xintercept = 1, color = "red") +
  scale_x_log10() +
  labs(x="Vidurkis", y = "Variacijos koeficientas")
ggplotly(p7)
ggsave(p7, filename = "p7.jpg", width = 5, height = 3, dpi = 600)

# mean vs coef of var realationship vst data
p8 <- ggplot(vars.tb, aes(mean_vst, coef_var_vst, text=paste(mir))) +
  geom_point() + geom_smooth(data = vars.tb, aes(mean_vst, coef_var_vst), inherit.aes = FALSE) +
  geom_vline(xintercept = 1.7, color = "red") +
  scale_x_log10() +
  labs(x="Vidurkis", y = "Variacijos koeficientas")
ggplotly(p8)
ggsave(p8, filename = "p8.jpg", width = 5, height = 3, dpi = 600)

# filter low-abundant and non-variable miRNAs
qced_mirs <- vars.tb %>% filter(mean_raw > 1) %>% .$mir

# re-run normalization
dds_se_qced <- DESeqDataSetFromMatrix(countData = counts(dds_se)[qced_mirs,], colData = colData(dds_se), design = ~Condition)
dds_se_qced <- DESeq(dds_se_qced)

# saveRDS   
saveRDS(dds_se_qced, file = "dds_se_qced.rda")
dds_se_qced <- readRDS(file = "dds_se_qced.rda")

# re-calculate variance and mean for every miRNA in vst counts
vars.tb.qced <- getVarianceStabilizedData(dds_se_qced) %>% as.data.frame() %>% rownames_to_column("mir") %>%
  gather(Sample_ID, value, -mir) %>% group_by(mir) %>%
  summarise(var_vst = var(value), mean_vst = mean(value), coef_var_vst = sd(value)/mean(value))

# qced miRNA final data set
p9 <- ggplot(vars.tb.qced, aes(mean_vst, coef_var_vst, text=paste(mir))) +
  geom_point() + geom_smooth(data = vars.tb.qced, aes(mean_vst, coef_var_vst), inherit.aes = FALSE) +
  scale_x_log10() +
  labs(x="Vidurkis", y = "Variacijos koeficientas")
ggplotly(p9)
ggsave(p9, filename = "p9.jpg", width = 5, height = 3, dpi = 600)

# MDS of final miRNA-arm dataset
p10 <- getVarianceStabilizedData(dds_se_qced) %>% t() %>% dist() %>% cmdscale() %>% as.data.frame() %>%
  rownames_to_column("Sample_ID") %>% dplyr::rename("MDS1" = "V1", "MDS2" = "V2") %>%
  left_join(., colData(ids_se_qced) %>% as.data.frame(), by = "Sample_ID") %>%
  ggplot(., aes(MDS1, MDS2, color = Condition)) + geom_point() + 
  stat_ellipse()
ggplotly(p10)

p11 <- p10 + scale_colour_brewer(palette = "Dark2", name = "Tiriamuju grupe", labels = c("AG", "SV-nav", "SV-?alia", "K"))
ggsave(p11, filename = "p11.jpg", width = 6, height = 4, dpi = 600)
summarise(dds_se_qced)

# extract the results from model (normal lfc shrinkage)
res.gc.ag <- results(dds_se_qced,contrast=c("Condition", "GC", "AG")) %>%
  as.data.frame() %>% rownames_to_column("mir") %>% mutate(comparison="GC_vs_AG")
res.gc.gcaj <- results(dds_se_qced,contrast=c("Condition", "GC", "GCaj")) %>%
  as.data.frame() %>% rownames_to_column("mir") %>% mutate(comparison="GC_vs_GCaj")
res.gc.hc <- results(dds_se_qced,contrast=c("Condition", "GC", "HC")) %>%
  as.data.frame() %>% rownames_to_column("mir") %>% mutate(comparison="GC_vs_HC")
res.ag.gcaj <- results(dds_se_qced,contrast=c("Condition", "AG", "GCaj")) %>%
  as.data.frame() %>% rownames_to_column("mir") %>% mutate(comparison="AG_vs_GCaj")
res.ag.hc <- results(dds_se_qced,contrast=c("Condition", "AG", "HC")) %>%
  as.data.frame() %>% rownames_to_column("mir") %>% mutate(comparison="AG_vs_HC")
res.gcaj.hc <- results(dds_se_qced,contrast=c("Condition", "GCaj", "HC")) %>%
  as.data.frame() %>% rownames_to_column("mir") %>% mutate(comparison="GCaj_vs_HC")

# threshold of dea: FDR adjusted p-value < 0.05
res.mir <- bind_rows(res.ag.gcaj, res.ag.hc, res.gcaj.hc, res.gc.ag, res.gc.gcaj, res.gc.hc) %>%
  mutate(dea = if_else(padj<0.05, "de", "non-de", "non-de")) %>%
  mutate(comparison = factor(comparison, levels = unique(comparison)))

p12 <- res.mir %>% ggplot(aes(baseMean, log2FoldChange, text = paste0("mirna: ", mir, "\n", "padj: ", padj))) +
  geom_point(aes(color = dea, alpha = dea)) + geom_hline(yintercept = 0, linetype=2, color = "red") +
  facet_wrap(~comparison) + scale_x_log10() + scale_alpha_discrete(range = c(0.6, 0.2)) +
  
  scale_y_continuous(limits = c(-3, 3), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  labs(x= "Vidutin? rai?ka", y= "Pokytis kartais (log2PK)")
ggplotly(p12)
p13 <- p12 + scale_colour_brewer(palette = "Dark2", name = "Rai?kos pokytis", labels = c("Pakitusi rai?ka", "Nepakitusi rai?ka")) +
  scale_fill_binned(name = "Rai?kos pokytis", labels = c("Pakitusi rai?ka", "Nepakitusi rai?ka"))
p14 <- ggsave(p13, filename = "p12.jpg", width = 6, height = 4, dpi = 600)

#Validavimo studija
library(readxl)
library(xlsx)
library(dplyr)
library(ggplot2)

#1. isikelti exel fila
con.vs.gc <- read_xlsx("C:/Users/win7/Desktop/PRIEDAI/results_miRNA_GCvsHC.xlsx", sheet=1, col_names = F)[-(1:2),]
colnames(con.vs.gc) <- con.vs.gc[1,]
con.vs.gc <- con.vs.gc[-1,]
con.vs.gc$case <- "CONvsGC"

con.vs.ag <- read_xlsx("C:/Users/win7/Desktop/PRIEDAI/results_miRNA_AGvsHC.xlsx", sheet=1, col_names = F)[-(1:2),]
colnames(con.vs.ag) <- con.vs.ag[1,]
con.vs.ag <- con.vs.ag[-1,]
con.vs.ag$case <- "CONvsAG"

ag.vs.gc <- read_xlsx("C:/Users/win7/Desktop/PRIEDAI/results_miRNA_GCvsAG.xlsx", sheet=1, col_names = F)[-(1:2),]
colnames(ag.vs.gc) <- ag.vs.gc[1,]
ag.vs.gc <- ag.vs.gc[-1,]
ag.vs.gc$case <- "AGvsGC"

#2. patikrinti ar vertes yra numeric
write.xlsx(con.vs.gc, file = "C:/Users/win7/Desktop/PRIEDAI/con.vs.gc.xlsx", col.names = T, row.names = T)
con.vs.gc <- read_xlsx("C:/Users/win7/Desktop/PRIEDAI/con.vs.gc.xlsx")
write.xlsx(con.vs.ag, file = "C:/Users/win7/Desktop/PRIEDAI/con.vs.ag.xlsx", col.names = T, row.names = T)
con.vs.ag <- read_xlsx("C:/Users/win7/Desktop/PRIEDAI/con.vs.ag.xlsx")
write.xlsx(ag.vs.gc, file = "C:/Users/win7/Desktop/PRIEDAI/ag.vs.gc.xlsx", col.names = T, row.names = T)
ag.vs.gc <- read_xlsx("C:/Users/win7/Desktop/PRIEDAI/ag.vs.gc.xlsx")

str(ag.vs.gc)
con.vs.gc$log2FoldChange <- as.numeric(con.vs.gc$log2FoldChange)
con.vs.gc$padj <- as.numeric(con.vs.gc$padj)

con.vs.ag$log2FoldChange <- as.numeric(con.vs.ag$log2FoldChange)
con.vs.ag$padj <- as.numeric(con.vs.ag$padj)

ag.vs.gc$log2FoldChange <- as.numeric(ag.vs.gc$log2FoldChange)
ag.vs.gc$padj <- as.numeric(ag.vs.gc$padj)

#3. duomenu filtravimas: log2FC < -0 arba > 0; padj < 0.05
con.vs.gc.filtered <- filter(con.vs.gc, log2FoldChange < -0 | log2FoldChange > 0) %>% filter(padj < 0.05)
con.vs.ag.filtered <- filter(con.vs.ag, log2FoldChange < -0 | log2FoldChange > 0) %>% filter(padj < 0.05)
ag.vs.gc.filtered <- filter(ag.vs.gc, log2FoldChange < -0 | log2FoldChange > 0) %>% filter(padj < 0.05)


#4. save xlsx
mir.seq.filtered.data.all.de <- rbind(con.vs.gc.filtered, con.vs.ag.filtered, ag.vs.gc.filtered)
mir.seq.filtered.data.all.de$absolute_FC <- abs(mir.seq.filtered.data.all.de$log2FoldChange)
#mir.seq.filtered.data.2 <- filter(mir.seq.filtered.data, ...2 %in% mircount2)
write.xlsx(mir.seq.filtered.data.all.de, "C:/Users/win7/Desktop/PRIEDAI/mir.seq.filtered.data.all.der.20201216.xlsx") # failas su visomis statistiskai 
#reiksmingai pakitusiomis mir exelyje paimti top 10 
# dereguliuotu (pagal absolute FC)

mircount3 <- as.data.frame(table(mir.seq.filtered.data.all.de$...2)) %>% filter(Freq == 3) #miR-196a-5p pakitusi con vs ag; con vs. gc; cag vs. gc
mircount2 <- as.data.frame(table(mir.seq.filtered.data.all.de$...2)) %>% filter(Freq == 2)

##miR-196a-5p paveikslas###
mir196a <- filter(mir.seq.filtered.data.all.de, ...2 == "hsa-miR-196a-5p")
mir196a.FC <- ggplot(mir196a, aes(x=factor(case, levels=c("CONvsAG", "AGvsGC", "CONvsGC")), y=log2FoldChange)) + geom_col() + xlab("Grupė") + labs(title = "hsa-miR-196a-5p")
ggsave("C:/Users/win7/Desktop/PRIEDAI/mir196a.tiff", plot = last_plot(), dpi = 300) # vienintele mikro RNR kuri pakitusi visais trim atvejais 
