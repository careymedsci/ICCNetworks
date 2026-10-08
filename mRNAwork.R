
# ============================================================================
# Copyright (c) 2025 Thomas Broggini and Liu Xiao. All rights reserved.
#
# This code was jointly written by:
# Thomas Broggini (Frankfurt University, Germany)
# Liu Xiao (Xiangyang First Hospital, China)
#
# The way how the mathematical models and methodologies implemented in this code have been
# individually customized and are not intended for generic use. 
#
# For permission requests or inquiries, please contact:
# Thomas Broggini  : broggini@med.uni-frankfurt.de
# Liu Xiao         : careyneurosurgery@gmail.com  /  carey-lau@foxmail.com
#
# Unauthorized use will be considered a violation of intellectual property rights.
# ============================================================================



# 作者：Thomas Broggini 和 刘晓；法兰克福大学医院；湖北医药学院襄阳市第一人民医院
# By Thomas Broggini (Supervisor) 
# Xiao Liu (Student)
# Department of Neurosurgery
# Neuroscience Centre
# University Hospital Frankfurt
# Goethe University Frankfurt, Germany
# Frankfurt Cancer Institute, Germany

# Xiangyang No.1 people's Hospital, China
# Also an python version coded 
# xiao.liu@stud.uni-frankfurt.de
# 2025.03.09

########################################################### first part DEG 
library(pathview)
library(org.Mm.eg.db)
library(Cairo)
library(DESeq2)


data <- read.csv("raw_counts.csv", header = TRUE, row.names = "ensembl_id")
data <- data[, sort(colnames(data))]

F10 <- c("sample1", "sample2", "sample3")  
BrM3 <- c("sample4", "sample5", "sample6")
data <- data[, c(F10, BrM3)]
Group <- factor(
  c(rep("B16-F10", length(F10)), 
    rep("B16-BrM3", length(BrM3))),
  levels = c("B16-F10", "B16-BrM3")
)

my_colData <- data.frame(Group, row.names = colnames(data))

library(edgeR)

cpm_counts <- cpm(data)  
keep <- rowSums(cpm_counts > 1) >= 0.75 * ncol(data)
counts_filtered <- data[keep, ]


# create DESeq2 obj -----------------------------------------------------------
dds <- DESeqDataSetFromMatrix(
  countData = data,
  colData = my_colData,
  design = ~ Group
)
# DESeq2,has the ability handle with 0 counts samples---------
dds <- DESeq(dds)

# method 1
# res <- results(
#  dds,
#  contrast = c("Group", "B16-BrM3", "B16-F10"),
#  alpha = 0.01
#)

#res$padj <- p.adjust(res$pvalue, method = "fdr")
#res$padj[res$padj == 0] <- 1e-300

# method 2 -------------------------------------------------------------
res_shrink <- lfcShrink(
  dds,
  contrast = c("Group", "B16-BrM3", "B16-F10"),
  type = "ashr" # for visualization
)

# result for dds -------------------------------------------------------------
library(dplyr)
library(tibble)

res_df <- as.data.frame(res_shrink) %>%
  rownames_to_column("gene_name") %>%
  mutate(
    Direction = case_when(
      log2FoldChange > 1.5 & padj < 0.05 ~ "Upregulated",   
      log2FoldChange < -1.5 & padj < 0.05 ~ "Downregulated", 
      TRUE ~ "Not significant"  #  "Not significant"
    )
  )

# 
head(res_df)



############################################################## second part
library(rtracklayer)
library(dplyr)


gtf <- import("mm39en.gtf")
gtf_meta <- as.data.frame(mcols(gtf))

colnames(gtf_meta)

gene_annot <- gtf_meta %>%
  filter(type == "gene") %>%
  distinct(gene_id, gene_name)

# add gene id for accurate caculation
res_annotated <- res_df %>%
  left_join(gene_annot, by = "gene_name")


############################################################### third part
library(ggplot2)
library(ggrepel)
library(KEGGREST)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(GO.db)
library(dplyr)

# ------------------- 1) set cutoff -------------------
log2FC_cutoff <- 1.5
padj_cutoff <- 0.05
res_annotated$padj[res_annotated$padj == 0] <- 1e-300

# ------------------- 2) KEGG Calcium -------------------
kegg_raw <- tryCatch(keggGet("mmu04020")[[1]]$GENE, error = function(e) NULL)
kegg_entrez <- character(0)
if (!is.null(kegg_raw)) {
  kegg_entrez <- kegg_raw[seq(1, length(kegg_raw), by = 2)]
  kegg_entrez <- gsub(" ", "", kegg_entrez)
}
kegg_map <- character(0)
if (length(kegg_entrez) > 0) {
  kegg_map_df <- AnnotationDbi::select(org.Mm.eg.db,
                                       keys = kegg_entrez,
                                       keytype = "ENTREZID",
                                       columns = c("SYMBOL"))
  kegg_map <- unique(na.omit(kegg_map_df$SYMBOL))
}

# ------------------- 3) GO calcium -------------------
go_ids <- keys(GO.db, keytype = "GOID")
go_terms <- sapply(go_ids, function(g) {
  tryCatch(GOTERM[[g]]@Term, error = function(e) NA)
}, USE.NAMES = TRUE)
has_calcium <- names(go_terms)[grepl("calcium", go_terms, ignore.case = TRUE)]
has_calcium <- unique(na.omit(has_calcium))

go_genes_df <- AnnotationDbi::select(org.Mm.eg.db,
                                     keys = has_calcium,
                                     keytype = "GO",
                                     columns = c("SYMBOL","GO"))
go_genes <- unique(na.omit(go_genes_df$SYMBOL))

# ------------------- 4) GPCR genes -------------------
all_symbols <- keys(org.Mm.eg.db, keytype = "SYMBOL")
gene_go <- AnnotationDbi::select(org.Mm.eg.db,
                                 keys = all_symbols,
                                 columns = c("SYMBOL","GO"),
                                 keytype = "SYMBOL")
gpcr_go_terms <- c("GO:0004930","GO:0007186")
gpcr_genes <- gene_go %>%
  filter(GO %in% gpcr_go_terms) %>%
  distinct(SYMBOL)

# ------------------- 5) all related genes -------------------
highlight_genes <- unique(c(kegg_map, go_genes, gpcr_genes$SYMBOL))

# ------------------- 6) prepare data for plotting -------------------
res_plot <- res_annotated %>%
  mutate(
    sig = (padj < padj_cutoff) & (abs(log2FoldChange) > log2FC_cutoff),
    need_highlight = sig & (gene_name %in% highlight_genes),
    label_color = case_when(
      need_highlight & log2FoldChange > 0 ~ "red",
      need_highlight & log2FoldChange < 0 ~ "blue",
      TRUE ~ "bg"
    ),
    r = sqrt(log2FoldChange^2 + (-log10(padj))^2),
    alpha = case_when(
      label_color == "bg" ~ exp(-r * 0.6),
      TRUE ~ 1
    ),
    point_size = case_when(
      label_color == "bg" ~ 3,
      TRUE ~ 3
    )
  )

# ------------------- 7) plotting -------------------
volcano_plot <- ggplot(res_plot, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(data = filter(res_plot, label_color == "bg"),
             aes(alpha = alpha, size = point_size),
             color = "black") +
  geom_point(data = filter(res_plot, label_color == "red"),
             aes(size = point_size),
             color = "red") +
  geom_point(data = filter(res_plot, label_color == "blue"),
             aes(size = point_size),
             color = "blue") +
  geom_text_repel(data = filter(res_plot, need_highlight),
                  aes(label = gene_name),
                  size = 4,
                  max.overlaps = 20) +
  geom_vline(xintercept = c(-log2FC_cutoff, log2FC_cutoff),
             linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_cutoff),
             linetype = "dashed") +
  theme_bw() +
  theme(legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  xlab("log2 Fold Change") +
  ylab("-log10 Adjusted p-value") +
  ggtitle("Volcano Plot (Calcium + GPCR Highlight)")

print(volcano_plot)

# ------------------- 8) save plotting by SVG -------------------
if (!requireNamespace("svglite", quietly = TRUE)) install.packages("svglite")
ggsave("volcano_plot.svg", plot = volcano_plot, device = "svg", width = 8, height = 6, dpi = 1200)

# ------------------- 9) save gene list to Excel -------------------
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)

highlight_df <- res_plot %>%
  filter(need_highlight) %>%
  select(gene_name, log2FoldChange, padj)

# sorting
highlight_df <- highlight_df %>%
  mutate(
    source = case_when(
      gene_name %in% kegg_map ~ "KEGG_calcium",
      gene_name %in% go_genes ~ "GO_calcium",
      gene_name %in% gpcr_genes$SYMBOL ~ "GPCR",
      TRUE ~ "Other"
    )
  )

#  Excel
wb <- createWorkbook()

sources <- unique(highlight_df$source)
for (s in sources) {
  addWorksheet(wb, s)
  writeData(wb, sheet = s, highlight_df %>% filter(source == s))
}

saveWorkbook(wb, file = "highlight_genes.xlsx", overwrite = TRUE)



########################################################

# KEGG Enrichment Analysis   fourth part

library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(ggplot2)
library(DOSE)
library(dplyr)

# ------------------- 1) Filter significant DEGs -------------------

sig_genes_kegg <- res_annotated %>%
  filter(Direction != "Not significant") %>%
  filter(!is.na(gene_id))

# ------------------- 2) Map ENSEMBL -> ENTREZID -------------------

gene_map <- bitr(
  unique(sig_genes_kegg$gene_id),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
) %>%
  distinct(ENSEMBL, .keep_all = TRUE)

sig_genes_kegg <- sig_genes_kegg %>%
  left_join(gene_map, by = c("gene_id" = "ENSEMBL")) %>%
  filter(!is.na(ENTREZID))

###############################################################

# ------------------- 3A) Try online KEGG enrichment -------------------

kegg_sig <- tryCatch({
  enrichKEGG(
    gene = sig_genes_kegg$ENTREZID,
    organism = "mmu",
    pvalueCutoff = 0.05
  )
}, error = function(e) {
  message("⚠️ enrichKEGG failed")
  return(NULL)
})


# ------------------- 3B) Offline KEGG enrichment (added) -------------------

if (is.null(kegg_sig)) {
  
  # ---- Load local KEGG files ----
  
  pathway_list <- read.table(
    "mmu_pathway_list.txt",
    sep = "\t",
    header = FALSE,
    col.names = c("pathway", "description"),
    stringsAsFactors = FALSE
  )
  
  pathway_link <- read.table(
    "mmu_pathway_link.txt",
    sep = "\t",
    header = FALSE,
    col.names = c("gene", "pathway"),
    stringsAsFactors = FALSE
  )
  
  # ---- Clean IDs ----
  
  pathway_list$pathway <- sub("path:", "", pathway_list$pathway)
  pathway_link$pathway <- sub("path:", "", pathway_link$pathway)
  pathway_link$gene <- sub("mmu:", "", pathway_link$gene)
  
  # ---- Construct TERM2GENE & TERM2NAME ----
  
  PATHID2NAME <- as.data.frame(pathway_list)
  PATHID2GENE <- pathway_link[, c("pathway", "gene")]
  
  # ---- Offline enrichment ----
  
  kegg_sig <- enricher(
    gene = sig_genes_kegg$ENTREZID,
    TERM2GENE = PATHID2GENE,
    TERM2NAME = PATHID2NAME,
    pvalueCutoff = 0.05
  )
  
  # ---- Convert ENTREZID → SYMBOL ----
  
  kegg_sig <- setReadable(kegg_sig, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
}


# ------------------- 4) Filter virus-related pathways -------------------


virus_pathways <- c(
  "Herpes simplex virus 1 infection",
  "Influenza A",
  "Human papillomavirus infection",
  "Epstein-Barr virus infection",
  "Hepatitis C"
)

kegg_df <- as.data.frame(kegg_sig)
kegg_df_filtered <- kegg_df %>%
  filter(!Description %in% virus_pathways)

kegg_sig_filtered <- kegg_sig
kegg_sig_filtered@result <- kegg_sig_filtered@result %>%
  filter(ID %in% kegg_df_filtered$ID)


# ------------------- 5) Dotplot -------------------

p <- dotplot(kegg_sig_filtered, showCategory = 25, title = "KEGG Enrichment: Significant DEGs") +
  scale_color_gradient2(
    low = rgb(0,113,188, maxColorValue = 255),
    mid = "white",
    high = rgb(163,30,50, maxColorValue = 255),
    midpoint = 0
  )

print(p)


# ------------------- 6) Save filtered KEGG enrichment results -------------------
write.csv(as.data.frame(kegg_sig_filtered), 
          file = "KEGG_sig_genes_results_filtered.csv", 
          row.names = FALSE)

cat("KEGG enrichment for significant DEGs completed (virus-related pathways removed)!\n")



######################################################## fifth part
# KEGG enrichment relationship visualization
library(enrichplot)
library(clusterProfiler)
library(igraph)
library(ggraph)

# ------------------- 1) caculation -------------------
#Jaccard similarity / gene overlap
kegg_sim <- pairwise_termsim(kegg_sig_filtered) # please first map all the genes

# ------------------- 2) emapplot -------------------

emap <- emapplot(
  kegg_sim,
  layout = "nicely",  # "kk", "fr", "nicely"
  showCategory = 100
)
print(emap)

# ------------------- 3) cnetplot: -------------------

cnet <- cnetplot(
  kegg_sig_filtered,
  categorySize = "pvalue",
  foldChange = NULL,
  showCategory = 10,
  node_label = "gene",
  layout = "kk"
)
print(cnet)

# ------------------- 4) save -------------------
# emap
ggsave("KEGG_enrichment_relationship_map.svg", plot = emap, device = "svg", width = 8, height = 6, dpi = 1200)

# cnet
ggsave("KEGG_enrichment_gene_network.svg", plot = cnet, device = "svg", width = 10, height = 8, dpi = 1200)

cat("KEGG enrichment relationship visualization completed!\n")



########################################################  sixth part
# GO enrichment (BP / CC / MF) dotplot for significant DEGs
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(ggplot2)
library(dplyr)

# ------------------- 1) Filter significant DEGs -------------------
sig_genes_go <- res_annotated %>%
  filter(Direction != "Not significant") %>%
  filter(!is.na(gene_id))  # remove NA

# ------------------- 2) Map ENSEMBL -> ENTREZID -------------------
gene_map_go <- bitr(
  unique(sig_genes_go$gene_id),
  fromType = "ENSEMBL",
  toType = "ENTREZID",
  OrgDb = org.Mm.eg.db
) %>% distinct(ENSEMBL, .keep_all = TRUE)

sig_genes_go <- sig_genes_go %>%
  left_join(
    gene_map_go,
    by = c("gene_id" = "ENSEMBL")
  ) %>%
  filter(!is.na(ENTREZID))

# ------------------- 3) GO enrichment for BP / CC / MF -------------------
go_ontologies <- c("BP", "CC", "MF")
go_results_list <- list()

for (ont in go_ontologies) {
  go_enrich <- enrichGO(
    gene = sig_genes_go$ENTREZID,
    OrgDb = org.Mm.eg.db,
    keyType = "ENTREZID",
    ont = ont,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  go_results_list[[ont]] <- go_enrich
}

# ------------------- 4) Dotplot for each ontology -------------------
for (ont in go_ontologies) {
  go_enrich <- go_results_list[[ont]]
  
  if (nrow(as.data.frame(go_enrich)) == 0) {
    cat(paste0("No significant GO terms found for ", ont, "\n"))
    next
  }
  
  p <- dotplot(go_enrich, showCategory = 20, title = paste0("GO Enrichment: ", ont)) +
    scale_color_gradient2(
      low = rgb(0,113,188, maxColorValue = 255),
      mid = "white",
      high = rgb(163,30,50, maxColorValue = 255),
      midpoint = 0
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  print(p)
  
  # Save plot as SVG
  ggsave(
    filename = paste0("GO_dotplot_", ont, ".svg"),
    plot = p,
    device = "svg",
    width = 8,
    height = 6,
    dpi = 1200
  )
  
  # Save GO enrichment table
  write.csv(
    as.data.frame(go_enrich),
    file = paste0("GO_enrichment_", ont, ".csv"),
    row.names = FALSE
  )
}

cat("GO enrichment dotplots for BP / CC / MF completed!\n")

writeLines(svg_lines, "your_image_fixed.svg")






########################################################  seventh part


library(pathview)
library(org.Mm.eg.db)

## 1. gene_fc
gene_fc <- sig_genes_all$log2FoldChange
names(gene_fc) <- sig_genes_all$ENTREZID

## 2. KEGG list
kegg_list <- c(
  "mmu04020","mmu04310","mmu04530","mmu00190","mmu01100",
  "mmu04080","mmu04110","mmu04151","mmu04330","mmu04360",
  "mmu04361","mmu04512","mmu04514","mmu04540","mmu04722",
  "mmu04721","mmu04724","mmu04725","mmu04726","mmu04727",
  "mmu04728","mmu04820","mmu04916","mmu05218"
)

## 3.
if (!dir.exists("pathview_pdf")) dir.create("pathview_pdf")
setwd("pathview_pdf")

## 4. 
for (pid in kegg_list) {
  pathview(
    gene.data   = gene_fc,
    pathway.id  = pid,
    species     = "mmu",
    kegg.native = FALSE,  
    out.suffix  = "fc",
    output      = "pdf",   
    same.layer  = TURE,
    limit       = list(gene = 5, cpd = 1)
  )
}










