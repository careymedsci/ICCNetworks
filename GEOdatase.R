
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



########################################################
# part1  Multi-cell Wnt/Fzd FC comparison and heatmap (GSE244678)
########################################################

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(data.table)
library(dplyr)
library(tibble) # for rownames_to_column / column_to_rownames

# -------------------- 0) Prepare online raw data files --------------------
exon_file <- "GSE244678_HGS_KD_cnts_exon.csv.gz"
exon_counts <- fread(exon_file) %>% as.data.frame()
rownames(exon_counts) <- exon_counts[[1]]
exon_counts <- exon_counts[,-1]

# gene ID extraction
gene_ids <- sapply(strsplit(rownames(exon_counts), "_"), `[`, 1)

# counts to gene-level
gene_counts <- exon_counts %>%
  mutate(gene_id = gene_ids) %>%
  group_by(gene_id) %>%
  summarise(across(everything(), sum)) %>%
  as.data.frame()

rownames(gene_counts) <- gene_counts$gene_id
gene_counts <- gene_counts[,-1]

# -------------------- 1) Select raw count files --------------------
b16_file <- file.choose()      # B16 counts
e0771_wt_file <- file.choose() # E0771 WT counts
e0771_inv_file <- file.choose()# E0771 Invivo counts

# -------------------- 2) Read counts --------------------
# B16
b16_counts <- read.table(b16_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
wt_cols <- grep("/wt/", colnames(b16_counts))
inv_cols <- grep("/brm/", colnames(b16_counts))
b16_counts <- b16_counts[, c(wt_cols, inv_cols)]
colnames(b16_counts) <- c("B16_WT1","B16_WT2","B16_WT3","B16_Inv1","B16_Inv2","B16_Inv3")

# E0771
e0771_wt <- read.table(e0771_wt_file, header=TRUE, row.names=1, sep="\t")
e0771_inv <- read.table(e0771_inv_file, header=TRUE, row.names=1, sep="\t")
e0771_counts <- cbind(e0771_wt, e0771_inv)
colnames(e0771_counts) <- c("E0771_WT1","E0771_WT2","E0771_Inv1","E0771_Inv2")

# -------------------- 3) B16 differential analysis --------------------
group_b16 <- factor(c(rep("WT",3), rep("Invivo",3)))
colData_b16 <- data.frame(Group=group_b16, row.names=colnames(b16_counts))
dds_b16 <- DESeqDataSetFromMatrix(countData=b16_counts, colData=colData_b16, design=~Group)
dds_b16 <- DESeq(dds_b16)
res_b16 <- lfcShrink(dds_b16, contrast=c("Group","Invivo","WT"), type="ashr")
res_b16_df <- as.data.frame(res_b16) %>%
  rownames_to_column("gene_name") %>%
  mutate(Direction = case_when(
    log2FoldChange > 1.5 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < -1.5 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not significant"
  ))

# -------------------- 4) Filter Fzd and Wnt genes --------------------
b16_sig_wntfzd <- res_b16_df %>%
  filter(Direction != "Not significant") %>%
  filter(grepl("^Fzd|^Wnt", gene_name, ignore.case=TRUE)) %>%
  select(gene_name, log2FoldChange)

wntfzd_genes <- b16_sig_wntfzd$gene_name

# -------------------- 5) Calculate FC for E0771 --------------------
group_e0771 <- factor(c(rep("WT",2), rep("Invivo",2)))
colData_e0771 <- data.frame(Group=group_e0771, row.names=colnames(e0771_counts))
dds_e0771 <- DESeqDataSetFromMatrix(countData=e0771_counts,
                                    colData=colData_e0771,
                                    design=~Group)
dds_e0771 <- DESeq(dds_e0771)
res_e0771 <- lfcShrink(dds_e0771, contrast=c("Group","Invivo","WT"), type="ashr")
res_e0771_df <- as.data.frame(res_e0771) %>%
  rownames_to_column("gene_name") %>%
  select(gene_name, log2FoldChange) %>%
  filter(gene_name %in% wntfzd_genes)

# ----------- 6) Calculate FC for RM (online data for B16-F10) -----------------
rm_counts_subset <- gene_counts[, c(3,4,5,7,9,10,11,12,13)]
colnames(rm_counts_subset) <- c(paste0("WT",1:4), paste0("Inv",1:5))
group_rm <- factor(c(rep("WT",4), rep("Invivo",5)))
colData_rm <- data.frame(Group=group_rm, row.names=colnames(rm_counts_subset))
dds_rm <- DESeqDataSetFromMatrix(countData=rm_counts_subset,
                                 colData=colData_rm,
                                 design=~Group)
dds_rm <- DESeq(dds_rm)
res_rm <- lfcShrink(dds_rm, contrast=c("Group","Invivo","WT"), type="ashr")
res_rm_df <- as.data.frame(res_rm) %>%
  rownames_to_column("gene_name") %>%
  select(gene_name, log2FoldChange) %>%
  filter(gene_name %in% wntfzd_genes)

# -------------------- 7) Merge FCs --------------------
fc_mat <- b16_sig_wntfzd %>%
  rename(B16 = log2FoldChange) %>%
  left_join(res_e0771_df %>% rename(E0771=log2FoldChange), by="gene_name") %>%
  left_join(res_rm_df %>% rename(RM=log2FoldChange), by="gene_name")

fc_mat_df <- fc_mat %>% column_to_rownames("gene_name")

# -------------------- 8) Heatmap --------------------
expr_mat <- as.matrix(fc_mat_df)
expr_mat[is.na(expr_mat)] <- 0

hm_colors <- colorRampPalette(c("blue","white","red"))(100)
pheatmap(expr_mat,
         cluster_rows=TRUE, cluster_cols=TRUE,
         color=hm_colors,
         main="Wnt/Fzd FC: B16 vs E0771 vs RM")

# -------------------- 9) Scatter plots with correlation --------------------
cor_e0771 <- cor(fc_mat_df$B16, fc_mat_df$E0771, method="spearman")
ggplot(fc_mat, aes(x=B16, y=E0771, label=rownames(fc_mat_df))) +
  geom_point() +
  geom_text_repel() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ggtitle(paste0("B16 vs E0771 FC, Spearman=", round(cor_e0771,2))) +
  xlab("B16 log2FC") + ylab("E0771 log2FC")

cor_rm <- cor(fc_mat_df$B16, fc_mat_df$RM, method="spearman")
ggplot(fc_mat, aes(x=B16, y=RM, label=rownames(fc_mat_df))) +
  geom_point() +
  geom_text_repel() +
  geom_abline(slope=1, intercept=0, linetype="dashed") +
  ggtitle(paste0("B16 vs RM FC, Spearman=", round(cor_rm,2))) +
  xlab("B16 log2FC") + ylab("RM log2FC")

# -------------------- 10) Save FC table --------------------
write.csv(fc_mat, "WntFzd_FC_B16_E0771_RM.csv", row.names=TRUE)



# ~~~~~~~~~~~~~~~~~~~~~part 2 single sample from online~~~~~~~~~~~~~~~~~~~
# ---------------------GSE294100
library(ggplot2)
library(dplyr)

file_path <- file.choose()
data <- read.table(file_path, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(data)[1:3] <- c("refgene","invivo","wt")

data <- data %>%
  mutate(
    refgene_clean = toupper(trimws(refgene)),   # clean names
    invivo = as.numeric(invivo),
    wt = as.numeric(wt),
    log2FC = log2((invivo + 1) / (wt + 1))
  )

# Genes of interest
genes_of_interest <- c("FZD9","WNT10A","WNT6","FZD3")
data_subset <- data %>%
  filter(refgene_clean %in% genes_of_interest)

# Fix order without losing any gene
data_subset$refgene <- factor(data_subset$refgene,
                              levels = data_subset$refgene[match(genes_of_interest, data_subset$refgene_clean)])

# Horizontal barplot
p <- ggplot(data_subset, aes(x = refgene, y = log2FC, fill = log2FC > 0)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(log2FC,2)), hjust = ifelse(data_subset$log2FC>0,-0.1,1.1)) +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="firebrick","FALSE"="steelblue"),
                    labels=c("Downregulated","Upregulated")) +
  theme_minimal() +
  labs(title="Log2 Fold Change of Selected Genes (Invivo vs WT)",
       x="Gene", y="log2(Invivo/WT)", fill="Direction") +
  theme(legend.position="top")

print(p)


