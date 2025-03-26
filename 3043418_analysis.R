#source functions###
source("3043418_functions.r")

#load necessary librarys####
library(ggplot2)
library(reshape2)
library(amap)
library(clusterProfiler)
library(eulerr)
library(ggrepel)
library(STRINGdb)
#BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

#Loading necessary tables####
em = read.table("EM.csv",sep="\t",header=TRUE,row.names = 1)
de_SvP = read.table("DE_Senes_vs_Prolif.csv",sep="\t",header=TRUE,row.names = 1)
de_SMvS = read.table("DE_Senes_MtD_vs_Senes.csv",sep="\t",header=TRUE,row.names = 1)
de_SMvP = read.table("DE_Senes_MtD_vs_Prolif.csv",sep="\t",header=TRUE,row.names = 1)
annotation = read.table("Human_Background_GRCh38.p13.csv",sep="\t",header=TRUE,row.names = 1)
ss = read.table("sample_sheet.csv",sep="\t",header=TRUE)

# Parse data####
#create master table for each DE comparisons 
names(annotation)[1] = "symbols"
master_SvP = create_master(em,de_SvP,annotation)
master_SMvS = create_master(em,de_SMvS,annotation)
master_SMvP = create_master(em,de_SMvP,annotation)

#create scaled expression table
em_scaled =data.frame(t(scale(t(em))))
em_scaled = na.omit(em_scaled)

#PCA plot and density to asses quality of samples####

#PCA plot
PCA_samples = plot_pca(em_scaled,ss$SAMPLE_GROUP)
PCA_samples


# plot Density plot
density_plot = plot_exp_density(em,3)
density_plot



#make mean_vale and log10p and sig columns ####
master_SvP$mean_value = rowMeans(master_SvP[,2:10])
master_SvP$mlog10p = -log10(master_SvP $p)
master_SvP$sig = as.factor(master_SvP$p.adj < 0.001 & abs(master_SvP$log2fold) > 2)

master_SMvS$mean_value = rowMeans(master_SMvS[,2:10])
master_SMvS$mlog10p = -log10(master_SMvS$p)
master_SMvS$sig = as.factor(master_SMvS$p.adj < 0.001 & abs(master_SMvS$log2fold) > 2)

master_SMvP$mean_value = rowMeans(master_SMvP[,2:10])
master_SMvP$mlog10p = -log10(master_SMvP$p)
master_SMvP$sig = as.factor(master_SMvP$p.adj < 0.001 & abs(master_SMvP$log2fold) > 2)

#make a list of significant genes#### 
sig_genes_SvP = row.names(master_SvP[master_SvP$sig == "TRUE", ])
sig_genes_SvP_up = subset(master_SvP, sig == "TRUE" & p.adj < 0.01 & log2fold >2)$symbols[1:10]
sig_genes_SvP_down =  subset(master_SvP, sig == "TRUE" & p.adj < 0.01 & log2fold <2)$symbols[1:10]
sig_genes_SvP_15 = rownames(master_SvP[order(master_SvP$p.adj), ])[1:15]

sig_genes_SMvS = row.names(master_SMvS[master_SMvS$sig == "TRUE",])
sig_genes_SMvS_up = subset(master_SMvS, sig == "TRUE" & p.adj < 0.01 & log2fold >2)$symbols[1:10]
sig_genes_SMvS_down =  subset(master_SMvS, sig == "TRUE" & p.adj < 0.01 & log2fold <2)$symbols[1:10]


sig_genes_SMvP = row.names(master_SMvP[master_SMvP$sig == "TRUE",])
sig_genes_SMvP_up = subset(master_SMvP, sig == "TRUE" & p.adj < 0.01 & log2fold >2)$symbols[1:10]
sig_genes_SMvP_down =  subset(master_SMvP, sig == "TRUE" & p.adj < 0.01 & log2fold <2)$symbols[1:10]


#volcano plots for each comparison
volcano_svp = plot_volcano(master_SvP,0.001,2)
volcano_svp

volcano_smvs =plot_volcano(master_SMvS,0.001,2)
volcano_smvs

volcano_smvp = plot_volcano(master_SMvP,0.001,2)
volcano_smvp

#MA plot for each comparison
Ma_senes_prolif = plot_MA(master_SvP,0.001,2)
Ma_senes_prolif

Ma_senesmtd_senes = plot_MA(master_SMvS,0.001,2)
Ma_senesmtd_senes

Ma_senesmtd_prolif = plot_MA(master_SMvP,0.001,2)
Ma_senesmtd_prolif


#create a master table with DEGs in  data####
# create master table with  svp and senes mtd vs prolif comparison 
master_SvSM= merge(master_SvP, master_SMvS, by = "symbols", suffixes = c(".SvP", ".SMvS"), all = TRUE)
rownames(master_SvSM) = master_SvSM$symbols
master_SvSM = na.omit(master_SvSM)

#plot fold vs fold plot comparisons for Senes vs senes mtd comparison
fold_changes = data.frame(
  gene = master_SvSM$symbols,
  log2FC_SvP = master_SvSM$log2fold.SvP,
  log2FC_SMvS = master_SvSM$log2fold.SMvS,
  sig_SvP = master_SvSM$sig.SvP,
  sig_SMvS = master_SvSM$sig.SMvS)


fold_changes$Sig_Status = ifelse(fold_changes$sig_SvP == TRUE & fold_changes$sig_SMvS == TRUE , "All Significant",
                                 ifelse(fold_changes$sig_SvP == TRUE, "Senes vs Prolif Only",
                                        ifelse(fold_changes$sig_SMvS == TRUE, "Senes_MtD vs Senes Only", "Non-sig")))

fold_plot = ggplot(fold_changes, aes(x = log2FC_SvP, y = log2FC_SMvS)) +
  geom_point(aes(color = Sig_Status), size = 3) +  
  scale_color_manual(values = c("All Significant" = "red3", "Senes vs Prolif Only" = "darkblue", 
                                "Senes_MtD vs Senes Only" = "darkolivegreen", "Non-sig" = "black")) +
  
  labs(title= "Fold vs Fold Plot",x = "Log2 Fold Change (Senes vs prolif)", y = "Log2 Fold Change (Senes_MtD v Senes)", color = "Significance") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
        axis.text.x = element_text(size = 12, hjust = 1), 
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        strip.text = element_text(size = 14, face = "bold"), 
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        panel.grid.major = element_line(linewidth = 0.4, color = "gray80"), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(size = 1, fill = NA)) 


fold_plot


#create expression scaled table for genes in master_SvSM table 
em_SvSM = master_SvSM[,c(3:11)]

#remove suffixes from sample names 
colnames(em_SvSM) = gsub("\\.SvP$|\\.SMvS$|\\.SMvP$", "", colnames(em_SvSM))

#scale expression table 
em_SvSM = data.frame(t(scale(t(em_SvSM))))
em_SvSM = na.omit(em_SvSM)

#make a master with all comparisons####
master_SMvP_renamed = master_SMvP

#give suffix .SMvP for all columns except symbols
colnames(master_SMvP_renamed)[-11] = paste0(colnames(master_SMvP_renamed)[-11], ".SMvP")  

master_all = merge(master_SvSM, master_SMvP_renamed, by = "symbols",  all = TRUE)
rownames(master_all) = master_SvSM$symbols
master_all = na.omit(master_all)

#venn diagram of comparisons #### 
venn_data = list("Senes vs prolif" = sig_genes_SvP, "SenesMtd vs Senes" = sig_genes_SMvS, "SenesMtd vs Prolif" = sig_genes_SMvP) 
plot(euler(venn_data, shape = "circle"), quantities = TRUE) 

venn_data_svp_vs_smvs = list("Senes vs prolif" = sig_genes_SvP, "SenesMtD vs Senes" = sig_genes_SMvS)
plot(euler(venn_data_svp_vs_smvs, shape = "circle"), quantities = TRUE) 

venn_data_svp_vs_smvp =  list("Senes vs prolif" = sig_genes_SvP, "SenesMtd vs Prolif" = sig_genes_SMvP) 
plot(euler(venn_data_svp_vs_smvp, shape = "circle"), quantities = TRUE) 

venn_data_smvs_vs_smvp = list("Senesmtd vs senes" = sig_genes_SMvS, "SenesMtD vs Prolif" = sig_genes_SMvP)
plot(euler(venn_data_smvs_vs_smvp, shape = "circle"), quantities = TRUE) 


#get all significant genes in atleast one comparisons
sig_genes_all = subset(master_all, sig.SvP == TRUE|sig.SMvS == TRUE|sig.SMvP == TRUE )$symbols

# Perform hypergeometric test
group1 = length(sig_genes_SvP)  
group2 = length(sig_genes_SMvS)
group3 = length(sig_genes_SMvP)
total = nrow(master_all) 
overlap = length(sig_genes_all) 

p_value =  phyper(overlap - 1, total, total - (group1 + group2 + group3), overlap, lower.tail = FALSE)
# Print result
print(paste("Overall Overlap p-value:", p_value))

p_value_1_vs_2 = phyper(overlap-1, group2, total-group2, group1,lower.tail= FALSE)
print(paste("Overall Overlap p-value:", p_value_1_vs_2))

p_value_1_vs_3 = phyper(overlap-1, group3, total-group3, group1,lower.tail= FALSE)
print(paste("Overall Overlap p-value:", p_value_1_vs_2))


#heat map for only significant genes atleast one comparison####
heatmap_sig = plot_heatmap(em_SvSM[sig_genes_all,],row_cluster = TRUE,col_cluster = FALSE)
heatmap_sig

# find  and plot signatures in heatmap####

#1. highly down in proliferation and senes highly up in senes mitochndria depleted 
signature_1 = subset(master_all,sig.SMvP == T & sig.SMvS == T &log2fold.SMvP > 0  & log2fold.SMvS > 0)$symbols
signature_1_plots= plot_signature(signature_1,em_SvSM,ss$SAMPLE_GROUP,org.Hs.eg.db,"up_senesmtdvProlif")
signature_1_plots$bar_plot
signature_1_plots$cnet_plot
ora_upSMvP = signature_1_plots$ora_results_table

#2. highly up mitochondria depleted highly down in senes 
signature_2 = subset(master_all,sig.SMvS ==T & log2fold.SMvS > 0 )$symbols
signature_2_plots= plot_signature(signature_2,em_SvSM,ss$SAMPLE_GROUP,org.Hs.eg.db,"down_senes_vSenesmtd")
signature_2_plots$bar_plot
ora_up_SmvS = signature_2_plots$ora_results_table
#here the cellular response to oxygen is increased 

#3. highly up in senes highly down in mitochondria depleted and prolif 
signature_3 = subset(master_all,sig.SMvS ==T & sig.SvP == T & log2fold.SMvS < 0 & log2fold.SvP > 0)$symbols
signature_3_plots= plot_signature(signature_3,em_SvSM,ss$SAMPLE_GROUP,org.Hs.eg.db,"up_senes_vSenesmtd")
signature_3_plots$bar_plot
ora_down_SvSM = signature_3_plots$ora_results_table
# these are ones involed in ECM 

#4.highly up in senes 
signature_4 = subset(master_all,sig.SvP == T & log2fold.SvP > 0 )$symbols
signature_4_plots = plot_signature(signature_4,em_SvSM,ss$SAMPLE_GROUP,org.Hs.eg.db,"up_senes_vprolif")
signature_4_plots$bar_plot
signature_4_plots$metagene
ora_upSvP = signature_4_plots$ora_results_table
genes_ora_upSvP_1 = enriched_heat_boxplot_string(ora_upSvP,em_SvSM,ss$SAMPLE_GROUP,"apoptotic cell clearance")
genes_ora_up_SvP_2 = enriched_heat_boxplot_string(ora_upSvP,em_SvSM,ss$SAMPLE_GROUP,"phagocytosis")
#this involves immune and ECM pathways 

#5. highly up in prolif  compared to both   
signature_5 = subset(master_all,sig.SMvP == TRUE & sig.SvP == T & log2fold.SMvP <0 & log2fold.SvP <0)$symbols
signature_5_plots = plot_signature(signature_5,em_SvSM,ss$SAMPLE_GROUP,org.Hs.eg.db,"up_prolif_senesmtd")
signature_5_plots$bar_plot
ora_down_SMvP= signature_5_plots$ora_results_table
#These involve cell cycle pathways 

#6. highly up in prolif highly down in senes 
signature_6 = subset(master_all,sig.SvP == T & log2fold.SvP < 0)$symbols
signature_6_plots = plot_signature(signature_6,em_SvSM,ss$SAMPLE_GROUP,org.Hs.eg.db,"up_prolif_vsenes")
signature_6_plots$bar_plot
ora_down_SvP= signature_6_plots$ora_results_table
#this is also cell division pathways 

#do gsea for each des####
#gsea for SvP
# we want the log2 fold change  
gsea_input_SvP = master_all$log2fold.SvP
# add gene names the vector 
names(gsea_input_SvP) = row.names(master_all) 
# omit any NA values  
gsea_input_SvP = na.omit(gsea_input_SvP) 
# sort the list in decreasing order (required for clusterProfiler) 
gsea_input_SvP = sort(gsea_input_SvP, decreasing = TRUE) 


gse_results_SvP = gseGO(geneList=gsea_input_SvP,  
                         ont ="BP",  
                         keyType = "SYMBOL",  
                         nPerm = 10000,  
                         minGSSize = 3,  
                         maxGSSize = 800,  
                         pvalueCutoff = 0.001,  
                         verbose = TRUE,  
                         OrgDb = org.Hs.eg.db,  
                         pAdjustMethod = "none")

#plot ridge plot 
ggp_SvP = ridgeplot(gse_results_SvP, showCategory = 20)
ggp_SvP 

# create data frame from gsea results 
description = gse_results_SvP$Description
p.adj = gse_results_SvP$p.adjust
NES = gse_results_SvP$NES
gene_sets = gse_results_SvP$core_enrichment

gse_SvP = data.frame(description, p.adj, NES,gene_sets, row.names = description)
gse_SvP = gse_SvP[order(gse_SvP$p.adj), ]

#analyse interested gene sets in the gsea result btwn senes v prolif ####

#5 up highly positive NES of interest heat map and box plots  
genes1_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"humoral immune response")
genes2_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"G protein-coupled receptor signaling pathway")
genes3_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"leukocyte chemotaxis")
genes4_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"B cell mediated immunity")
genes5_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"acute inflammatory response")
genes6_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"complement activation, classical pathway")

#most immune genes are upregulated in senes compared to prolif

#5 down highly genative NES of interest heat map and box plots  
#genes associated with cell cycle shows downregulation in senes resulting in a cell proliferation arrest
genes7_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"cell cycle DNA replication initiation")
genes8_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"positive regulation of cell cycle checkpoint")
genes9_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"meiotic cell cycle phase transition")
genes10_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"positive regulation of mitotic cytokinesis")
genes11_SvP = enriched_heat_boxplot_string(gse_SvP,em_SvSM,ss$SAMPLE_GROUP,"regulation of mitotic sister chromatid separation")


#GSEA for SM v S ####
# we want the log2 fold change  
gsea_input_SMvS = master_all$log2fold.SMvS 
# add gene names the vector 
names(gsea_input_SMvS) = row.names(master_all) 
# omit any NA values  
gsea_input_SMvS = na.omit(gsea_input_SMvS) 
# sort the list in decreasing order (required for clusterProfiler) 
gsea_input_SMvS = sort(gsea_input_SMvS, decreasing = TRUE) 


gse_results_SMvS = gseGO(geneList=gsea_input_SMvS,  
                         ont ="BP",  
                         keyType = "SYMBOL",  
                         nPerm = 10000,  
                         minGSSize = 3,  
                         maxGSSize = 800,  
                         pvalueCutoff = 0.001,  
                         verbose = TRUE,  
                         OrgDb = org.Hs.eg.db,  
                         pAdjustMethod = "none")


ggp_SMvS = ridgeplot(gse_results_SMvS, showCategory = 20)
ggp_SMvS

# obtain objects from GSEA results 
description = gse_results_SMvS$Description
p.adj = gse_results_SMvS$p.adjust
NES = gse_results_SMvS$NES
gene_sets = gse_results_SMvS$core_enrichment

# Create a data frame
gse_SMvS = data.frame(description, p.adj, NES, gene_sets, row.names = description)
gse_SMvS = gse_SMvS[order(gse_SMvS$p.adj), ]

#most pathways involved in immune such as chemotaxis, humoral immune response, leucocyte chemotaxis etc are seem enriched in this comparison as well
#This suggest strong immune regulation in cells that become senes

#analyse interested gene sets in the gsea result btwn senesmtd v senes ####
genes_1_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"positive regulation of cytokine production")
genes2_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"inflammatory response")
genes3_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"integrin-mediated signaling pathway")
genes4_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"glycoprotein metabolic process")
genes5_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"regulation of apoptotic cell clearance")
genes6_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"leukocyte migration")
genes7_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"external encapsulating structure organization")

genes8_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"NADH metabolic process")
genes9_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"stress granule assembly")
genes10_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"nuclear-transcribed mRNA catabolic process, nonsense-mediated decay")
genes11_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"cell fate specification")
genes12_SMvS = enriched_heat_boxplot_string(gse_SMvS,em_SvSM,ss$SAMPLE_GROUP,"GABAergic neuron differentiation")


#the metabolic processes aree down regulated in senes mtd due to inavailability of energy producing mitochondrias 
#cytokine production is highly downregulated in senes mtd comparatively cells become les inflamatory
#whereas NADH metabolic process is incresed in senes mtd to compensate normal ATP energy production
#The absence of mitochondria may result in metabolic shifts and a reduction in inflammatory signaling, 
#while the cell prioritizes survival in the face of compromised energy production and function.
#when mitochondria is removed cells shifts to a more intensive metabolic state looking for lternate pathways to get energy and becomes more adaptive and survival state 
#senescels cells give more inflamatory signals that indicate damage to cells  


#GSEA for SM v prolif ####
# we want the log2 fold change  
gsea_input_SMvP = master_all$log2fold.SMvP
# add gene names the vector 
names(gsea_input_SMvP) = row.names(master_all) 
# omit any NA values  
gsea_input_SMvP = na.omit(gsea_input_SMvP) 
# sort the list in decreasing order (required for clusterProfiler) 
gsea_input_SMvP = sort(gsea_input_SMvP, decreasing = TRUE) 


gse_results_SMvP = gseGO(geneList=gsea_input_SMvP,  
                         ont ="BP",  
                         keyType = "SYMBOL",  
                         nPerm = 10000,  
                         minGSSize = 3,  
                         maxGSSize = 800,  
                         pvalueCutoff = 0.001,  
                         verbose = TRUE,  
                         OrgDb = org.Hs.eg.db,  
                         pAdjustMethod = "none")


ggp_SMvP = ridgeplot(gse_results_SMvP,showCategory = 20)
ggp_SMvP


description = gse_results_SMvP$Description
p.adj = gse_results_SMvP$p.adjust
NES = gse_results_SMvP$NES
gene_sets = gse_results_SMvP$core_enrichment

# Create a data frame
gse_SMvP = data.frame(description, p.adj, NES, gene_sets, row.names = description)
gse_SMvP = gse_SMvP[order(gse_SMvP$p.adj), ]

#analyse interested gene sets in the gsea result btwn senesmtd  v prolif ####
#up
genes_1_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"cellular response to hypoxia")
genes_2_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"secretion")
genes_3_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"epidermis development")
genes_4_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"cellular response to decreased oxygen levels")
genes_5_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"intracellular chemical homeostasis")
genes_11_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"organ or tissue specific immune response")
genes_12_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"cellular response to decreased oxygen levels")


#down
genes_6_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"regulation of meiotic cell cycle phase transition")
genes_7_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"positive regulation of mitotic cell cycle spindle assembly checkpoint")
genes_8_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"positive regulation of chromosome separation")
genes_9_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"cell cycle DNA replication")
genes_10_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"cell cycle G2/M phase transition")

genes_13_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"chromosome segregation")

genes_13_SMvP = enriched_heat_boxplot_string(gse_SMvP,em_SvSM,ss$SAMPLE_GROUP,"nuclear division")


#The cell does not regain its proliferative function but become more adaptive and become differentiated 

#venn diagram to show common pathways between conditions 
pathways_SvP = row.names(gse_SvP)
pathways_SMvS = row.names(gse_SMvS)
pathways_SMvP = row.names(gse_SMvP)
venn_data = list("Senes vs prolif pathways" = pathways_SvP, "SenesMtd vs Senes pathways" = pathways_SMvS, "SenesMtd vs Prolif pathways" = pathways_SMvP) 
plot(euler(venn_data, shape = "circle"), quantities = TRUE) 

#analysing some genes of interests 

#chemokines and cytokines such as :
candidate_genes_1 = as.vector(c("C7","C3","CXCL14","CXCL10","CCL13","CCL5","SAA1","IL16","IL6","IL4","IL6R","SAA2","IL20RB"))
multi_box_plot_candidate_genes = multi_boxplot_facet(candidate_genes_1,em_SvSM,ss$SAMPLE_GROUP)
multi_box_plot_candidate_genes 
#-The cells become more inflammatory as they become senense and as the mitochondria is removed they become more in a survival state and in a adaptive state with increased secretory and differenciation 
#and also tries to reverese senesce inflamatiory signalling and imbalenses such as ion imbalance -#





