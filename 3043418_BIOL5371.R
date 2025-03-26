#HYPOTHESIS
#Ho:The blood can not be used be used to differentiate between GOUT and SEPSIS based on gene expression profile 
#Ha:The blood can be used to differentiate between GOUT and SEPSIS based on gene expression profile 

#import the datasets into R
expr_value=read.table("C:/Users/liyan/Desktop/Glasgow/ststistics for bioinformatics 2425/data for assesment report/Expression_Table.csv",header=TRUE,row.names=1,sep="\t")
sample_info=read.table("C:/Users/liyan/Desktop/Glasgow/ststistics for bioinformatics 2425/data for assesment report/Sample_Information.csv",header=TRUE,row.names=1,sep="\t")
annotation=read.table("C:/Users/liyan/Desktop/Glasgow/ststistics for bioinformatics 2425/data for assesment report/Annotations.csv",header=TRUE,row.names=1,sep="\t")
DE_gout_vs_healthy=read.table("C:/Users/liyan/Desktop/Glasgow/ststistics for bioinformatics 2425/data for assesment report/DE_GOUT_vs_HC.csv",header=TRUE,row.names=1,sep="\t")
De_sa_vs_healthy=read.table("C:/Users/liyan/Desktop/Glasgow/ststistics for bioinformatics 2425/data for assesment report/DE_SA_vs_HC.csv",header=TRUE,row.names=1,sep="\t")

#annotating differential and expression tables 
expression_annotated=merge(annotation,expr_value,by.x=0,by.y=0)
names(expression_annotated)[1]="GENE_ID"
DE_gout_vs_healthy_annotated=merge(annotation,DE_gout_vs_healthy,by.x=0,by.y=0)
names(DE_gout_vs_healthy_annotated)[1]="GENE_ID"
De_sa_vs_healthy_annotated=merge(annotation,De_sa_vs_healthy,by.x=0,by.y=0)
names(De_sa_vs_healthy_annotated)[1]="GENE_ID"


#summary of clinical data 

sample_info$SEX=as.factor(sample_info$SEX )#converting 'sex' and 'sample group' to factor 
sample_info$SAMPLE_GROUP=as.factor(sample_info$SAMPLE_GROUP)
summary(sample_info)

#There are 14 females nd 13 males in the study and from each group 9 individuals are sampled resulting in a total of 27 observations  in the total sample min neutrophil count is 1.8 median is 7.6 max is 14.9 

#ANOVA test to see if neutrophil count is well matched between healthy,gout and SA samples and TukeyHSD test to understand the difference between each group

aggregate(cbind(NEUTROPHILS) ~ SAMPLE_GROUP, data = sample_info, summary)
anova_neutrophils=aov(NEUTROPHILS ~ SAMPLE_GROUP, data = sample_info)
summary(anova_neutrophils)
TukeyHSD(anova_neutrophils)

#Chi-squared test to see is sex is well matched between H,G&SA

table(sample_info$SEX, sample_info$SAMPLE_GROUP)
chisq.test(table(sample_info$SEX, sample_info$SAMPLE_GROUP))

#the result says the neutrophil count have significant difference between each group this may be due to immune responses due to the disease environment the neutrophil count particularly varies between SA vs healthy group.

#identifying significant genes from the differentially expressed genes between gout and healthy and sa and healthy
sig_gout=subset(DE_gout_vs_healthy_annotated,p.adj<0.05& abs(log2Fold)>1)
#arranging significant genes in ascending order of p value to get most significant gene
sig_gout=sig_gout[order(sig_gout$p.adj), ]

sig_SA=subset(De_sa_vs_healthy_annotated,p.adj<0.05&abs(log2Fold)>1)
sig_SA=sig_SA[order(sig_SA$p.adj),]

row.names(sig_gout)=sig_gout[,"GENE_ID"]#setting row names to"GENE_ID"
row.names(sig_SA)=sig_SA[,"GENE_ID"]

#combining all significant genes 
sig_genes=rbind(sig_gout,sig_SA)
samples=row.names(expr_value)
genes=names(expr_value)
#creating expression table for al significant genes 
all_sig_genes_expr=expr_value[samples,genes]
all_sig_S_G=all_sig_genes_expr[, grepl("GOUT|SEPSIS", colnames(all_sig_genes_expr))]
all_sig_S_Gt=as.data.frame(t(all_sig_S_G))  

#Check for zero variance
non_constant_cols = apply(all_sig_S_Gt, 2, function(x) var(x, na.rm = TRUE) > 1e-8)
#create data without non zero variance
filtered_data = all_sig_S_Gt[, non_constant_cols]

#creating PCA data

PCA_data = prcomp(filtered_data, scale. = TRUE)
PCA_data= as.data.frame(PCA_data$x)

#collecting only Gout and Sepsis samples 
sample_names = rownames(filtered_data)
group_vector = ifelse(grepl("GOUT", rownames(filtered_data)), "GOUT",      
                      ifelse(grepl("SEPSIS", rownames(filtered_data)), "SEPSIS", NA))
PCA_data$Group= group_vector

# Plot PCA 
library(ggplot2)
ggplot(PCA_data, aes(x = PC1, y = PC2, color = Group)) + geom_point(size = 3) +
  labs(title = "PCA of Gene Expression Data for all sig genes Hc vs G & HC vs SA ", x = "Principal Component 1", y = "Principal Component 2")

#find overlapping genes between significant genes 
overlap_genes=intersect(sig_gout$GENE_ID,sig_SA$GENE_ID)
length(overlap_genes)
print(overlap_genes)
#THERE ARE only 8 GENES THAT ARE SIMILIAR IN DEGS BETWEEN HC VS GOUT AND HC VS SA THIS TELLS THAT THERE ARE SIGNIFICANT TRANSCRIPTIONAL DIFFERENCE BETWEEN GOUT AND SA THE SIGNIFICANT GENES ARE NOT SIMILIAR BETWEEN THEM 

#the most significantly differential genes in HC VS SA and HC VS GOUT by filtering top 10 genes in the significant data set
top10_gout = head(sig_gout[order(sig_gout$p.adj), ], 10)
top10_sa = head(sig_SA[order(sig_SA$p.adj), ], 10)

#intersecting to see if there are any common genes 
intersect(top10_gout$GENE_ID,top10_sa$GENE_ID)

#there are no common top 10 genes we can see if these genes differ in their expression levels in each group gout and sepsis to see if we can differentiate between these conditions 

# Add a 'group' column to the top 10 data 
top10_gout$group= "Gout_vs_HC"
top10_sa$group="SA_vs_HC"

# Combine top 10 data 
top10_combined = rbind(top10_gout, top10_sa)

#creating an expression data of most significant genes 
#create neccessary objects 
object2=names(expr_value)
object1.1=top10_gout$GENE_ID
object3.3=top10_sa$GENE_ID
object4.4=top10_combined$GENE_ID

#Combining objects to create necessary data expression data 
top10_merged=expr_value[object4.4,object2]
top10_gout_exp=expr_value[object1.1,object2]
top10_sa_exp=expr_value[object3.3,object2]

#create annotated data 
top_gout_annonated=merge(annotation,top10_gout_exp,by.x=0,by.y=0)
top10_sa_annotated=merge(annotation,top10_sa_exp,by.x=0,by.y=0)

#PCA ON TOP 10 GENES OF HC vsGOUT and HC vs SEPSIS

#create PCA data
pca_data_1=prcomp(t(top10_merged), scale. = TRUE)
pca_data_1= as.data.frame(pca_data_1$x)
pca_data_1$Group= sample_info$SAMPLE_GROUP

# Plot PCA 
ggplot(pca_data_1, aes(x = PC1, y = PC2, color = Group)) + geom_point(size = 3) +
  labs(title = "PCA of Gene Expression Data for top 10 genes ", x = "Principal Component 1", y = "Principal Component 2") 

#we can see there are distinct clustering of these significant genes in each groups

#preparing files to plot box plots showing expression between groups HC, GOUT and  SEPSIS for all top significant genes in hc vs gout and hc vs sa 

#transforming the table
top10_merged_trans=as.data.frame(t(top10_merged)) 

#combining  columns for sex and neutrophil and sample group
top10_merged_trans$sex=sample_info$SEX
top10_merged_trans$NEUTROPHILS=sample_info$NEUTROPHILS 
top10_merged_trans$group=sample_info$SAMPLE_GROUP

#creating an annotation file for top 10 significant genes
top10_merged_annotation=merge(annotation,top10_merged,by.x=0,by.y=0) 
row.names(top10_merged_annotation)=top10_merged_annotation$symbol
names(top10_merged_annotation)="GENE_ID"

#plotting box plot
library(ggplot2)
for (gene in names(top10_merged_trans)[-c(ncol(top10_merged_trans)-1,ncol(top10_merged_trans)-2,ncol(top10_merged_trans))]) {
  boxplot(top10_merged_trans[[gene]] ~ top10_merged_trans$group,
          main=paste("Gene Expression of", gene, "by Group"),
          xlab="Group",
          ylab="Expression Level")
  system.time(2)
}

#most the significant genes in HC vs sepsis are up-regulated in sepsis group these genes are almost not expressed in gout and healthy control

#To see if gene expression is effected by clinical variables like sex and neutrophil count
#1.Main hypothesis 
#Ho:gene expression is not effected by and neutrophil count and sex 
#Ha: gene expression is effected by neutrophil count and sex 
#2.Hypothesis for sex
#Ho:gene expression is not effected by sex
#Ha:Gene expression is effected by sex 
#3.Hypothesis for neutrophils
#Ho:gene expression is not effected by neutrophil count 
#Ha:gene expression is effected by neutrophil count 

#creating data to fit GLM
top10_merged_trans=as.data.frame(t(top10_merged))
top10_merged_trans$sex=sample_info$SEX
top10_merged_trans$NEUTROPHILS=sample_info$NEUTROPHILS
top10_merged_trans$group=sample_info$SAMPLE_GROUP

#fitting linear model  all significant top20 genes 
#creating a list to store results 
results_lm = data.frame(GENE_ID = character(),
                        p_value_sex = numeric(),
                        p_value_neutrophils = numeric(),
                        stringsAsFactors = FALSE)

#Looping thriugh all genes 
for (gene in names(top10_merged_trans)[-c(ncol(top10_merged_trans)-1, ncol(top10_merged_trans)-2, ncol(top10_merged_trans))]) {
  
  # Linear Model for sex and NEUTROPHILS
  model_sig = lm(as.numeric(top10_merged_trans[[gene]]) ~ sex + NEUTROPHILS, data = top10_merged_trans)
  summary_model_sig = summary(model_sig)
 
  
  # Print results for the current gene
  cat("Gene ID:", gene, "\n")
  print(summary_model_sig)
  
  
  # Extract p-values for sex and neutrophil counts
  p_value_sex = coef(summary_model_sig)[2, 4]  
  p_value_neutrophils = coef(summary_model_sig)[3, 4]  # p-value for NEUTROPHILS
  
   # Combine results
  results_lm = rbind(results_lm, data.frame(GENE_ID = gene,
                                            p_value_sex = p_value_sex,
                                            p_value_neutrophils = p_value_neutrophils))
  
 }

# Apply multiple testing correction using BH method  
results_lm$corrected_p_value_sex = p.adjust(results_lm$p_value_sex, method = "BH")
results_lm$corrected_p_value_neutrophils = p.adjust(results_lm$p_value_neutrophils, method = "BH")


# Display the results
print(results_lm)

#we can see the neutrophil counts effect the gene expressions in most groups in sepsis and some in gout and sex does not effect gene expression rejecting main and neutrophil null hypothesis while accepting sex null hypothesis
#collecting only gout and sepsis sample expressions 
top10_merged_G_S=top10_merged[, grepl("GOUT|SEPSIS", colnames(top10_merged))] 

#creating transformed data 
top10_merged_G_St=as.data.frame(t(top10_merged_G_S))

#Creating a vector of sample names 'GOUT' & 'SEPSIS'
sample_names = rownames(top10_merged_G_St)
group_vector = ifelse(grepl("GOUT", rownames(top10_merged_G_St)), "GOUT",     
                      ifelse(grepl("SEPSIS", rownames(top10_merged_G_St)), "SEPSIS", NA))

#assigning this vector for the samples in a new column 'group' to compare between them
top10_merged_G_St$group = group_vector                                        
top10_merged_G_St$group = as.factor(top10_merged_G_St$group)

#performing t-test for these top20  genes  between gout and sepsis samples to find most significant genes in gout vs sepsis 
#Ho:There are no genes varying significantly in their expression between gout and sepsis
#Ha:There are genes varying significantly in their expression between gout and sepsis
#create dataframe  to store results
t_test_results =data.frame(Gene_ID = character(),
                           mean_Gout = numeric(),
                           mean_sepsis = numeric(),
                           p_value_ttest = numeric(),
                           tatistic = numeric(),
                           conf_int_lower = numeric(),
                           conf_int_upper = numeric(),
                           stringsAsFactors = FALSE) 
#loop through all 20 genes 
for (gene in colnames(top10_merged_G_St)[-ncol(top10_merged_G_St)]) {
  # Extract the expression data for the current gene
  gout_expr =top10_merged_G_St[top10_merged_G_St$group == "GOUT", gene]
  sepsis_expr = top10_merged_G_St[top10_merged_G_St$group == "SEPSIS", gene]
  
  # Perform the t-test
  t_test_result = t.test(gout_expr, sepsis_expr)
  
  #calculate mean
  mean_gout = mean(gout_expr, na.rm = TRUE)
  mean_sepsis = mean(sepsis_expr, na.rm = TRUE)
  
  # Store the results
  t_test_results= rbind(t_test_results, data.frame(
    Gene_ID = gene,
    mean_Gout = mean_gout,
    mean_Sepsis = mean_sepsis,
    p_value_ttest = t_test_result$p.value,
    statistic = t_test_result$statistic,
    conf_int_lower = t_test_result$conf.int[1],
    conf_int_upper = t_test_result$conf.int[2]
  ))
}

#correcting P value 
t_test_results$adj_p_value = p.adjust(t_test_results$p_value_ttest, method = "BH")
# View the results
print(t_test_results)

#selecting top 10 genes from the t-test
top_genes_SA_GOUT = t_test_results[order(t_test_results$adj_p_value), ][1:10, ]

#since we see a significant difference in gene expression between SEPSIS and GOUT we re just the null hypothesis 

#CREATING HEATMAP AND BOXPLOT FOR EXPRESSION BETWEEEN SEPSIS AND GOUT FOR TOP 10 GENES 

#creating expression data
objects_sig=top_genes_SA_GOUT$Gene_ID
sig_G_SA_EXP=expr_value[objects_sig,object2]
sig_G_SA_EXP_anotated=merge(annotation,sig_G_SA_EXP,by.x=0,by.y=0)

#getting only samples for GOUT and SEPSIS in them 
gout_sepsis_exp = sig_G_SA_EXP[, grepl("GOUT|SEPSIS", colnames(sig_G_SA_EXP))]
gout_sepsis_expt=as.data.frame(t(gout_sepsis_exp))

#Creating a vector of sample names 'GOUT' & 'SEPSIS'
sample_names = rownames(gout_sepsis_expt)
group_vector = ifelse(grepl("GOUT", rownames(gout_sepsis_expt)), "GOUT",      
                      ifelse(grepl("SEPSIS", rownames(gout_sepsis_expt)), "SEPSIS", NA))

#assigning this vector for the samples in a new column group to compare between them
gout_sepsis_expt$group = group_vector                                        
gout_sepsis_expt$group = as.factor(gout_sepsis_expt$group) #converting group from character into a factor 


#plotting box plot for all genes significant in GOUT vs SA 
library(ggplot2)

for (gene in names(gout_sepsis_expt)[-ncol(gout_sepsis_expt)]) {  
  boxplot(gout_sepsis_expt[[gene]] ~ gout_sepsis_expt$group,
          main=paste("Gene Expression of", gene, "by Group"),
          xlab="Group",
          ylab="Expression Level")
  Sys.sleep(1)
}

# Heatmap for genes significant in GOUT VS SA 
library(pheatmap)
heatmap_data_2= as.matrix(gout_sepsis_exp)

# standardize the data 
heatmap_data_scaled_2 = t(scale(t(heatmap_data_2)))

# Create a heatmap
pheatmap(heatmap_data_scaled_2,
         cluster_rows = TRUE,  # Cluster genes
         cluster_cols = TRUE,  # Cluster samples
         show_rownames = TRUE,  # Show gene names
         show_colnames = TRUE,  # Show sample names
         main = "Heatmap of top 10 Significant Gene Expression",
         color = colorRampPalette(c("blue", "white", "red"))(50))


#plotting a PCA for Gout vs sepsis



#create PCA data
pca_data_2=prcomp(t(top10_merged_G_S), scale. = TRUE)
pca_data_2= as.data.frame(pca_data_2$x)

sample_names = rownames(top10_merged_G_St)
group_vector = ifelse(grepl("GOUT", rownames(top10_merged_G_St)), "GOUT",      
                      ifelse(grepl("SEPSIS", rownames(top10_merged_G_St)), "SEPSIS", NA))
pca_data_2$Group= group_vector

# Plot PCA 
ggplot(pca_data_2, aes(x = PC1, y = PC2, color = Group)) + geom_point(size = 3) +
  labs(title = "PCA of Gene Expression Data for top 10 genes G vs SA ", x = "Principal Component 1", y = "Principal Component 2") 


cat('\n____________________________________NULL HYPOTHESIS REJECTED_________________________________________________\n')








 


