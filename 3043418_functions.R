

#function to create master table####
create_master = function(em_table,de_table,anotation_table){
  
  #create master temp
  master_temp = merge(em_table,anotation_table,by.x=0,by.y=0)
  master=merge(master_temp,de_table,by.x=1,by.y=0)
  
  #make the row names as gene symbols
  row.names(master) = master[,"symbols"]
  names(master)[1]="Gene_ID"
  
  #omit NA values
  master = na.omit(master)
  return(master)
  
}

#function to plot volcano####

plot_volcano = function(de_table, p_threshold, fold_threshold) 
{ 
  # Obtain necessary tables 
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) > fold_threshold)  
  de_table_sig_up = subset(de_table_sig, p.adj < p_threshold & (log2fold) >fold_threshold) 
  de_table_sig_down = subset(de_table_sig, p.adj < p_threshold & (log2fold) < (fold_threshold))
  de_table_sig_up_top5 =   de_table_sig_up[1:5,] 
  de_table_sig_down_top5 =  de_table_sig_down[1:5,] 
  
  # Plot Volcano
  ggp = ggplot(de_table,aes(x=log2fold, y=mlog10p)) + 
    geom_point(aes(colour="a")) +  
    geom_point(data= de_table_sig_up, aes(colour="b")) +  
    geom_point(data= de_table_sig_down, aes(colour="c")) +
    labs(title="", x="-log10p", y="log2fold") +
    geom_label_repel(data=de_table_sig_up_top5, aes(label = symbols,colour ="b"), show.legend = FALSE,max.overlaps = Inf) +
    geom_label_repel(data=de_table_sig_down_top5,aes(label = symbols,colour = "c"), show.legend = FALSE,max.overlaps = Inf) + 
    scale_color_manual(values = c("black","red3","darkblue"),labels =c("No change","Upregulated","Downregulated"),name ="")+
    geom_vline(xintercept= -fold_threshold) +
    geom_vline(xintercept= fold_threshold) +
    geom_hline(yintercept=-log10(p_threshold))+
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 12, hjust = 1),  # Rotate x-axis labels for clarity
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          panel.grid.major = element_line(linewidth = 0.4, color = "gray80"), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1, fill = NA)) 
  return(ggp) 
}

#function to plot MA plot####
plot_MA= function(de_table, p_threshold, fold_threshold) 
{ 
  de_table_sig = subset(de_table, p.adj < p_threshold & abs(log2fold) > fold_threshold)  
  de_table_sig_up = subset(de_table_sig, p.adj < p_threshold & (log2fold)> fold_threshold) 
  de_table_sig_down = subset(de_table_sig, p.adj < p_threshold & (log2fold) < (-fold_threshold)) 
  de_table_sig_up_top5 =   de_table_sig_up[1:5,] 
  de_table_sig_down_top5 =  de_table_sig_down[1:5,] 
  
  ggp1 = ggplot(de_table, aes(x=log10(mean_value), y=log2fold)) +  
    geom_point(aes(colour="a")) +  
    geom_point(data= de_table_sig_up, aes(colour="b")) +  
    geom_point(data= de_table_sig_down, aes(colour="c")) +
    labs(title="", x="-log10(mean_expression)", y="log2fold") +
    geom_label_repel(data=de_table_sig_up_top5, aes(label = symbols,colour ="b"), show.legend = FALSE) +
    geom_label_repel(data=de_table_sig_down_top5,aes(label = symbols,colour = "c"), show.legend = FALSE) + 
    scale_color_manual(values = c("black","red3","darkblue"),labels =c("No change","Up","Down"),name ="")+
    
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 12, hjust = 1),  # Rotate x-axis labels for clarity
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          panel.grid.major = element_line(linewidth = 0.4, color = "gray80"), 
          panel.grid.minor = element_blank(),
          panel.border = element_rect(size = 1, fill = NA))
  return (ggp1)
}


#function to create PCA plot####
plot_pca = function(expression_scaled, groups) 
  
{
  # Converet scaled expression data into a matrix 
  numeric_matrix = as.matrix(sapply(expression_scaled, as.numeric))
  pca = prcomp(t(numeric_matrix))
  pca_coordinates = data.frame(pca$x)
  pca_coordinates$sample = groups
  
  # To apply percentage variation 
  vars = apply(pca$x, 2, var) 
  prop_x = round(vars["PC1"] / sum(vars),4) * 100 
  prop_y = round(vars["PC2"] / sum(vars),4) * 100 
  x_axis_label = paste("PC1 ", " (",prop_x, "%)",sep="") 
  y_axis_label = paste("PC2 ", " (",prop_y, "%)",sep="") 
  
  # PLOT PCA 
  ggp = ggplot(pca_coordinates, aes(x = PC1, y = PC2, colour = groups)) + 
    geom_point(size = 5) +
    labs(x = x_axis_label, y = y_axis_label) +  # Add x and y axis labels
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
  return(ggp) 
}

#function to create density plot####
plot_exp_density = function(expression_table, columns) 
{
  # Melt expression table 
  expression_table_melted = melt(expression_table)
  
  # Plot Density 
  ggp = ggplot(expression_table_melted,aes(x = log10(value),fill = variable))+
    facet_wrap(~variable,ncol = columns)+
    xlab("Expression")+
    ylab("Density")+
    labs(title = "Density Plots")+
    geom_density() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.text.x = element_text(size = 12, hjust = 1), 
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"), 
      legend.position = "none",
      legend.title = element_blank(),
      legend.text = element_blank(),
      panel.grid.major = element_line(linewidth = 0.4, color = "gray80"), 
      panel.grid.minor = element_blank(),
      panel.border = element_rect(size = 1, fill = NA)) 
    
  return(ggp)
}


#function to create box plot for single gene####
plot_boxplot = function(gene_name,expression_tabe,group)
{
  # Create expression table of gene 
  gene = expression_tabe[gene_name,] 
  gene_data = data.frame(t(gene)) 
  gene_data$sample_group = group
  gene_data$sample_group = factor(gene_data.$sample_group , levels = level)
  names(gene_data) = c("expression","sample_group")
  
  # Plot the gene expression 
  ggp= ggplot(gene_data,aes(x=sample_group,y = expression,fill = sample_group)) +
    geom_boxplot(size = 1,  outlier.size = 0,  alpha = 0.5)  +
    geom_jitter(width = 0.2, size = 1, color = "black", alpha = 0.7) +
    labs(title = paste("Expression of", gene_name),
         x = "Sample Group",
         y = "Expression Level",
         fill = "Group") +
    theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
          axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          strip.text = element_text(size = 14, face = "bold"), 
          legend.title = element_text(size = 13),
          legend.text = element_text(size = 12),
          panel.grid.major = element_blank(),   
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),  
          axis.line = element_line(size = 0.5))
  return(ggp)
}

#function to create box plot for multiple genes####
plot_multi_boxplot = function(genes,expression_table,group)
{
  #only include 30 genes if list is more than 30 
  genes = head(genes, 30)
  gene_data = expression_table[genes,] 
  gene_data = data.frame(t(gene_data)) 
  gene_data$sample_group = factor(group) 
  
  gene_data.m = melt(gene_data, id.vars="sample_group",variable.name = "gene", value.name = "expression") 
  gene_data.m$sample_group = factor(gene_data.m$sample_group)
  
  
  
  ggp = ggplot(gene_data.m, aes(x=gene, y=expression, fill=sample_group)) +
    geom_boxplot(position=position_dodge(width=0.75), width=0.6) +
  
    labs(
      title = "",
      x = "Sample",
      y = "Expression Level",
      fill = "Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),   # Remove grid lines
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),  
      axis.line = element_line(size = 0.5)
    ) 
  
  return(ggp)
}

#Function to plot heat map with row or column clustering choice ####
plot_heatmap = function(data,row_cluster = TRUE,col_cluster = TRUE)
{
  hm.matrix = as.matrix(data)
  # If row clustering is true 
  if (row_cluster){
    y.dist = Dist(hm.matrix, method="spearman")
    y.cluster = hclust(y.dist, method="average")
    y.dd = as.dendrogram(y.cluster)
    y.dd.reorder = reorder(y.dd,0,FUN="average")
    y.order = order.dendrogram(y.dd.reorder)
    hm.matrix = hm.matrix[y.order,]
  }
  # If col clustering is true 
  if (col_cluster){
    x.dist = Dist(t(hm.matrix), method = "spearman")
    x.cluster = hclust(x.dist, method = "average")
    x.dd = as.dendrogram(x.cluster)
    x.dd.reorder = reorder(x.dd, 0, FUN = "average")
    x.order = order.dendrogram(x.dd.reorder)
    hm.matrix = hm.matrix[, x.order]
  }
  
  hm.matrix_clustered = melt(hm.matrix)
  colnames(hm.matrix_clustered) =  c("Var1", "Var2", "value")
  
  # Assign colours and palette
  colours = c("blue","pink" ,"red")
  palette = colorRampPalette(colours)(100)
  
  # Plot heat map 
  ggp =  ggplot(hm.matrix_clustered, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colours = palette,name="Expression Value") +
    ylab("Genes") +
    xlab("Samples") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_blank(),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1, face = "bold"),  # Improve X-axis label readability
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      legend.title = element_text(size = 12),  # Make the legend title more prominent
      legend.text = element_text(size = 12),  # Increase the size of the legend text
      legend.spacing.x = unit(0.25, 'cm'),  # Adjust the spacing of the legend
      
    )
  
  
  
  return(ggp)
}

#function for pathway ####
do_pathway = function(sig_genes, organism.db, ont_value, pvalue) {
  results_list = list()
  
  # Convert gene symbols to Entrez IDs
  sig_genes_entrez = bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism.db)
  sig_genes_entrez = na.omit(sig_genes_entrez)
 
  # Perform GO enrichment analysis
  ora_results = enrichGO(gene = sig_genes_entrez$ENTREZID, OrgDb = organism.db, readable = T, ont = ont_value, pvalueCutoff = pvalue, qvalueCutoff = 0.10)
  
  # Create bar plot
  ggp.bar = barplot(ora_results, showCategory = 10)
  
  # Create cnet plot
  ggp.cnet = cnetplot(ora_results, categorySize = "pvalue")
  
  
  #get the objects 
  gene_sets = ora_results$geneID 
  description = ora_results$Description 
  p.adj = ora_results$p.adjust
 
  
  
  # Create the results table
  ora_results_table = data.frame(gene_sets = ora_results$geneID, 
                                 p.adj = ora_results$p.adjust,
                                 row.names = ora_results$Description)
  
  
  
  # Return the list of results
  results_list$ggp.bar = ggp.bar
  results_list$ggp.cnet = ggp.cnet
  results_list$ora_results_table = ora_results_table
  
  return(results_list)
}

#function for plotting heatmap and boxplot for enriced genes in pathways####

enriched_heat_boxplot_string = function(ora_table, expression_table,group,pathway_number) {
  
  
  list = list()
  
  # Get enriched gene set from pathway table 
  enriched_gene_set = as.character(ora_table[pathway_number,"gene_sets"]) 
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))       
  
  # Create heatmap data 
  heatmap_data = expression_table[candidate_genes, ]
  heatmap_data = na.omit(heatmap_data)
  
  # Perform clustering by y 
  hm.matrix = as.matrix(heatmap_data)
  y.dist = Dist(hm.matrix, method="spearman")
  y.cluster = hclust(y.dist, method="average")
  y.dd = as.dendrogram(y.cluster)
  y.dd.reorder = reorder(y.dd,0,FUN="average")
  y.order = order.dendrogram(y.dd.reorder)
  hm.matrix = hm.matrix[y.order,]
  
  hm.matrix_clustered = melt(hm.matrix)
  colnames(hm.matrix_clustered) =  c("Var1", "Var2", "value")
  
  colours = c("blue","pink" ,"red")
  palette = colorRampPalette(colours)(100)
  
  # Plot heatmap 
  heat_map =  ggplot(hm.matrix_clustered, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colours = palette,name="Expression") +
    ylab("Samples") +
    xlab("Genes") +
    theme(
      plot.title = element_text(face = "bold"),
      axis.text.y = element_text(size= 10, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "bold"),  
      axis.title = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 12), 
      legend.text = element_text(size = 12),  
      legend.spacing.x = unit(0.25, 'cm'))
  
  
  # plot boxplot 
  boxplot_plot = plot_multi_boxplot(candidate_genes, expression_table,group)
  
  
  # Create string 
  
  #Create dataframe of genes from candidate genes 
  candidate_genes_table = data.frame(candidate_genes) 
  names(candidate_genes_table) = "gene"
  
  # load the database
  string_db = STRINGdb$new( version="11.5", species=9606, score_threshold=200, 
                            network_type="full", input_directory="") 
  
  # Map the candidate genes 
  string_mapped = string_db$map(candidate_genes_table, "gene", removeUnmappedRows = TRUE ) 
  
  # Plot the string netwek
  string_db$plot_network(string_mapped)
  
  
  # Display heat map and box plot 
  print(heat_map)
  print(boxplot_plot)
  
  list$string_mapped= string_mapped
  list$heat_map = heat_map
  list$boxplot_plot = boxplot_plot
  
  
  return(list)
}

# Function to make metagene boxplot####

make_metagene_boxplot = function (gene_list, expression_table, groups,signature_name)  
{ 
  metagene = data.frame(colMeans(expression_table[gene_list, ], na.rm = TRUE))
  colnames(metagene) = "metagene"
  
  # Convert 'metagene' into a vector for ggplot
  metagene_data = data.frame(metagene = as.vector(metagene), group = groups)
  
  # Ensure 'group' is a factor 
  metagene_data$group = factor(metagene_data$group , levels = unique(groups))
  
  # Plot metagene box plot
  ggp = ggplot(metagene_data, aes(x = group, y = metagene,fill = group)) +
    geom_boxplot() +
    ylab("gene Expression") +
    xlab("Group") +
    theme( plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
           axis.text.x = element_text(size = 12, angle = 45, hjust = 1),  # Rotate x-axis labels for clarity
           axis.text.y = element_text(size = 12),
           axis.title = element_text(size = 14, face = "bold"),
           strip.text = element_text(size = 14, face = "bold"), 
           legend.title = element_text(size = 13),
           legend.text = element_text(size = 12),
           panel.grid.major = element_line(linewidth = 0.4, color = "gray80"), 
           panel.grid.minor = element_blank(),
           panel.border = element_rect(size = 1, fill = NA))
    
  
  
  return(ggp)
} 

##function for signature analysis####

plot_signature = function(gene_list, expression_table, groups, organism_db, signature_name) {
  result_list = list()
  
  #Create and display the heatmap
  ggp_heatmap = plot_heatmap(expression_table[gene_list, ],row_cluster = TRUE, col_cluster = FALSE)
  print(ggp_heatmap)
  result_list$heatmap = ggp_heatmap
  
  #Create and display the metagene boxplot
  ggp_metagene = make_metagene_boxplot(gene_list, expression_table, groups, signature_name)
  result_list$metagene = ggp_metagene
  print(ggp_metagene)
  
  #Perform pathway analysis
  pathway_results = do_pathway(gene_list, organism_db, "BP", 0.05)
  # Store ORA table into results 
  result_list$ora_results_table = pathway_results$ora_results_table
  
  # Create and display the pathway bar plot
  ggp_bar = pathway_results$ggp.bar
  result_list$bar_plot = ggp_bar
  
  # Createand display the pathway c.net plot
  ggp_cnet = pathway_results$ggp.cnet
  result_list$cnet_plot = ggp_cnet
 
  
 
  
  return(result_list)
}

#some extra functions ####

#to do upregulated and downregulated pathway analysis
do_pathway_up_down = function(sig_genes_up, sig_genes_down, organism.db, ont_value, pvalue) {
  results_list = list()
  
  # Convert gene symbols to Entrez IDs for upregulated and downregulated genes
  upregulated_entrez = bitr(sig_genes_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism.db)
  downregulated_entrez = bitr(sig_genes_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = organism.db)
  
  # Perform GO enrichment analysis for upregulated genes
  ora_results_up = enrichGO(gene = upregulated_entrez$ENTREZID, OrgDb = organism.db, readable = T, 
                                       ont = ont_value, pvalueCutoff = pvalue, qvalueCutoff = 0.10)
  
  # Perform GO enrichment analysis for downregulated genes
  ora_results_down = enrichGO(gene = downregulated_entrez$ENTREZID, OrgDb = organism.db, readable = T, 
                                          ont = ont_value, pvalueCutoff = pvalue, qvalueCutoff = 0.10)
  
  # Create bar plots for both upregulated and downregulated pathways
  ggp.bar_up = barplot(ora_results_up, showCategory = 10)
  
  ggp.bar_down = barplot(ora_results_down, showCategory = 10)
  
  # Create cnet plots for both upregulated and downregulated pathways
  ggp.cnet_up = cnetplot(ora_results_up, categorySize = "pvalue")
  
  ggp.cnet_down = cnetplot(ora_results_down, categorySize = "pvalue")
 
  # Create the results tables for both upregulated and downregulated pathways
  ora_results_table_up = data.frame(gene_sets = ora_results_up$geneID, 
                                    p.adj = ora_results_up$p.adjust,
                                    row.names = ora_results_up$Description)
  ora_results_table_up = ora_results_table_up
  ora_results_table_down = data.frame(gene_sets = ora_results_down$geneID, 
                                      p.adj = ora_results_down$p.adjust,
                                      row.names = ora_results_down$Description)
  ora_results_table_down = ora_results_table_down
  
  # Return the list of results 
  results_list$ggp.bar_up = ggp.bar_up
  results_list$ggp.bar_down = ggp.bar_down
  results_list$ggp.cnet_up = ggp.cnet_up
  results_list$ggp.cnet_down = ggp.cnet_down
  results_list$ora_results_table_up = ora_results_table_up
  results_list$ora_results_table_down = ora_results_table_down
  
  return(results_list)
}

# Plot multi box plot with Facet wrap 
multi_boxplot_facet = function(genes,expression_table,group)
{
  #only include 30 genes if list is very big 
  genes = head(genes, 30)
  gene_data = expression_table[genes,] 
  gene_data = data.frame(t(gene_data)) 
  gene_data$sample_group = factor(group) 
  
  gene_data.m = melt(gene_data, id.vars="sample_group",variable.name = "gene", value.name = "expression") 
  gene_data.m$sample_group = factor(gene_data.m$sample_group)
  
  
  
  ggp = ggplot(gene_data.m, aes(x=sample_group, y=expression, fill=sample_group)) +
    facet_wrap(~gene, ncol = 5 , scales = "fixed")+
    geom_boxplot() +
    
    labs(
      title = "",
      x = "Sample",
      y = "Expression Level",
      fill = "Group"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      axis.text.x = element_text(size = 12, angle = 45, hjust = 1, face = "bold"),  
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold"), 
      legend.title = element_text(size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major = element_blank(),   
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),  
      axis.line = element_line(size = 0.5)) 
  
  return(ggp)
}


















