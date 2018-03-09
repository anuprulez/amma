get_perc = function(v, l){
    s = sum(v, na.rm=T)
    return(round(s*100/l, digits=2))
}

get_sign_padj = function(dge_res) {
    return(get_perc(dge_res$padj < 0.05, sum(!is.na(dge_res$log2FoldChange))))
}

get_pos_sign_padj = function(dge_res) {
    return(get_perc(dge_res$padj < 0.05 & dge_res$log2FoldChange > 0, sum(!is.na(dge_res$log2FoldChange))))
}

get_neg_sign_padj = function(dge_res) {
    return(get_perc(dge_res$padj < 0.05 & dge_res$log2FoldChange < 0, sum(!is.na(dge_res$log2FoldChange))))
}

get_stats_padj  = function(dge_res) {
    v = c(get_sign_padj(dge_res), get_pos_sign_padj(dge_res), get_neg_sign_padj(dge_res))
    names(v) = c("Wald padj < 0.05", "LFC > 0 (Wald padj < 0.05)", "LFC < 0 (Wald padj < 0.05)")
    return(v)
}

clean_mat = function(mat){
    new_mat = mat
    new_mat[is.na(new_mat)] = 0
    new_mat = new_mat[rowSums(new_mat) > 0,]
    return(new_mat)
}

get_interesting_cat = function(wall, data_type, cat_type){
    # extract the categories with a significant p-values (over or under represented pvalue)
    enriched_cat = lapply(wall, function(x) x[p.adjust(x[,data_type],method="BH")<.05,])
    # extract the category names for the significant categories and merge them for all the comparisons
    if(cat_type == "GO"){
        cat = unique(do.call(rbind,lapply(enriched_cat, function(x) x[,c("category","term","ontology")])))
        rownames(cat) = do.call(rbind,strsplit(rownames(cat), split='.', fixed=TRUE))[,2]
    }else{
        cat = unique(do.call(rbind,lapply(enriched_cat, function(x) cbind(rownames(x),x[,c("category")]))))
        names = cat[,1]
        cat = matrix(cat[,-1])
        rownames(cat) = names
    }
    # combine the significant categories for all comparisons
    full_mat = cbind(cat, sapply(enriched_cat, function(x) x[rownames(cat),data_type]))
    return(full_mat)
}

extract_diff_expr_genes = function(in_l, name){
    l = list()
    # extract the significant differentially expressed genes (all, upregulated, downregulated)
    l$deg = sapply(in_l, function(mat) return(mat$padj < 0.05))*1
    rownames(l$deg) = rownames(in_l[[1]])
    l$deg = clean_mat(l$deg)
    deg_names = rownames(l$deg)
    head(l$deg)
    # extract the significant differentially more expressed genes
    l$pos = sapply(in_l, function(mat) return(mat[deg_names, 'log2FoldChange'] > 0))*1
    l$pos[l$deg == 0] = 0
    rownames(l$pos) = deg_names
    l$pos = clean_mat(l$pos)
    # extract the significant differentially less expressed genes
    l$neg = sapply(in_l, function(mat) return(mat[deg_names, 'log2FoldChange'] < 0))*1
    l$neg[l$deg == 0] = 0
    rownames(l$neg) = deg_names
    l$neg = clean_mat(l$neg)
    # extract the log2FC of the significant differentially expressed genes
    l$fc_deg = sapply(in_l, function(mat) return(mat[deg_names, 'log2FoldChange']))
    rownames(l$fc_deg) = deg_names
    l$fc_deg[l$deg == 0] = NA
    write.table(l$fc_deg, paste("../results/dge/", name), sep = "\t", quote = FALSE)
    #system(paste("put -p", name, "-t tabular"), intern=T)
    ## GO and KEGG analysis                  
    assayed_genes = rownames(in_l[[1]])
    # extract the DE genes (1 in the deg matrix)
    de_genes = sapply(colnames(l$deg), function(x) names(which(abs(l$fc_deg[,x])>=1)))
    # extract in the full list of genes the DE ones                  
    gene_vector = sapply(de_genes, function(x) as.integer(assayed_genes%in%x))
    rownames(gene_vector) = assayed_genes
    # fit the probability weighting function
    pwf = lapply(1:dim(gene_vector)[2], function(x) nullp(gene_vector[,x], 'mm10', 'geneSymbol', plot.fit=F))
    names(pwf) = colnames(gene_vector)
    # calculate the over and under expressed GO categories among the DE genes
    l$GO_wall = lapply(pwf, function(x) goseq(x,'mm10', 'geneSymbol'))
    # calculate the over and under expressed KEGG pathways among the DE genes
    l$KEGG_wall = lapply(pwf, function(x) goseq(x,'mm10', 'geneSymbol', test.cats="KEGG"))
    # extract interesting pathways/categories and export them
    l$over_represented_GO = get_interesting_cat(l$GO_wall, "over_represented_pvalue", "GO")
    write.table(l$over_represented_GO, paste("../results/dge/", name, "over_represented_GO"), sep = "\t", quote = FALSE)                
    l$under_represented_GO = get_interesting_cat(l$GO_wall, "under_represented_pvalue", "GO")
    write.table(l$under_represented_GO, paste("../results/dge/", name, "under_represented_GO"), sep = "\t", quote = FALSE)
    l$over_represented_KEGG = get_interesting_cat(l$KEGG_wall, "over_represented_pvalue", "KEGG")
    write.table(l$over_represented_KEGG, paste("../results/dge/", name, "over_represented_KEGG"), sep = "\t", quote = FALSE)    
    l$under_represented_KEGG = get_interesting_cat(l$KEGG_wall, "under_represented_pvalue", "KEGG")
    write.table(l$under_represented_KEGG, paste("../results/dge/", name, "under_represented_KEGG"), sep = "\t", quote = FALSE)                       
    return(l)
}

get_number = function(bool_expr) return(sum(bool_expr*1))
get_repartition = function(dge_mat) return(cbind(get_number(dge_mat[,1]>0 & dge_mat[,2]>0), get_number(dge_mat[,1]>0 & dge_mat[,2]<0), get_number(dge_mat[,1]<0 & dge_mat[,2]>0), get_number(dge_mat[,1]<0 & dge_mat[,2]<0)))
get_repartition_3col = function(dge_mat) return(cbind(
    get_number(dge_mat[,1]>0 & dge_mat[,2]>0 & dge_mat[,3]>0),
    get_number(dge_mat[,1]>0 & dge_mat[,2]>0 & dge_mat[,3]<0), 
    get_number(dge_mat[,1]>0 & dge_mat[,2]<0 & dge_mat[,3]>0), 
    get_number(dge_mat[,1]>0 & dge_mat[,2]<0 & dge_mat[,3]<0), 
    get_number(dge_mat[,1]<0 & dge_mat[,2]>0 & dge_mat[,3]>0), 
    get_number(dge_mat[,1]<0 & dge_mat[,2]>0 & dge_mat[,3]<0), 
    get_number(dge_mat[,1]<0 & dge_mat[,2]<0 & dge_mat[,3]>0), 
    get_number(dge_mat[,1]<0 & dge_mat[,2]<0 & dge_mat[,3]<0)))

quantile_breaks = function(xs, n = 10) {
    breaks = quantile(xs, probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
}
plot_count_heatmap = function(genes, samples, annot){
    data = norm_counts[genes,samples]
    breaks = quantile_breaks(data, n = 11)
    pheatmap(data, cluster_rows=F, cluster_cols=F, show_rownames=F, show_colnames=F, annotation_col=annot, breaks=breaks, color=inferno(10))  
}

get_list = function(mapping){
    mapped_genes = mappedkeys(mapping)
    return(as.list(mapping[mapped_genes]))
}

# search on Entrez the correct name
search_name = function(name){
    search = entrez_search(db="gene",term=name)
    names = c()
    for(id in search$ids){
        sum = entrez_summary(db="gene", id=id)
        if(sum$organism$scientificname == 'Mus musculus' & grepl(name,sum$otheraliases)){
            names = c(names, sum$name)
        }
    }
    return(names)
}

capFirst = function(s) {
    paste(substring(s, 1, 1), tolower(substring(s, 2)), sep = "")
}


plot_net_with_fr_layout = function(net, colors, pal2){
    l = layout_with_fr(net)
    plot(net, 
        vertex.label=NA,
        vertex.size=4,
        vertex.color=pal2[colors],
        layout=l)
    legend(x=-1.5,
        y=-1,
        c("response to endoplasmic reticulum stress/transport", "chromatine/chromosome organization/regulation", "metabolic process", "translation","immune system"),
        pch=21,
        col="#777777",
        pt.bg=pal2[c(1,2,3,5,9)],
        pt.cex=2,
        cex=.8,
        bty="n",
        ncol=1)
}

get_deg_colors = function(comp_deg, comp, connected_gene_colors, module_nb){
    deg_col = connected_gene_colors
    sign_pos_FC = rownames(comp_deg$pos)[comp_deg$pos[,comp] == 1]
    deg_col[sign_pos_FC[sign_pos_FC %in% names(deg_col)]] = module_nb + 1
    sign_neg_FC = rownames(comp_deg$neg)[comp_deg$neg[,comp] == 1]
    deg_col[sign_neg_FC[sign_neg_FC %in% names(deg_col)]] = module_nb + 2
    return(deg_col)
}

get_smallest_value_id = function(matrix, column, nb){
    return(rownames(matrix)[order(matrix[,column])][1:nb])
}
