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
    pvalue_threshold = 0.05
    # 
    adj_pvalues = sapply(wall, function(x) p.adjust(x[,data_type], method="BH"))
    # extract the categories with a significant p-values (over or under represented pvalue)
    enriched_cat = which(apply(adj_pvalues <= pvalue_threshold & adj_pvalues > 0, 1, any))
    # 
    adj_pvalues[adj_pvalues > pvalue_threshold] = NA
    # extract the category names for the significant categories and merge them for all the comparisons
    if(cat_type == "GO"){
        to_extract = c("category","term","ontology")
    }else{
        to_extract = c("category")
    }
    cat = wall[[1]][enriched_cat, to_extract]
    # combine the significant categories for all comparisons
    adj_pvalues = adj_pvalues[enriched_cat,]
    if(is.matrix(adj_pvalues)){
        full_mat = cbind(cat, adj_pvalues)
    }else{
        full_mat = t(c(cat, adj_pvalues))
    }
    
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


plot_top_go = function(go, go_wall, ont, comp, top_nb){
    # extract the GO for the ontology
    ont_go = go[go$ontology == ont,c(2,4:6)]
    # conserve the top 20 for every comparison
    to_conserve = unique(c(sapply(comp, function(i) return(get_smallest_value_id(ont_go, i, top_nb)))))
    ont_go = ont_go[to_conserve,]
    ont_go = ont_go[order(ont_go[,comp[1]]),]
    # extract the ratios of the conserved GO
    ratios = sapply(go_wall, function(x) x[to_conserve,"numDEInCat"]/x[to_conserve,"numInCat"])
    # create a matrix for the ggplot
    go_mat = expand.grid(x = ont_go$term, y = comp)
    colnames(go_mat) = c("term", "comparison")
    go_mat$p_value = c(sapply(comp, function(i) return(ont_go[, i])))
    go_mat$ratio = c(sapply(comp, function(i) return(ratios[, i])))
    # plot                 
    require(ggplot2)
    p = ggplot(go_mat, aes(factor(comparison), factor(term))) + labs(x = "", y = "")
    p + geom_point(aes(size=ratio,col=p_value)) + scale_colour_gradient(low = "blue", high="red")
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

get_interesting_GO = function(go, ont){
    col_nb = dim(go)[2]
    not_na = !is.na(go[,4:col_nb])
    if(!is.matrix(not_na)) not_na = t(not_na)
    return(go[go[,"ontology"] == ont & apply(not_na, 1, any), "category"])
}


create_GO_network = function(deg, ont, ont_go){
    l = list()
    # Get GO term for ontology and that are at least once over-represented or under-represented
    l$interesting_GO = c(get_interesting_GO(deg$over_represented_GO, ont),
                       get_interesting_GO(deg$under_represented_GO, ont))
    # Get similarity between GO terms
    l$go_sim = mgoSim(l$interesting_GO, l$interesting_GO, semData=ont_go, measure="Wang", combine=NULL)
    # Extract adjency matrix by keeping only the distance > 0.5
    l$go_sim = (l$go_sim > 0.5)*1
    # Replace the diagonal values by 0
    diag(l$go_sim) = 0
    # Keep only the GO with at least one connection
    #l$go_sim = l$go_sim[rowSums(l$go_sim)>0,rowSums(l$go_sim)>0]
    # Build the network
    l$go_net = graph_from_adjacency_matrix(l$go_sim, diag = FALSE, weighted = TRUE, mode="undirected")
    l$go_net_layout = layout_with_fr(l$go_net)
    l$go_net_layout_3d = layout_with_fr(l$go_net, dim=3)
    #plot(l$go_net, 
    #     vertex.label=NA,
    #     vertex.size=4,
    #     layout=l$go_net_layout)
    # Update list of interesting GO
    l$interesting_GO = rownames(l$go_sim)
    return(l)
}


get_ont_GO = function(go, ont){
    col_nb = dim(go)[2]
    return(go[go$ontology == ont, c(1,4:col_nb)])
}

get_col_ramp = function(values, min_col, max_col){
    if(dim(values)[1] == 0){
        return(c())
    }else{
        all_values = unlist(c(values[,2:dim(values)[2]]))
        all_values = all_values[!is.na(all_values)]
        ii = cut(all_values,
            breaks = exp(log(10)*seq(log10(min(all_values)),log10(0.05),len = 100)),
            include.lowest = TRUE)
        colors = cbind(sort(all_values), colorRampPalette(c(min_col, max_col))(99)[ii])
        return(colors)
    }
}

get_GO_network_col = function(over_represented, under_represented, comp, all_GO){
    #
    over_repr_colors = get_col_ramp(over_represented, "red", "lightpink")
    under_repr_colors = get_col_ramp(under_represented, "blue", "lightblue")
    # extract under and over represented GO
    over_represented_GO = over_represented[!is.na(over_represented[,comp]),1]
    under_represented_GO = under_represented[!is.na(under_represented[,comp]),1]
    # colors for the nodes
    col = rep("white", length(all_GO))
    names(col) = all_GO
    if(length(over_represented_GO) > 0){
        val = over_represented[names(col) %in% over_represented_GO,comp]
        col[names(col) %in% over_represented_GO] = sapply(val, function(i) return(over_repr_colors[which(over_repr_colors[,1] == i)[1],2]))
    }
    if(length(under_represented_GO) > 0){
        val = under_represented_GO[names(col) %in% under_represented_GO,comp]
        col[names(col) %in% under_represented_GO] = sapply(val, function(i) return(under_repr_colors[which(under_repr_colors[,1] == i)[1],2]))
    }
    col = unlist(col)
    return(col)
}


plot_GO_network = function(over_repr, under_repr, network, comp){
    # get color based on adjusted p-value
    col = get_GO_network_col(over_repr, under_repr, comp, network$interesting_GO)
    # plot network
    plot(network$go_net, 
         vertex.label=NA,
         vertex.size=4,
         vertex.color=col,
         layout=network$go_net_layout)
}

plot_interactive_GO_network = function(over_repr, under_repr, net, comp, full_go_desc){
    col = get_GO_network_col(over_repr, under_repr, comp, net$interesting_GO)
    # get term for selected GO terms
    go_desc = full_go_desc[names(col)]
    # plot the graph
    graphjs(net$go_net,
            layout=net$go_net_layout_3d,
            height=500,
            vertex.size=1,
            vertex.color=as.vector(col),
            vertex.shape="sphere",
            vertex.label=as.vector(go_desc),
            edge.color="grey")
}
