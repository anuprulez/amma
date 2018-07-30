get_perc = function(v, l){
    s = sum(v, na.rm=T)
    return(s)
    #return(round(s*100/l, digits=2))
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

plot_stat_mat = function(stat_mat){
    #colnames(stat_mat) = c("DEG", "Up-regulated", "Down-regulated")
    mat = melt(stat_mat)
    colnames(mat) = c("comp", "type", "value")
    mat$comp = factor(mat$comp)
    mat$type = factor(mat$type)
    ggplot(mat, aes(x = reorder(comp, desc(comp)), y = reorder(type, desc(type)))) +
        labs(x = "", y = "") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.title.y = element_text(size = rel(1.8)), panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_point(aes(size=value,col=value)) +
        scale_colour_gradient(low = "blue", high="red")
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
    if(!is.matrix(adj_pvalues)){
        if(length(adj_pvalues)==0){
            full_mat = c()
        }
        else if(length(names(wall))==1){
            full_mat = cbind(cat, adj_pvalues)
            colnames(full_mat) = c(to_extract, names(wall))
        }else{
            full_mat = t(c(cat, adj_pvalues))
            colnames(full_mat) = c(to_extract, names(wall))
        }
    }else{
        full_mat = cbind(cat, adj_pvalues)
        colnames(full_mat) = c(to_extract, names(wall)) 
    }
    return(full_mat)
}


get_cat_de_genes = function(cat_id, cat2genes, de_genes){
    cat_genes = cat2genes[cat2genes$values == cat_id, "ind"]
    return(as.character(cat_genes[cat_genes %in% toupper(de_genes)]))
}

extract_cat_de_genes = function(selected_cat, cat, sign_fc_deg, cat2genes, file_prefix){
    de_genes = rownames(sign_fc_deg)[!is.na(sign_fc_deg[,cat])]
    cat_de_genes = sapply(selected_cat, function(y) get_cat_de_genes(y, cat2genes, de_genes))
    capture.output(cat_de_genes, file = paste(file_prefix, gsub("[(),]", "", gsub(" ", "_", cat)), sep=""))
}

get_interesting_categories = function(deg, interesting_cat){
    l = sapply(colnames(deg), function(x) interesting_cat[which(!is.na(interesting_cat[,x])),"category"])
    if(dim(deg)[2]>1){
        l = as.list(l)
    }else{
        l2 = list(l)
        l2[[colnames(deg)]] = l
        l = l2
    }
    names(l) = gsub(".category", "", names(l))
    return(l)
}

extract_diff_expr_genes = function(in_l, dir_path){
    # create dir if it does not exist
    full_dir_path = paste("../results/dge/", dir_path, sep="")
    dir.create(full_dir_path, showWarnings = FALSE)
    go_dir_path = paste(full_dir_path, "go/", sep="")
    dir.create(go_dir_path, showWarnings = FALSE)
    kegg_dir_path = paste(full_dir_path, "kegg/", sep="")
    dir.create(kegg_dir_path, showWarnings = FALSE)
    l = list()
    # extract the significant differentially expressed genes (all, upregulated, downregulated)
    l$deg = sapply(in_l, function(mat) return(mat$padj < 0.05))*1
    rownames(l$deg) = rownames(in_l[[1]])
    l$deg = clean_mat(l$deg)
    l$deg = as.data.frame(l$deg)
    colnames(l$deg) = names(in_l)
    deg_names=rownames(l$deg)
    # extract the log2FC of the significant differentially expressed genes
    l$fc_deg = sapply(in_l, function(mat) return(mat[deg_names, 'log2FoldChange']))
    rownames(l$fc_deg) = deg_names
    l$fc_deg[l$deg == 0] = NA
    write.table(l$fc_deg, paste(full_dir_path, "fc_deg", sep=""), sep = "\t", quote = FALSE)
    # extract the sign log2FC (> log2(1.5))
    l$sign_fc_deg = l$fc_deg
    l$sign_fc_deg[abs(l$sign_fc_deg) < log2(1.5)] = NA
    l$sign_fc_deg = l$sign_fc_deg[!apply(is.na(l$sign_fc_deg),1,all),]
    l$sign_fc_deg = as.data.frame(l$sign_fc_deg)
    colnames(l$sign_fc_deg) = names(in_l)
    write.table(l$sign_fc_deg, paste(full_dir_path, "sign_fc_deg", sep=""), sep = "\t", quote = FALSE)
    # extract stats 
    deg_stat = colSums(!is.na(l$fc_deg), na.rm = TRUE)
    pos_deg_stat = colSums(l$fc_deg > 0, na.rm = TRUE)
    neg_deg_stat = colSums(l$fc_deg < 0, na.rm = TRUE)
    sign_deg_stat = colSums(!is.na(l$sign_fc_deg), na.rm = TRUE)
    sign_pos_deg_stat = colSums(l$sign_fc_deg > 0, na.rm = TRUE)
    sign_neg_deg_stat = colSums(l$sign_fc_deg < 0, na.rm = TRUE)
    l$stat = cbind(deg_stat, pos_deg_stat, neg_deg_stat, sign_deg_stat, sign_pos_deg_stat, sign_neg_deg_stat)
    colnames(l$stat) = c("All DEG (Wald padj < 0.05)",
                         "All over-expressed genes (Wald padj < 0.05 & FC > 0)",
                         "All under-expressed genes (Wald padj < 0.05 & FC < 0)",
                         "DEG (Wald padj < 0.05 & abs(FC) >= 1.5)",
                         "Over-expressed genes (Wald padj < 0.05 & FC >= 1.5)",
                         "Under-expressed genes (Wald padj < 0.05 & FC <= -1.5)")
    plot_stat_mat(l$stat)
    ## GO and KEGG analysis
    # get gene vector                  
    gene_vector = matrix(0, ncol = length(colnames(l$sign_fc_deg)), nrow = length(rownames(in_l[[1]])))
    rownames(gene_vector) = rownames(in_l[[1]])
    colnames(gene_vector) = colnames(l$sign_fc_deg)
    assayed_genes_in_sign_genes = rownames(gene_vector)[rownames(gene_vector) %in% rownames(l$sign_fc_deg)]
    gene_vector[assayed_genes_in_sign_genes,] = 1*(!is.na(l$sign_fc_deg[assayed_genes_in_sign_genes,]))
    # fit the probability weighting function
    pwf = lapply(1:dim(gene_vector)[2], function(x) suppressMessages(nullp(gene_vector[,x], 'mm10', 'geneSymbol', plot.fit=F,  bias.data=gene_length)))
    names(pwf) = colnames(gene_vector)
    ### GO analysis
    # calculate the over and under expressed GO categories among the DE genes
    l$GO_wall = lapply(pwf, function(x) suppressMessages(goseq(x, 'mm10', 'geneSymbol')))
    # extract interesting pathways/categories and export them
    l$over_represented_GO = get_interesting_cat(l$GO_wall, "over_represented_pvalue", "GO")
    write.table(l$over_represented_GO, paste(go_dir_path, "full_over_represented_GO", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    l$under_represented_GO = get_interesting_cat(l$GO_wall, "under_represented_pvalue", "GO")
    write.table(l$under_represented_GO, paste(go_dir_path, "full_under_represented_GO", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    # exclude GO terms
    go_term_to_exclude = read.table("../data/go_term_to_exclude.csv", h = F)$V1
    over_represented_GO_to_keep = !l$over_represented_GO[, "category"] %in% go_term_to_exclude
    col_names = colnames(l$over_represented_GO)
    l$over_represented_GO = as.data.frame(l$over_represented_GO[over_represented_GO_to_keep,])
    colnames(l$over_represented_GO) = col_names                
    write.table(l$over_represented_GO, paste(go_dir_path, "over_represented_GO", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    col_names = colnames(l$under_represented_GO)
    under_represented_GO_to_keep = !l$under_represented_GO[, "category"] %in% go_term_to_exclude               
    l$under_represented_GO = as.data.frame(l$under_represented_GO[under_represented_GO_to_keep,])
    colnames(l$under_represented_GO) = col_names                    
    write.table(l$under_represented_GO, paste(go_dir_path, "under_represented_GO", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    # extract list of genes involved in the over and under represented GO
    full_go_genes = stack(getgo(rownames(l$sign_fc_deg), 'mm10', 'geneSymbol'))
    over_represented_GO = get_interesting_categories(l$sign_fc_deg, l$over_represented_GO)
    under_represented_GO = get_interesting_categories(l$sign_fc_deg, l$under_represented_GO)
    for(x in colnames(l$sign_fc_deg)){
        extract_cat_de_genes(over_represented_GO[[x]], x, l$sign_fc_deg, full_go_genes, paste(go_dir_path, "over_repr_", sep = ""))
        extract_cat_de_genes(under_represented_GO[[x]], x, l$sign_fc_deg, full_go_genes, paste(go_dir_path, "under_repr_", sep = ""))
    }
    # KEGG analysis                   
    # calculate the over and under expressed KEGG pathways among the DE genes
    l$KEGG_wall = lapply(pwf, function(x) suppressMessages(goseq(x, 'mm10', 'geneSymbol', test.cats="KEGG")))
    # extract interesting pathways/categories and export them
    l$over_represented_KEGG = get_interesting_cat(l$KEGG_wall, "over_represented_pvalue", "KEGG")
    write.table(l$over_represented_KEGG, paste(kegg_dir_path, "over_represented_KEGG", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)    
    l$under_represented_KEGG = get_interesting_cat(l$KEGG_wall, "under_represented_pvalue", "KEGG")
    write.table(l$under_represented_KEGG, paste(kegg_dir_path, "under_represented_KEGG", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
    # extract list of genes involved in the over and under represented KEGG
    full_kegg_genes = stack(getgo(rownames(as.data.frame(l$deg)), 'mm10', 'geneSymbol', fetch.cats="KEGG"))
    over_represented_KEGG = get_interesting_categories(l$deg, l$over_represented_KEGG)
    under_represented_KEGG =  get_interesting_categories(l$deg, l$under_represented_KEGG)
    for(x in colnames(l$deg)){
        extract_cat_de_genes(over_represented_KEGG[[x]], x, l$sign_fc_deg, full_kegg_genes, paste(kegg_dir_path, "over_repr_", sep = ""))
        extract_cat_de_genes(under_represented_KEGG[[x]], x, l$sign_fc_deg, full_kegg_genes, paste(kegg_dir_path, "under_repr_", sep = ""))
    }
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
    breaks = quantile(unlist(xs), probs = seq(0, 1, length.out = n))
    breaks[!duplicated(breaks)]
}

plot_count_heatmap = function(genes, samples, annot){
    plot_heatmap(norm_counts, genes, samples, annot)
}

plot_heatmap = function(count, genes, samples, annot){
    data = count[genes,samples]
    breaks = quantile_breaks(data, n = 11)
    pheatmap(data,
             cluster_rows=F,
             cluster_cols=F,
             show_rownames=F,
             show_colnames=F,
             annotation_col=annot,
             breaks=breaks,
             color=inferno(10))  
}

plot_fc_heatmap = function(sign_fc_deg, annot){
    data_vector = unlist(sign_fc_deg)
    data_vector = data_vector[!is.na(data_vector)]
    breaks = quantile_breaks(data_vector, n = 11)
    data = sign_fc_deg
    data[is.na(data)] = 0
    pheatmap(data,
        cluster_rows=T,
        cluster_cols=F,
        show_rownames=F,
        show_colnames=F,
        annotation_col=annot,
        breaks=breaks,
        color=rev(brewer.pal(10, "RdBu")))
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


plot_net_with_layout = function(net, colors, pal2, layout, add_legend = TRUE){
    plot(net, 
        vertex.label=NA,
        vertex.size=4,
        vertex.color=pal2[colors],
        layout=layout)
    if(add_legend){
       module_annot = list()
        module_annot$"1" = "chromatine/chromosome organization/RNA metabolic process"
        module_annot$"2" = "response to endoplasmic reticulum stress/transport"
        module_annot$"3" = "metabolic process (ATP, ribonucleoside)"
        module_annot$"6" = "metabolic process (primary, cellular)"
        module_annot$"7" = "immune system"
        module_annot$"8" = "translation/rRNA"
        module_annot$"10" = "organelle"
        module_annot$"11" = "localization"
        legend(x=-1.5,
            y=-.9,
            unlist(module_annot),
            pch=21,
            col="#777777",
            pt.bg=pal2[as.integer(names(module_annot))],
            pt.cex=2,
            cex=.8,
            bty="n",
            ncol=1) 
    }
}

get_top_go = function(go, top_nb, ont, comp){
    ont_go = go[go$ontology == ont,c(2,4:dim(go)[2])]
    if(length(ont_go)==0){
        return(ont_go)
    }else{
        # conserve the top 20 for every comparison
        to_conserve = unique(c(sapply(comp, function(i) return(get_smallest_value_id(ont_go, i, top_nb)))))
        to_conserve = to_conserve[!is.na(to_conserve)]
        return(ont_go[to_conserve,])
    }
}


plot_top_go = function(deg, ont, top_nb){
    over_repr_go = deg$over_represented_GO
    under_repr_go = deg$under_represented_GO
    go_wall = deg$GO_wall
    # extract the top GO for the ontology                                  
    over_repr_top_go = get_top_go(over_repr_go, top_nb, ont, names(go_wall))
    over_repr_top_go$id = rownames(over_repr_top_go)
    under_repr_top_go = get_top_go(under_repr_go, top_nb, ont, names(go_wall))
    under_repr_top_go$id = rownames(under_repr_top_go)
    if(dim(over_repr_top_go)[1] > 0 || dim(under_repr_top_go)[1] > 0){
        # extract the ratios of the conserved GO
        cons_GO = c(over_repr_top_go$id, under_repr_top_go$id)
        ratios = data.frame(sapply(go_wall, function(x) x[cons_GO,"numDEInCat"]/x[cons_GO,"numInCat"]))
        rownames(ratios) = cons_GO
        # create a matrix for the ggplot
        if(dim(over_repr_top_go)[1] > 0){
            over_mat = cbind(melt(over_repr_top_go), "over")
            colnames(over_mat) = c("term", "id", "comparison", "p_value", "type")
        }else{
            over_mat = NULL
        }
        if(dim(under_repr_top_go)[1] > 0){
            under_mat = cbind(melt(under_repr_top_go), "under")
            colnames(under_mat) = c("term", "id", "comparison", "p_value", "type")
        }else{
            under_mat = NULL
        }
        mat = rbind(over_mat, under_mat)
        # add ratios
        mat$ratio = sapply(1:dim(mat)[1], function(i) ratios[mat$id[i], mat$comp[i]])
        # reformate the terms to have less than 8 words
        mat$term = sapply(strsplit(as.character(mat$term), " "), function(i) paste(i[1:ifelse(length(i)>8,8,length(i))],collapse = " "))                
        ## plot
        size_scale_lim = c(min(mat$ratio),max(mat$ratio))
        col_scale_lim = c(min(mat$p_value, na.rm = T),max(mat$p_value, na.rm = T))
        plot1 <- mat %>%
          filter(type == "over") %>%
          ggplot() +
          geom_point(aes(x = factor(comparison), y = factor(term), size=ratio,col=p_value)) +
          labs(x = "", y = "Over represented") +
          theme(axis.text.x = element_blank(), axis.text.y = element_text(size = rel(.75))) +
          scale_colour_gradient(limits = col_scale_lim, low = "red", high="blue", guide = 'none') +
          scale_size(limits = size_scale_lim, range=c(0.2,2))
        plot2 <- mat %>%
          filter(type == "under") %>%
          ggplot() +
          geom_point(aes(x = factor(comparison), y = factor(term), size=ratio,col=p_value)) +
          labs(x = "", y = "Under represented") +
          theme(axis.text.x = element_text(angle = 90, hjust = 1), axis.text.y = element_text(size = rel(.75))) +
          scale_colour_gradient(limits = col_scale_lim, low = "red", high="blue") +
          scale_size(limits = size_scale_lim, range=c(0.2,2), guide = 'none')
        grid.newpage()
        #grid.arrange(arrangeGrob(plot1, plot2, nrow=2), heights=c(4,1))
        #print(plot_both)
        grid.draw(rbind(ggplotGrob(plot1), ggplotGrob(plot2), size = "first"))
    }
}

get_deg_colors = function(comp_deg, comp, connected_gene_colors, module_nb){
    deg_col = connected_gene_colors
    deg_col = connected_gene_colors
    sign_pos_FC = rownames(comp_deg$sign_fc_deg)[comp_deg$sign_fc_deg[,comp] > 0]
    sign_neg_FC = rownames(comp_deg$sign_fc_deg)[comp_deg$sign_fc_deg[,comp] < 0]
    deg_col[which(names(deg_col) %in% sign_pos_FC)] = module_nb + 1
    deg_col[which(names(deg_col) %in% sign_neg_FC)] = module_nb + 2
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
    # extract over and under represented GO
    l$over_repr = get_ont_GO(deg$over_represented_GO, ont)
    l$under_repr = get_ont_GO(deg$under_represented_GO, ont)
    # 
    return(l)
}


get_ont_GO = function(go, ont){
    col_nb = dim(go)[2]
    val = which(go$ontology == ont)
    ont_go = go[val, 4:col_nb]
    if(!is.matrix(ont_go) && !is.data.frame(ont_go)){
        ont_go = matrix(ont_go, ncol = 1)
        colnames(ont_go) = colnames(go)[4:col_nb]
    }
    rownames(ont_go) = go[val, 1]
    return(ont_go)
}

get_col_ramp = function(values, min_col, max_col){
    if(dim(values)[1] == 0){
        return(c())
    }else{
        all_values = unlist(values)
        all_values = all_values[!is.na(all_values)]
        ii = cut(all_values,
            breaks = exp(log(10)*seq(log10(min(all_values)),log10(0.05),len = 100)),
            include.lowest = TRUE)
        colors = cbind(sort(all_values), colorRampPalette(c(min_col, max_col))(99)[ii])
        return(colors)
    }
}

get_GO_network_col = function(net, comp, all_GO = NULL){
    #
    over_repr_colors = get_col_ramp(net$over_repr, "red", "lightpink")
    under_repr_colors = get_col_ramp(net$under_repr, "blue", "lightblue")
    # extract under and over represented GO and their value
    v = which(!is.na(net$over_repr[,comp]))
    over_repr_GO = net$over_repr[v,comp]
    names(over_repr_GO) = rownames(net$over_repr)[v]
    v = which(!is.na(net$under_repr[,comp]))
    under_repr_GO = net$under_repr[v,comp]
    names(under_repr_GO) = rownames(net$under_repr)[v]
    # colors for the nodes
    col = rep("white", length(net$interesting_GO))
    names(col) = net$interesting_GO
    if(length(over_repr_GO) > 0){
        col[names(over_repr_GO)] = sapply(over_repr_GO, function(i) return(over_repr_colors[which(over_repr_colors[,1] == i)[1],2]))
    }
    if(length(under_repr_GO) > 0){
        col[names(under_repr_GO)] = sapply(under_repr_GO, function(i) return(under_repr_colors[which(under_repr_colors[,1] == i)[1],2]))
    } 
    if(!is.null(all_GO)){
        names(col) = all_GO[net$interesting_GO]
    }
    return(col)
}


plot_GO_network = function(net, col){
    plot(net$go_net, 
         vertex.label=NA,
         vertex.size=4,
         vertex.color=col,
         layout=net$go_net_layout)
}


plot_interactive_GO_network = function(net, col, go_desc){
    graphjs(net$go_net,
            layout=net$go_net_layout_3d,
            height=500,
            vertex.size=1,
            vertex.color=as.vector(col),
            vertex.shape="sphere",
            vertex.label=as.vector(go_desc),
            edge.color="grey")
}

plot_GO_networks = function(net, comp, full_go_desc, plot_non_interactive = TRUE, plot_interactive = TRUE){
    # get color based on adjusted p-value
    col = get_GO_network_col(net, comp, full_go_desc)
    # get names for the vertext: term for selected GO terms
    go_desc = full_go_desc[names(col)]
    # plots
    if(plot_non_interactive){
        plot_GO_network(net, col)
    }
    if(plot_interactive){
        plot_interactive_GO_network(net, col, names(col))
    }
}


plot_kegg_pathways = function(kegg_cats, fc_deg, dir_path){
    dir.create(dir_path, showWarnings = FALSE)
    for(cat in kegg_cats){
        if(cat!="01100"){
            suppressMessages(pathview(gene.data=fc_deg,
                                      pathway.id=cat,
                                      species="Mus musculus",
                                      gene.idtype="Symbol"))
            if(dim(fc_deg)[2] > 1){
                file.rename(from=paste("mmu", cat, ".pathview.multi.png",sep=""),
                            to=paste(dir_path,"mmu", cat, ".pathview.multi.png",sep=""))
            }else{
                file.rename(from=paste("mmu", cat, ".pathview.png",sep=""),
                            to=paste(dir_path,"mmu", cat, ".pathview.png",sep=""))
            }   
        }
    }
}

investigate_gene_set = function(mat){
    print(dim(mat)[1])
    print(cor.test(mat[,1],mat[,2]))     
}

investigate_enrichement = function(set, all_genes){
    deg = 1*(all_genes %in% set)
    names(deg) = all_genes
    pwf = nullp(deg, 'mm10', 'geneSymbol', plot.fit=F)
    res = matrix(0,nrow=2,ncol=2, dimnames=list(c("over", "under"),c("GO","KEGG")))
    # GO
    GO_wall = goseq(pwf,'mm10', 'geneSymbol')
    over_represented_GO = GO_wall[GO_wall$over_represented_pvalue < 0.05,c("category","term","ontology")]
    under_represented_GO = GO_wall[GO_wall$under_represented_pvalue < 0.05,c("category","term","ontology")]
    res[1,1] = dim(over_represented_GO)[1]
    res[2,1] = dim(under_represented_GO)[1]
    # plot ontology barplot
    GO_ontology_counts = merge(count(over_represented_GO, var="ontology"), count(under_represented_GO, var="ontology"), by="ontology")
    rownames(GO_ontology_counts) = GO_ontology_counts[,1]
    GO_ontology_counts = GO_ontology_counts[,-1]
    colnames(GO_ontology_counts) = c("over_represented_GO", "under_represented_GO")
    GO_ontology_counts = as.matrix(GO_ontology_counts)
    barplot(t(GO_ontology_counts), beside = TRUE, col=c("green4", "red4"))
    legend("topright", c("over", "under"), fill=c("green4", "red4"))
    # KEGG pathways
    KEGG_wall = goseq(pwf,'mm10', 'geneSymbol', test.cats="KEGG")
    over_represented_KEGG = KEGG_wall[KEGG_wall$over_represented_pvalue < 0.05,]
    under_represented_KEGG = KEGG_wall[KEGG_wall$under_represented_pvalue < 0.05,]
    res[1,2] = dim(over_represented_KEGG)[1]
    res[2,2] = dim(under_represented_KEGG)[1]
    print(res)
}


order_by_min_na = function(data){
    na_nb = apply(is.na(data), 2, sum)
    return(data[order(data[,which.min(na_nb)]),])
}

plot_fc_heatmap_with_modules = function(data, fc_annot, connected_gene_colors){
    module_nb = length(unique(connected_gene_colors))
    # get genes per modules and order inside module by log2FC of the comparison with the less NA
    module_gene_fc_deg = matrix(0, ncol = 0, nrow = 0)
    for(i in 1:module_nb){
        genes_in_mod = names(connected_gene_colors)[connected_gene_colors == i & names(connected_gene_colors) %in% rownames(data)]
        module_gene_fc_deg = rbind(module_gene_fc_deg, order_by_min_na(data[genes_in_mod,]))
    }
    # add genes not in modules (also ordered)
    genes_in_mod_to_keep = rownames(module_gene_fc_deg)
    genes_not_in_mod = rownames(data)[!rownames(data) %in% genes_in_mod_to_keep]
    module_gene_fc_deg = rbind(module_gene_fc_deg, order_by_min_na(data[genes_not_in_mod,]))
    # add annotation
    annot_row = data.frame(module = as.factor(c(connected_gene_colors[genes_in_mod_to_keep], rep("No module", length(genes_not_in_mod)))))
    rownames(annot_row) = c(genes_in_mod_to_keep, genes_not_in_mod)
    mod_pal = c(head(pal2, -2),'white')
    names(mod_pal) = as.factor(c(1:module_nb,"No module"))
    annot_colors = list(
        module = mod_pal
    )
    # get breaks for the colors
    data_vector = unlist(module_gene_fc_deg)
    data_vector = data_vector[!is.na(data_vector)]
    breaks = quantile_breaks(data_vector, n = 11)
    # plot heatmap
    pheatmap(module_gene_fc_deg,
            cluster_rows=F,
            cluster_cols=F,
            show_rownames=F,
            show_colnames=F,
            annotation_col=annot_col,
            annotation_row=annot_row,
            annotation_colors = annot_colors,
            breaks=breaks,
            color=rev(brewer.pal(10, "RdBu")))
}


order_rows = function(z){
    # perform hierarchical clustering inside the module to order the genes
    hc = hclust(dist(z), method = "complete")
    return(z[hc$order,])
}

get_genes_in_mod = function(mod_id, connected_gene_colors, gene_subset){
    # Get the genes in a given module and in the subset of genes
    return(names(connected_gene_colors)[connected_gene_colors == mod_id & 
                                        names(connected_gene_colors) %in% gene_subset])
}


plot_z_score_heatmap_with_modules = function(z_scores, deg, col_order, annot_col, title){
    # Plot the Z-score heatmap with module on the left
    
    # get z_score for the DE genes and correct column order
    data = z_scores[deg, col_order]
    # get module nb
    module_nb = length(unique(connected_gene_colors))
    # get genes per modules
    module_gene_scores = matrix(0, ncol = dim(norm_counts)[2], nrow = 0)
    for(mod_id in unique(connected_gene_colors)){
        genes_in_mod = get_genes_in_mod(mod_id, connected_gene_colors, rownames(data))
        # add to all genes after a hierarchical clustering inside the module to order the genes
        if(length(genes_in_mod) != 0){
            z = data[genes_in_mod,]
            if(length(genes_in_mod) < 2){
                z = as.matrix(t(z))
                rownames(z) = genes_in_mod
            }else{
                z = order_rows(z)
            }
            module_gene_scores = rbind(module_gene_scores, z)
        }
    }
    # add genes not in modules (also ordered)
    genes_in_mod_to_keep = rownames(module_gene_scores)
    genes_not_in_mod = rownames(data)[!rownames(data) %in% genes_in_mod_to_keep]
    module_gene_scores = rbind(module_gene_scores, order_rows(data[genes_not_in_mod,]))
    # add annotation
    annot_row = data.frame(module = as.factor(c(paste("ME", connected_gene_colors[genes_in_mod_to_keep], sep=""), 
                                                rep("No module", length(genes_not_in_mod)))))
    rownames(annot_row) = c(genes_in_mod_to_keep, genes_not_in_mod)
    # plot heatmap
    pheatmap(module_gene_scores,
            cluster_rows=F,
            cluster_cols=F,
            show_rownames=F,
            show_colnames=F,
            annotation_col=annot_col,
            annotation_row=annot_row,
            annotation_colors = annot_colors,
            color=rev(brewer.pal(11, "RdBu")),
            breaks = seq(-3.5, 3.5, length=11),
            main = title)
}


plot_z_score_heatmap = function(z_scores, de_genes, col_order, annot_col, title, col_for_clust){
    # get z_score for the DE genes and correct column order
    data = z_scores[de_genes,]
    # cluster rows
    hc = hclust(dist(data[,col_for_clust]), method = "complete")
    # plot
    pheatmap(data[hc$order, col_order],
             cluster_rows=F,
             cluster_cols=F,
             show_rownames=F,
             show_colnames=F,
             annotation_col=annot_col,
             annotation_row=NULL,
             annotation_colors = NULL,
             color=rev(brewer.pal(11, "RdBu")),
             breaks = seq(-3.5, 3.5, length=11),
             main = title)
}


plot_module_groups = function(groups, vertsep){
    # Calculate correlation
    moduleTraitCor = cor(MEs, trait, use = "p")
    moduleTraitPvalue = corPvalueStudent(moduleTraitCor, dim(filtered_norm_counts)[2])
    # perform hierarchical clustering of the modules
    hc = hclust(dist(moduleTraitCor), method = "complete")
    moduleTraitCor = moduleTraitCor[hc$order,]
    moduleTraitPvalue = moduleTraitPvalue[hc$order,]
    # Will display correlations and their p-values
    textMatrix = paste(signif(moduleTraitCor, 2), "\n(",signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) = dim(moduleTraitCor)
    # Row cols
    row_col = paste("ME", pal2, sep="")
    names(row_col) = paste("ME", names(pal2), sep="")
    row_col = row_col[match(rownames(moduleTraitCor), names(row_col))]
    # Row names
    row_names = sapply(names(row_col), function(x){
        y = if(x == 'ME0') 'No module' else x 
        paste(y, " \n(", mod_sizes[x], " genes)", sep = "")})
    # Plot correlation
    par(mar = c(6, 8.5, 3, 3))
    labeledHeatmap(Matrix = moduleTraitCor,
                    xLabels = colnames(trait),
                    yLabels = row_col,
                    ySymbols = row_names,
                    colorLabels = FALSE,
                    yColorLabels = TRUE,
                    colors = rev(brewer.pal(11, "RdBu")),
                    textMatrix = textMatrix,
                    setStdMargins = FALSE,
                    cex.text = 0.5,
                    cex.lab.y = .75,
                    zlim = c(-1,1),
                    verticalSeparator.x = vertsep)
}

plot_top_deg_in_modules = function(fc, comp, connected_gene_colors){
    gene_subset = rownames(fc)[!is.na(fc[,comp])]
    # 
    modules = unique(connected_gene_colors)
    mod_genes = lapply(modules, function(mod){
        # get DE genes in module
        mod_de_gene = get_genes_in_mod(mod, connected_gene_colors, gene_subset)
        # get FC for the genes
        mod_fc = fc[mod_de_gene, comp]
        names(mod_fc) = mod_de_gene
        mod_fc = sort(mod_fc, decreasing = TRUE)
        # get most up and down DEG
        neg_genes = names(tail(mod_fc[mod_fc < 0], n = 10))
        pos_genes = names(head(mod_fc[mod_fc > 0], n = 10))
        return(list('pos_genes' = pos_genes, 'neg_genes' = neg_genes, 'max' = max(c(length(pos_genes), length(neg_genes)))))
    })
    names(mod_genes) = modules
    # plot empty plot              
    div = 2
    top_neg = 10*1/div
    top_pos = 2*top_neg + .5
    top = top_pos + 1
    w_offset = .5
    w = 2*length(modules) + w_offset
    h = top
    prev_w = 7
    prev_h = 7
    options(repr.plot.width=w, repr.plot.height=h)
    par(mar = c(0,0,0,0))
    plot(c(0, w), c(0, h), type= "n", axes=FALSE, ann=FALSE)
    text(w/2, top, paste("Top DE genes in", comp), cex = 2.8, font = 2)
    text(0, top_neg + .5 + top_neg/2,"FC > 0", srt = 90, cex = 2)
    text(0, top_neg/2,"FC < 0", srt = 90, cex = 2)                       
    # parse modules
    t = sapply(1:length(mod_genes), function(i){
        mod = names(mod_genes)[i]
        mod_info = mod_genes[[mod]]
        # plot rect on the top
        rect((i-1)*2 + w_offset, top-.8, i*2 + w_offset, top - 1, col = pal2[mod], border = NA)
        text((i-1)*2 + 1 + w_offset, top-.7, paste("ME", mod, sep = ''), cex = 1.5, adj = c(0.5, 0))
        # add pos gene names to plot
        pos_genes = mod_info[['pos_genes']]
        if(length(pos_genes) > 0){
            t = sapply(1:length(pos_genes), function(x){
                text((i-1)*2 + w_offset, top_pos - x/div, pos_genes[x], cex = 1.2, adj = c(0, 0.5), offset = c(1,NA))
            })
        }
        # add neg gene names to plot
        neg_genes = mod_info[['neg_genes']]
        if(length(neg_genes) > 0){
            t = sapply(1:length(neg_genes), function(x){
                text((i-1)*2 + w_offset, top_neg - x/div, neg_genes[x], cex = 1.2, adj = c(0, 0.5), offset = c(1,NA))
            })
        }
    })
}