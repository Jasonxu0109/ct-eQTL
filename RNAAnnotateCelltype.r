RNAAnnotateCelltype <- function (RNA, genes, signatures = "human.immune.CIBERSORT", 
    min.score = 0, orig.ident = NULL, outdir = ".") 
{
    if (is.null(orig.ident)) {
        if (class(signatures) == "character") {
            signatures_name = signatures
            data(list = signatures)
            signatures = get(signatures)
        }
        else {
            signatures_name = ""
        }
        celltypes <- as.character(unique(signatures[, 1]))
        signature_list <- sapply(1:length(celltypes), function(x) {
            return(toupper(as.character(signatures[which(signatures[, 
                1] == celltypes[x]), 2])))
        })
        names(signature_list) <- celltypes
        cluster_celltype_score = sapply(as.integer(unique(RNA@meta.data$seurat_clusters)) - 
            1, function(x) {
            idx = genes$cluster == x
            avglogFC = genes$avg_logFC[idx]
            names(avglogFC) = toupper(genes$gene[idx])
            score_cluster = sapply(signature_list, function(y) {
                score = sum(avglogFC[y], na.rm = TRUE)/log2(length(y))
                return(score)
            })
        })
		cluster_celltype_score[is.na(cluster_celltype_score)] <- 0
		cluster_celltype_score[is.infinite(cluster_celltype_score)] <- 0
        colnames(cluster_celltype_score) = as.character(as.integer(unique(RNA@meta.data$seurat_clusters)) - 
            1)
        cellscore_max = apply(cluster_celltype_score, 2, max, 
            na.rm = TRUE)
        cellscore_max_celltype = apply(cluster_celltype_score, 
            2, function(x) {
                if (max(x) < min.score) {
                  return("Others")
                }
                else {
                  return(rownames(cluster_celltype_score)[which.max(x)])
                }
            })
        RNA@meta.data$labels = as.character(as.integer(RNA@meta.data$seurat_clusters) - 
            1)
        current.cluster.ids = as.character(as.integer(unique(RNA@meta.data$seurat_clusters)) - 
            1)
        new.cluster.ids = cellscore_max_celltype
        RNA@meta.data$assign.score = cellscore_max[RNA@meta.data$labels]
        RNA@meta.data$labels = plyr::mapvalues(x = RNA@meta.data$labels, 
            from = current.cluster.ids, to = new.cluster.ids)
        p = DimPlot(object = RNA, label = TRUE, pt.size = 0.2, 
            group.by = "labels", label.size = 3, repel = T)
        ggsave(file.path(outdir, paste0(RNA@project.name, "_annotated.png")), 
            p, width = 6, height = 4)
    }
    else {
        RNA$labels <- orig.ident
        p = DimPlot(object = RNA, label = TRUE, pt.size = 0.2, 
            group.by = "labels", label.size = 3, repel = T)
        ggsave(file.path(outdir, paste0(RNA@project.name, "_original.png")), 
            p, width = 5.5, height = 4)
    }
    return(RNA)
}