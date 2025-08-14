lis_genes <- as.list(all.genes)
matches <- grep("Ha412HOChr10g0435441", lis_genes, value = TRUE)
print(matches)

seurat_XRQ.markers_split <- seurat_XRQ.markers %>%
  separate(gene, into = c("gene_prefix", "gene_suffix"), sep = 17)

seurat_XRQ.markers_split <- seurat_XRQ.markers_split %>%
  separate(gene_prefix, into = c("gene_constant", "chr_num"), sep = 15)

seurat_XRQ.markers_split_10 <- seurat_XRQ.markers_split[seurat_XRQ.markers_split$chr_num == "10",]

unique(seurat_XRQ.markers_split_10$chr_num)

table(seurat_XRQ.markers_split_10$cluster)

clusters <- as.list(as.character(unique(seurat_XRQ.markers_split$cluster)))

for (i in 1:length(clusters)) {
  tmp = seurat_XRQ.markers_split[seurat_XRQ.markers_split$cluster == as.character(clusters[i]),]
  print(paste0("table chr occurances for cluster ", as.character(clusters[i])))
  print(table(tmp$chr_num))
}


