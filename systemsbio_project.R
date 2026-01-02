# ==============================================================================
# Biyoinformatik Ödevi - TAM ENTEGRE PROFESYONEL ANALİZ
# Özellikler: 
# 1. Performans için 20+20 örnekleme (Subsampling)
# 2. clusterProfiler ile detaylı GO Görselleştirmesi (Dotplot, Cnetplot, Treeplot)
# 3. Network, SVM ve Heatmap analizleri
# ==============================================================================

# --- 0. Kütüphaneler ---
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")

# Temel ve Görselleştirme Paketleri
packages <- c("GEOquery", "e1071", "pheatmap", "ppcor", "hgu219.db")
# GO Analizi Paketleri (ClusterProfiler)
packages_go <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE")

all_packages <- c(packages, packages_go)

for (p in all_packages) {
  if (!require(p, character.only = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
    library(p, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 1. VERİ İNDİRME VE SEÇME (20 Healthy + 20 Tumor)
# ------------------------------------------------------------------------------
cat("Veri indiriliyor (GSE44076)...\n")
gse <- getGEO("GSE44076", GSEMatrix = TRUE, AnnotGPL = FALSE)
eset <- gse[[1]]
Full_ExpressionData <- exprs(eset)

# Grup Etiketlerini Oluşturma
pdata <- pData(eset)
Full_GroupVector <- rep(NA, nrow(pdata))
Full_GroupVector[grep("normal", pdata$source_name_ch1, ignore.case = TRUE)] <- "Normal"
Full_GroupVector[grep("tumor|cancer|adenocarcinoma", pdata$source_name_ch1, ignore.case = TRUE)] <- "Tumor"

valid_idx <- !is.na(Full_GroupVector)
Full_ExpressionData <- Full_ExpressionData[, valid_idx]
Full_GroupVector <- Full_GroupVector[valid_idx]

# --- RASTGELE SEÇİM (SUBSAMPLING) ---
cat("Bilgisayarın çökmemesi için 20 Healthy + 20 Tumor seçiliyor...\n")
set.seed(102) 

normal_indices <- which(Full_GroupVector == "Normal")
tumor_indices <- which(Full_GroupVector == "Tumor")

n_select <- 20
sel_normal <- sample(normal_indices, min(length(normal_indices), n_select))
sel_tumor <- sample(tumor_indices, min(length(tumor_indices), n_select))

selected_indices <- c(sel_normal, sel_tumor)
ExpressionData <- Full_ExpressionData[, selected_indices]
GroupVector <- factor(Full_GroupVector[selected_indices], levels = c("Normal", "Tumor"))

cat("Analiz toplam", ncol(ExpressionData), "örnek ile devam ediyor.\n")

# ------------------------------------------------------------------------------
# 2. PCA ANALİZİ
# ------------------------------------------------------------------------------
cat("PCA Grafiği çiziliyor...\n")
pca_res <- prcomp(t(ExpressionData), scale = TRUE)
pca_summary <- summary(pca_res)
pc1_var <- round(pca_summary$importance[2,1] * 100, 2)
pc2_var <- round(pca_summary$importance[2,2] * 100, 2)

pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = GroupVector)

plot(pca_df$PC1, pca_df$PC2, col = ifelse(pca_df$Group=="Normal", "green", "red"),
     pch = 19, main = "PCA: Normal vs Tumor (Subsampled)", 
     xlab = paste0("PC1 (", pc1_var, "%)"), ylab = paste0("PC2 (", pc2_var, "%)"))
legend("topright", legend=levels(GroupVector), col=c("green", "red"), pch=19)

# ------------------------------------------------------------------------------
# 3. T-TEST VE DEG ANALİZİ
# ------------------------------------------------------------------------------
cat("Diferansiyel Gen İfadesi (DEG) Analizi yapılıyor...\n")
p_values <- apply(ExpressionData, 1, function(x) {
  if(var(x) == 0) return(1)
  t.test(x[GroupVector == "Normal"], x[GroupVector == "Tumor"])$p.value
})

p_adj <- p.adjust(p_values, method = "BH")
deg_indices <- which(p_adj < 0.05)
DEG_Data <- ExpressionData[deg_indices, ]

# Genleri sırala
sorted_indices <- order(p_values[deg_indices])
DEG_Data_Sorted <- DEG_Data[sorted_indices, ] 
top5_genes <- rownames(DEG_Data_Sorted)[1:5]

cat("Anlamlı Gen Sayısı:", nrow(DEG_Data_Sorted), "\n")
cat("Top 5 Gen:", paste(top5_genes, collapse=", "), "\n")

# ------------------------------------------------------------------------------
# 4. DETAYLI GO ANALİZİ (ClusterProfiler ENTEGRASYONU)
# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# 4. DETAYLI GO ANALİZİ (DÜZELTİLMİŞ VERSİYON)
# ------------------------------------------------------------------------------
cat("\n--- DETAYLI GO ANALİZİ VE GÖRSELLEŞTİRME BAŞLIYOR ---\n")

# A. ID Dönüşümü (Prob ID -> Entrez ID)
# HATA ÇÖZÜMÜ: Dönüşüm için 'hgu219.db' kullanıyoruz.
# -------------------------------------
library(hgu219.db)

genes_of_interest <- rownames(DEG_Data_Sorted)
all_genes_universe <- rownames(ExpressionData)

# 1. Anlamlı Genler için Dönüşüm
deg_entrez <- mapIds(hgu219.db, 
                     keys = genes_of_interest, 
                     column = "ENTREZID", 
                     keytype = "PROBEID",
                     multiVals = "first")

# NA (eşleşmeyen) değerleri temizle
deg_entrez <- deg_entrez[!is.na(deg_entrez)]

# 2. Tüm Genler (Universe) için Dönüşüm
universe_entrez <- mapIds(hgu219.db, 
                          keys = all_genes_universe, 
                          column = "ENTREZID", 
                          keytype = "PROBEID",
                          multiVals = "first")

# NA değerleri temizle
universe_entrez <- universe_entrez[!is.na(universe_entrez)]

cat("Gen dönüşümü tamamlandı. GO Analizi başlatılıyor...\n")

# B. Zenginleştirme Analizi (Biological Process)
# Artık elimizde ENTREZ ID olduğu için org.Hs.eg.db kullanabiliriz.
# ----------------------------------------------
ego <- enrichGO(gene          = deg_entrez,
                universe      = universe_entrez,
                OrgDb         = org.Hs.eg.db,
                ont           = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                readable      = TRUE) # Sonuçlarda Gen Sembollerini göster

# Eğer sonuç boş dönerse uyarı ver
if(is.null(ego)) {
  cat("UYARI: GO analizi anlamlı bir sonuç bulamadı. P-value sınırını gevşetmeyi deneyebilirsiniz.\n")
} else {
  cat("GO Analizi tamamlandı. Grafikler çiziliyor...\n")
  
  # C. Grafikler
  # ------------
  # 1. Dotplot
  print(dotplot(ego, showCategory=15, title="Top 15 Biological Processes (Dotplot)"))
  
  # 2. Barplot
  print(barplot(ego, showCategory=15, title="Top 15 Biological Processes (Barplot)"))
  
  # 3. Cnetplot (Gen-Kavram Ağı)
  tryCatch({
    p_cnet <- cnetplot(ego, showCategory = 5, circular = FALSE, colorEdge = TRUE)
    print(p_cnet)
  }, error = function(e) { cat("Cnetplot çizilemedi.\n") })
  
  # 4. Treeplot
  tryCatch({
    ego_sim <- pairwise_termsim(ego)
    p_tree <- treeplot(ego_sim, showCategory = 30)
    print(p_tree)
  }, error = function(e) { cat("Treeplot çizilemedi.\n") })
}

cat("ClusterProfiler analizleri tamamlandı.\n")

# ------------------------------------------------------------------------------
# 5. HEATMAP (HEATMAP)
# ------------------------------------------------------------------------------
cat("\nHeatmap çiziliyor...\n")
max_genes_heatmap <- 1000 
DEG_Data_Heatmap <- if(nrow(DEG_Data_Sorted) > max_genes_heatmap) DEG_Data_Sorted[1:max_genes_heatmap, ] else DEG_Data_Sorted

pheatmap(DEG_Data_Heatmap, 
         annotation_col = data.frame(Group = GroupVector, row.names = colnames(DEG_Data_Heatmap)),
         scale = "row", show_rownames = FALSE, main = "Hierarchical Clustering (Top DEGs)")

# ------------------------------------------------------------------------------
# 6. SVM SINIFLANDIRMA
# ------------------------------------------------------------------------------
cat("SVM Modeli çalıştırılıyor...\n")
# Son 5'er örneği test için ayır
test_indices <- c(which(GroupVector=="Normal")[20], which(GroupVector=="Tumor")[20])

svm_model <- svm(t(DEG_Data[, -test_indices]), GroupVector[-test_indices], kernel = "linear")
preds <- predict(svm_model, t(DEG_Data[, test_indices]))

res_svm <- data.frame(Sample = rownames(t(DEG_Data[, test_indices])), Gerçek = GroupVector[test_indices], Tahmin = preds)
print(res_svm)

# ------------------------------------------------------------------------------
# 7. AĞ ANALİZİ (PCOR SKORLU)
# ------------------------------------------------------------------------------
cat("\n--- AĞ ANALİZİ (Top 5 Gen) ---\n")

build_network_with_score <- function(expression_matrix, gene_names) {
  results_list <- list()
  n <- length(gene_names)
  data <- t(expression_matrix)
  
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      edge_exists <- TRUE
      min_abs_score <- 1
      final_score <- 0
      
      other_genes <- setdiff(1:n, c(i, j))
      if(length(other_genes) == 0) {
        res <- cor.test(data[,i], data[,j])
        final_score <- res$estimate
        if(abs(final_score) <= 0.10) edge_exists <- FALSE
      } else {
        for (k in other_genes) {
          tryCatch({
            res <- pcor.test(data[,i], data[,j], data[,k])
            est <- res$estimate
            if (!is.na(est)) {
              if (abs(est) <= 0.10) { edge_exists <- FALSE; break }
              if (abs(est) < min_abs_score) { min_abs_score <- abs(est); final_score <- est }
            }
          }, error = function(e) { })
        }
      }
      if (edge_exists) {
        results_list[[length(results_list)+1]] <- data.frame(Gene1 = gene_names[i], Gene2 = gene_names[j], Score = round(final_score, 4))
      }
    }
  }
  if(length(results_list) > 0) return(do.call(rbind, results_list)) else return(NULL)
}

cat("Normal Grubu Ağı:\n")
data_normal_sub <- ExpressionData[top5_genes, GroupVector == "Normal"]
net_normal <- build_network_with_score(data_healthy_sub, top5_genes)
if(!is.null(net_healthy)) print(net_healthy) else print("Normal grubunda bağlantı yok.")

cat("\nTumor Grubu Ağı:\n")
data_tumor_sub <- ExpressionData[top5_genes, GroupVector == "Tumor"]
net_tumor <- build_network_with_score(data_tumor_sub, top5_genes)
if(!is.null(net_tumor)) print(net_tumor) else print("Tumor grubunda bağlantı yok.")

# ------------------------------------------------------------------------------
# 8. GEN İSİMLERİ (MAPPING)
# ------------------------------------------------------------------------------
cat("\n--- PROB ID -> GEN SEMBOLÜ ---\n")
gene_symbols <- mapIds(hgu219.db, keys = top5_genes, column = "SYMBOL", keytype = "PROBEID")
print(data.frame(ProbeID = top5_genes, GeneSymbol = gene_symbols))

cat("\nBütün analizler başarıyla tamamlandı.\n")

