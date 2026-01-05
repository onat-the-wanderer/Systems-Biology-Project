# ==============================================================================
# BIOINFORMATICS PROJECT - INTEGRATED FINAL VERSION
# ==============================================================================
# Overview:
# 1. Data Downloading & Preprocessing (GSE44076)
# 2. Statistical Analysis (PCA, DEG, Heatmap)
# 3. Functional Enrichment (GO Analysis)
# 4. Machine Learning (SVM Classification)
# 5. Network Analysis using KeyPathwayMineR (KPM)
# 6. Functional Analysis of the Identified Subnetwork (KEGG)
# ==============================================================================

# --- 0. SETTINGS & MEMORY MANAGEMENT ---
# CRITICAL: Increase Java memory to 8GB to prevent KPM crashes
options(java.parameters = "-Xmx8g")
options(timeout = 10000)

# Check and Install Libraries
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("rJava")) install.packages("rJava")
library(rJava)
.jinit() # Initialize Java

# Define required packages
packages <- c("GEOquery", "e1071", "pheatmap", "ppcor", "hgu219.db", "matrixStats", "igraph", "KeyPathwayMineR", "easycsv")
packages_go <- c("clusterProfiler", "org.Hs.eg.db", "enrichplot", "DOSE", "GOSemSim")
all_packages <- c(packages, packages_go)

# Install/Load packages
for (p in all_packages) {
  if (!require(p, character.only = TRUE)) {
    BiocManager::install(p, update = FALSE, ask = FALSE)
    library(p, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 1. DATA ACQUISITION & PREPROCESSING
# ------------------------------------------------------------------------------
cat("\n--- STEP 1: Data Acquisition (GSE44076) ---\n")

gse <- NULL
tryCatch({
  gse <- getGEO("GSE44076", GSEMatrix = TRUE, AnnotGPL = FALSE)
}, error = function(e) { stop("Error downloading data. Check internet connection.") })

eset <- gse[[1]]
Full_ExpressionData <- exprs(eset)
pdata <- pData(eset)

# Create Group Labels (Normal vs Tumor)
Full_GroupVector <- rep(NA, nrow(pdata))
Full_GroupVector[grep("normal", pdata$source_name_ch1, ignore.case = TRUE)] <- "Normal"
Full_GroupVector[grep("tumor|cancer|adenocarcinoma", pdata$source_name_ch1, ignore.case = TRUE)] <- "Tumor"

# Filter valid samples
valid_idx <- !is.na(Full_GroupVector)
Full_ExpressionData <- Full_ExpressionData[, valid_idx]
Full_GroupVector <- Full_GroupVector[valid_idx]

# Subsampling (20 Normal + 20 Tumor) for performance
set.seed(102) 
sel_normal <- sample(which(Full_GroupVector == "Normal"), 20)
sel_tumor <- sample(which(Full_GroupVector == "Tumor"), 20)
selected_indices <- c(sel_normal, sel_tumor)

ExpressionData <- Full_ExpressionData[, selected_indices]
GroupVector <- factor(Full_GroupVector[selected_indices], levels = c("Normal", "Tumor"))

cat("Data Ready: ", ncol(ExpressionData), " samples selected.\n")

# ------------------------------------------------------------------------------
# 2. PCA ANALYSIS
# ------------------------------------------------------------------------------
cat("\n--- STEP 2: PCA Analysis ---\n")
pca_res <- prcomp(t(ExpressionData), scale = TRUE)
pca_df <- data.frame(PC1 = pca_res$x[,1], PC2 = pca_res$x[,2], Group = GroupVector)

# Plot PCA
plot(pca_df$PC1, pca_df$PC2, col = ifelse(pca_df$Group=="Normal", "green", "red"),
     pch = 19, main = "PCA: Normal vs Tumor", xlab="PC1", ylab="PC2")
legend("topright", legend=levels(GroupVector), col=c("green", "red"), pch=19)

# ------------------------------------------------------------------------------
# 3. DIFFERENTIAL EXPRESSION & HEATMAP
# ------------------------------------------------------------------------------
cat("\n--- STEP 3: Differential Expression Analysis ---\n")

# T-Test for each gene
p_values <- apply(ExpressionData, 1, function(x) {
  if(var(x) == 0) return(1)
  t.test(x[GroupVector == "Normal"], x[GroupVector == "Tumor"])$p.value
})

# Adjust P-values (Benjamini-Hochberg)
deg_indices <- which(p.adjust(p_values, "BH") < 0.05)
DEG_Data <- ExpressionData[deg_indices, ]
DEG_Data_Sorted <- DEG_Data[order(p_values[deg_indices]), ]
top5_genes <- rownames(DEG_Data_Sorted)[1:5]

cat("Significant DEGs found: ", nrow(DEG_Data_Sorted), "\n")

# Heatmap Visualization
max_genes_heatmap <- 1000 
DEG_Heatmap <- if(nrow(DEG_Data_Sorted) > max_genes_heatmap) DEG_Data_Sorted[1:max_genes_heatmap, ] else DEG_Data_Sorted

pheatmap(DEG_Heatmap, annotation_col = data.frame(Group = GroupVector, row.names = colnames(DEG_Heatmap)),
         scale = "row", show_rownames = FALSE, main = "Heatmap (Top DEGs)")

# ==============================================================================
# 4. GENE ONTOLOGY (GO) & KEGG ENRICHMENT (GLOBAL)
# ==============================================================================
cat("\n--- STEP 4: Global GO & KEGG Enrichment Analysis ---\n")

# 1. ID Mapping (Probe -> Entrez ID)
# --------------------------------------------------
genes_of_interest <- rownames(DEG_Data_Sorted)
universe_genes <- rownames(ExpressionData)

# DEG Entrez IDs
deg_entrez <- mapIds(hgu219.db, keys = genes_of_interest, column = "ENTREZID", keytype = "PROBEID", multiVals="first")
deg_entrez <- deg_entrez[!is.na(deg_entrez)]

# Universe Entrez IDs
universe_entrez <- mapIds(hgu219.db, keys = universe_genes, column = "ENTREZID", keytype = "PROBEID", multiVals="first")
universe_entrez <- universe_entrez[!is.na(universe_entrez)]

# 2. GO Enrichment (Biological Process)
# --------------------------------------------------
cat("Running Global GO Analysis...\n")
ego <- enrichGO(gene = deg_entrez, 
                universe = universe_entrez, 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                readable = TRUE)

if(!is.null(ego)) {
  # GO Plots
  print(dotplot(ego, showCategory=15, title="Global GO Enrichment (Biological Process)"))
  cat("GO Dotplot created.\n")
}

# 3. KEGG Enrichment (Global) - YENİ EKLENEN KISIM
# --------------------------------------------------
cat("Running Global KEGG Analysis...\n")
kk_global <- enrichKEGG(gene = deg_entrez,
                        organism = 'hsa',
                        pvalueCutoff = 0.05)

if(!is.null(kk_global)) {
  # Entrez ID'leri sembole çevir (Okunabilirlik için)
  kk_global_readable <- setReadable(kk_global, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  
  # KEGG Plots
  print(dotplot(kk_global_readable, showCategory=15, title="Global KEGG Pathways (All DEGs)"))
  
  print(cnetplot(kk_global_readable, showCategory = 5, circular = FALSE, colorEdge = TRUE, cex_label_category=0.8))
  
  cat("Global KEGG Dotplot & Cnetplot created.\n")
} else {
  cat("No significant global KEGG pathways found.\n")
}
# ------------------------------------------------------------------------------
# 5. SVM CLASSIFICATION
# ------------------------------------------------------------------------------
cat("\n--- STEP 5: SVM Classification ---\n")
test_idx <- c(which(GroupVector=="Normal")[20], which(GroupVector=="Tumor")[20])
svm_m <- svm(t(DEG_Data[,-test_idx]), GroupVector[-test_idx], kernel="linear")
preds <- predict(svm_m, t(DEG_Data[, test_idx]))

cat("SVM Predictions for Test Set:\n")
print(data.frame(Sample=names(preds), Prediction=preds, Actual=GroupVector[test_idx]))

# ------------------------------------------------------------------------------
# 6. BASIC NETWORK (PCOR)
# ------------------------------------------------------------------------------
cat("\n--- STEP 6: Basic Network Construction (PCOR) ---\n")
# Function for correlation network
build_network <- function(expr_mat, genes) {
  res_list <- list(); data <- t(expr_mat)
  n <- length(genes)
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      cor_res <- cor.test(data[,i], data[,j])
      if(abs(cor_res$estimate) > 0.5) { # Threshold
        res_list[[length(res_list)+1]] <- data.frame(Node1=genes[i], Node2=genes[j], Weight=cor_res$estimate)
      }
    }
  }
  if(length(res_list)>0) return(do.call(rbind, res_list)) else return(NULL)
}

# Build network for Top 5 genes
tumor_net <- build_network(ExpressionData[top5_genes, GroupVector=="Tumor"], top5_genes)
if(!is.null(tumor_net)) print(head(tumor_net))

# ------------------------------------------------------------------------------
# 7. KEYPATHWAYMINER (KPM) ANALYSIS & VISUALIZATION
# ------------------------------------------------------------------------------
cat("\n--- STEP 7: KeyPathwayMineR (KPM) Analysis ---\n")
library(igraph)

# A. Data Preparation for KPM
# ----------------------------------------------------
cat("Preparing data for KPM...\n")
# Map Probes to Symbols
all_symbols <- mapIds(hgu219.db, keys = rownames(ExpressionData), column = "SYMBOL", keytype = "PROBEID", multiVals="first")
valid_idx <- !is.na(all_symbols)
Exp_Sym <- ExpressionData[valid_idx, ]
rownames(Exp_Sym) <- all_symbols[valid_idx]

# Aggregate duplicate symbols
Exp_Sym <- as.matrix(aggregate(Exp_Sym, by=list(rownames(Exp_Sym)), FUN=mean)[,-1])
rownames(Exp_Sym) <- unique(sort(all_symbols[valid_idx]))

# Calculate Z-Scores and Indicator Matrix
ctrl_mean <- rowMeans(Exp_Sym[, GroupVector == "Normal"])
ctrl_sd <- matrixStats::rowSds(Exp_Sym[, GroupVector == "Normal"])
z_score <- (Exp_Sym[, GroupVector == "Tumor"] - ctrl_mean) / ctrl_sd
indicator_df <- data.frame(Gene = rownames(Exp_Sym), ifelse(abs(z_score) > 2.0, 1, 0))

# B. BioGRID File Selection & Processing
# ----------------------------------------------------
message("\n-----------------------------------------------------------")
message("ACTION REQUIRED: Please select the BioGRID file (Zip or Txt)")
message("-----------------------------------------------------------")
Sys.sleep(1)
bg_file <- file.choose()

cat("Reading BioGRID file: ", basename(bg_file), "\n")

if (grepl("\\.zip$", bg_file, ignore.case = TRUE)) {
  temp_dir <- tempdir()
  unzip(bg_file, exdir = temp_dir)
  found <- list.files(temp_dir, pattern = "\\.txt$", full.names = TRUE, recursive = TRUE)
  if(length(found)>0) bg_file <- found[1]
}

# Read and Clean BioGRID Data
lines <- readLines(bg_file, warn=FALSE)
data_lines <- grep("\t", lines, value=TRUE)
data_lines <- data_lines[!grepl("^#", data_lines)] 
bg_data <- read.delim(textConnection(data_lines), header=FALSE, stringsAsFactors=FALSE, quote="")

# Identify columns
if(ncol(bg_data) >= 9) {
  if(any(grepl("Interactor", bg_data[1,], ignore.case=TRUE))) {
    colnames(bg_data) <- bg_data[1,]
    bg_data <- bg_data[-1,]
    ints <- bg_data[, c(grep("Symbol.*A", colnames(bg_data))[1], grep("Symbol.*B", colnames(bg_data))[1])]
  } else {
    ints <- bg_data[,c(8,9)]
  }
} else {
  ints <- bg_data[,c(1,2)]
}

# Create Graph
ints <- na.omit(ints)
ints <- ints[ints[,1]!="" & ints[,1]!="-" & ints[,2]!="" & ints[,2]!="-", ]
bg_graph <- graph_from_data_frame(ints, directed=FALSE)

# Prune Graph (Speed Optimization)
common_genes <- intersect(V(bg_graph)$name, rownames(indicator_df))
bg_small <- induced_subgraph(bg_graph, vids = common_genes)

# C. Run KPM Algorithm
# ----------------------------------------------------
if(vcount(bg_small) > 0) {
  cat("Running KPM Algorithm (INES/Greedy)...\n")
  kpm_options(execution="Local", strategy="INES", algorithm="Greedy", L_min=5, L_max=5, K_min=5, K_max=5, link_type="OR")
  
  tryCatch({
    kpm_result <- kpm(indicator_mat = indicator_df, graph = bg_small)
    cat("\n>>> KPM ANALYSIS COMPLETED! Results saved. <<<\n")
    
    # D. Visualization (Manual Reconstruction)
    # ----------------------------------------------------
    message("\n-----------------------------------------------------------")
    message("ACTION REQUIRED: Please select the 'node_pathway_hits-KPM.txt' result file")
    message("(Located in Documents/computational_genomics/results/...)")
    message("-----------------------------------------------------------")
    Sys.sleep(2)
    
    res_file <- file.choose()
    cat("Selected Result File: ", basename(res_file), "\n")
    
    # Parse Result File
    raw_data <- read.table(res_file, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors=FALSE)
    candidate_genes <- raw_data[,1]
    
    # Clean Gene List
    clean_genes <- candidate_genes[!grepl("node|gene|id|hit|score|p-value|l1|k-1", candidate_genes, ignore.case=TRUE)]
    clean_genes <- clean_genes[clean_genes != ""]
    clean_genes <- unique(clean_genes)
    
    cat("Identified Genes in Subnetwork: ", length(clean_genes), "\n")
    
    # Reconstruct Subnetwork from BioGRID
    valid_genes <- intersect(clean_genes, V(bg_graph)$name)
    
    if(length(valid_genes) > 0) {
      final_net <- induced_subgraph(bg_graph, vids = valid_genes)
      
      # FIX: Use igraph::simplify explicitly
      final_net_clean <- igraph::simplify(final_net, remove.multiple = TRUE, remove.loops = TRUE)
      
      # Remove isolated nodes
      deg <- degree(final_net_clean, mode="all")
      final_net_clean <- delete.vertices(final_net_clean, V(final_net_clean)[deg == 0]) 
      
      cat("\n>>> PLOTTING SUBNETWORK... <<<\n")
      
      if(vcount(final_net_clean) > 0) {
        # Visual Settings
        V(final_net_clean)$color <- "tomato"
        V(final_net_clean)$frame.color <- "darkred"
        V(final_net_clean)$size <- 20
        V(final_net_clean)$label.color <- "black"
        V(final_net_clean)$label.font <- 2
        V(final_net_clean)$label.cex <- 0.8
        E(final_net_clean)$color <- "gray40"
        
        # Plot
        plot(final_net_clean, layout=layout_with_fr(final_net_clean), 
             main = paste("Identified Subnetwork (Cleaned)\nGene Count:", vcount(final_net_clean)))
        
        cat("Subnetwork plot created successfully.\n")
        
      } else {
        cat("Error: No connections found between identified genes.\n")
      }
    } else {
      cat("Error: Identified genes not found in BioGRID network.\n")
    }
    
  }, error = function(e) { cat("KPM Error: ", e$message) })
  
} else {
  cat("Error: No common genes between data and network.\n")
}

# ------------------------------------------------------------------------------
# 8. SUBNETWORK FUNCTIONAL ANALYSIS (KEGG & GO)
# ------------------------------------------------------------------------------
cat("\n--- STEP 8: Functional Analysis of the Subnetwork ---\n")

if(exists("final_net_clean") && vcount(final_net_clean) > 0) {
  
  sub_genes <- V(final_net_clean)$name
  
  # Map to Entrez IDs
  sub_entrez <- mapIds(hgu219.db, keys = sub_genes, column = "ENTREZID", keytype = "SYMBOL", multiVals="first")
  sub_entrez <- sub_entrez[!is.na(sub_entrez)]
  
  # KEGG Enrichment
  cat("Running KEGG Enrichment...\n")
  kk_sub <- enrichKEGG(gene = sub_entrez, organism = 'hsa', pvalueCutoff = 0.05)
  
  if(!is.null(kk_sub) && nrow(kk_sub) > 0) {
    
    # Convert IDs to Symbols for Readability
    kk_readable <- setReadable(kk_sub, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
    
    # 1. Barplot
    print(barplot(kk_readable, showCategory=10, title="KEGG Pathways (Subnetwork)"))
    
    # 2. Cnetplot (Network View)
    print(cnetplot(kk_readable, showCategory = 5, circular = FALSE, colorEdge = TRUE, 
                   node_label="all", cex_label_category=0.8, cex_label_gene=0.8))
    
    cat("KEGG plots created.\n")
    print(head(kk_sub@result[, c("ID", "Description", "p.adjust")]))
    
  } else {
    cat("No significant KEGG pathways found.\n")
  }
  
} else {
  cat("Skipping Step 8 (No subnetwork defined).\n")
}

cat("\n--- PROJECT COMPLETED SUCCESSFULLY ---\n")