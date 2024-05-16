createSortedVector <- function(n) {
  char_numbers = as.character(1:n)
  sorted_char_numbers = sort(char_numbers, method = "radix")
  sorted_numbers = as.integer(sorted_char_numbers)
  return(sorted_numbers)
}

# Function to process each image and write to FCS
write_fcs <- function(data, image_id) {
  # Remove 'ImageID' and 'ImageName' column and ensure the rest is numeric
  subset_data <- as.matrix(data[, -c(1, (ncol(data)-1):ncol(data))])
  subset_data <- apply(subset_data, 2, as.numeric)
  # Ensure column names are set
  colnames(subset_data) <- make.names(colnames(subset_data))
  # Convert to flowFrame
  ff <- flowFrame(subset_data)
  # Write to FCS
  fcs_filename <- paste0("../FlowJo/", image_id, ".fcs")
  write.FCS(ff, file = fcs_filename)
}

## Subset cells, batch correct and  SOM cluster
harmonySOM <- function(spe, compartment_vector = NULL, celltype_column = NULL, celltype_values = NULL, clustMarkers = NULL, batchCol = NULL) {
  # Function processes single-cell data for UMAP, PCA, and SOM clustering.
  # - spe: Single-cell experiment object.
  # - celltype_column, celltype_values, compartment_vector: Optional parameters for subsetting.
  # - clustMarkers: Optional. Markers to use for clustering. Uses 'use_channel' if not provided.
  # - batchCol: Optional. Column name for batch correction groups. If NULL, skip batch correction.
  # Returns a list containing the SOM clustered data and subsetted spe data
  
  # Assign variable names
  speData = "norm_exprs"
  pcaName = "PCA_norm"
  umapName = "UMAP"  # Unified UMAP reduction name
  
  # Conditional subsetting
  subCells <- rep(TRUE, ncol(spe))
  if (!is.null(celltype_column) && !is.null(celltype_values)) {
    subCells <- subCells & (!is.na(colData(spe)[[celltype_column]])) & (colData(spe)[[celltype_column]] %in% celltype_values)
  }
  if (!is.null(compartment_vector)) {
    subCells <- subCells & (colData(spe)$Compartment %in% compartment_vector)
  }
  
  spe_subset <- spe[, subCells]
  
  # Update 'use_channel' based on clustMarkers
  if (!is.null(clustMarkers)) {
    if (!all(clustMarkers %in% rownames(spe_subset))) {
      stop("Some of the specified markers for clustering do not exist in rowData.")
    }
    rowData(spe_subset)$use_channel <- rownames(spe_subset) %in% clustMarkers
  } else {
    if (!"use_channel" %in% names(rowData(spe_subset))) {
      stop("The 'use_channel' does not exist in rowData.")
    }
  }
  
  # Run PCA on the subset
  spe_subset <- runPCA(spe_subset, subset_row = rowData(spe_subset)$use_channel, exprs_values = speData, ncomponents = 30, BSPARAM = ExactParam())
  reducedDims(spe_subset)[[pcaName]] <- reducedDims(spe_subset)[["PCA"]]
  
  # Harmony batch correction (if specified)
  if (!is.null(batchCol)) {
    set.seed(230616)
    out <- RunHarmony(spe_subset, group.by.vars = batchCol)
    stopifnot(all.equal(colnames(spe_subset), colnames(out)))
    reducedDim(spe_subset, umapName) <- reducedDim(out, "HARMONY")
  } else {
    # Run UMAP without batch correction
    set.seed(220228)
    spe_subset <- runUMAP(spe_subset, subset_row = rowData(spe_subset)$use_channel, exprs_values = speData, name = umapName)
  }
  
  # Perform SOM clustering
  set.seed(220410)
  som.out <- clusterRows(reducedDim(spe_subset, umapName), SomParam(100), full = TRUE)
  
  return(list(som.out = som.out, spe_subset = spe_subset))
}


#Example Useage
# clustOut = harmonySOM(spe, "LA", "celltype_flowjo", "Tcell")


## choose cluster number, generate graphs and output spe_subset with cluster numbers annotated
clustGraphs <- function(clustOutF, ccpF, ccp_clustNum, graph_type, annotate_labels, save_path,  assayName = "norm_exprs") {
  # Function clustGraphs: Generates and saves UMAP, Heatmap, and Z-score graphs based on clustering results.
  # - clustOutF: Output of harmonySOM function which contains in entry 1 the SOM clustered data and entry 2 the subsetted SingleCellExperiment object with clustering information.
  # - ccpF: is the ccp variable which holds the output of ConsensusClusterPlus Clustering.
  # - ccp_clustNum: Numeric value of clusters from ConsensusClusterPlus.
  # - graph_type: Type of graph to generate ('All', 'UMAP', 'Heatmap', 'Zscore').
  # - annotate_labels: Boolean to indicate whether to annotate labels on graphs.
  # - save_path: Folder path for saving the generated graphs.
  # Returns speF which has the cluster number annotations added as a metadata column
  
  set.seed(220228)
  
  somOutput = clustOutF[[1]]
  speF = clustOutF[[2]]
  # Check if save_path has a trailing slash and add one if it doesn't
  if (substr(save_path, nchar(save_path), nchar(save_path)) != "/") {
    save_path <- paste0(save_path, "/")
  }
  
  # Convert ccp_clustNum to string for file naming
  cluster_str <- as.character(ccp_clustNum)
  
  # Link ConsensusClusterPlus clusters with SOM codes and save in object
  som.cluster <- ccpF[[ccp_clustNum]][["consensusClass"]][somOutput$clusters]
  speF$som_clusters_corrected <- as.factor(som.cluster)
  
  # Determine graph types and generate accordingly
  if (graph_type %in% c("All", "UMAP")) {
    
    p1 <- dittoDimPlot(speF, var = "som_clusters_corrected", 
                       reduction.use = "UMAP", size = 1.2,
                       do.label = annotate_labels, labels.size = 2.5,
                       legend.size = 4) +
      ggtitle("SOM clusters on UMAP, integrated cells")
    print(p1)
    ggsave(paste0(save_path, cluster_str, "_UMAP_labelled.png"), p1, width = 6, height = 4.5)
  }
  
  if (graph_type %in% c("All", "Heatmap")) {
    # Choose the smaller number between the total number of cells and 2000
    cell_count <- min(ncol(speF), 2000)
    cur_cells <- sample(seq_len(ncol(speF)), cell_count)
    p1 <- dittoHeatmap(speF[,cur_cells], 
                       genes = rownames(speF)[rowData(speF)$use_channel],
                       assay = assayName, scale = "none",
                       heatmap.colors = viridis(100), 
                       annot.by = c("som_clusters_corrected", "DonorID"),
                       annot.colors = c(dittoColors(1)[1:length(unique(speF$som_clusters_corrected))],
                                        metadata(speF)$color_vectors$DonorID))
    print(p1)
    ggsave(paste0(save_path, cluster_str, "_Heatmap.png"), p1, width = 8.5, height = 6, dpi = 300)
  }
  
  if (graph_type %in% c("All", "Zscore")) {
    celltype_mean <- aggregateAcrossCells(as(speF, "SingleCellExperiment"),  
                                          ids = speF$som_clusters_corrected, 
                                          statistics = "mean",
                                          use.assay.type = "exprs")
    p1 <- dittoHeatmap(celltype_mean,
                       genes = rownames(speF)[rowData(speF)$use_channel], 
                       assay = "exprs", 
                       cluster_rows = TRUE, 
                       annot.by = c("som_clusters_corrected"))
    print(p1)
    ggsave(paste0(save_path, cluster_str, "_Zscore.png"), p1, width = 8.5, height = 6, dpi = 300)
  }
  
  return(speF)
}

# Example Useage
# spe_subset = clustGraphs(clustOut, ccp, 12, "All", TRUE, "LA_Tcell")


## Export new clusters for validation in Mantis Viewer
mantis_pop <- function(image_name, suffix, cellSubs) {
  df = celltype_by_image[[image_name]]
  rownames(df) = df$CellID
  df = as.data.frame(df) %>%
    select(cellSubs)
  write.table(df, paste0("mantis_CSVs/", image_name, suffix, ".csv"), sep=",",  col.names=FALSE)
}

## Visualise marker histogram to find cut of for filtering cells
visMarker <- function(speF, target, marker, value, celltype_col = "celltype", assayF = "norm_exprs", group_by = "DonorID", overall_title = NULL, target_title = NULL) {
  # Remove NAs
  plotData <- speF[, !is.na(colData(speF)[[celltype_col]]) & colData(speF)[[celltype_col]] == target]
  plotData$All <- rep("All", length(colnames(plotData)))
  
  # Create titles if not provided
  if (is.null(overall_title)) {
    overall_title <- paste("Overall Expression of", marker, "in", target, "cells")
  }
  if (is.null(target_title)) {
    target_title <- paste("Expression of", marker, "in", target, "cells by", group_by)
  }
  
  # Create the first plot by specified group
  plot1 <- dittoRidgePlot(plotData, var = marker, group.by = group_by, assay = assayF, max = 1) +
    ggtitle(target_title) +
    geom_vline(xintercept = value, col = "black", lwd = 1.5, linetype = "dashed")
  
  # Create a second plot for all cells together
  plot2 <- dittoRidgePlot(plotData, var = marker, group.by = "All", assay = assayF, max = 1) +
    ggtitle(overall_title) +
    geom_vline(xintercept = value, col = "red", lwd = 1.5, linetype = "dashed")
  
  # Print both plots on the same screen
  grid.arrange(plot1, plot2, ncol = 2)
}


# Example usage:
# visMarker(spe_subset, "celltype", "Cycling DC", "CD20", 0.5, group_by = "DonorID")


## Track population reassignments in real-time in a table written to working directory

# Check if the tracking_table already exists in the global environment
if (file.exists("clusterModTracker.csv")) {
  # Load the existing tracking table from the file
  tracking_table <<- read.csv("clusterModTracker.csv", stringsAsFactors = FALSE)
} else {
  # Create a new tracking table if the file doesn't exist
  tracking_table <<- data.frame(
    from = character(),
    to = character(),
    condition = character(),
    stringsAsFactors = FALSE
  )
}

update_tracking <- function(from_pop, to_pop, marker, value, direction) {
  
  # Create a new entry
  new_entry <- data.frame(from_pop = from_pop,
                          to_pop = to_pop,
                          condition = paste(marker, direction, as.character(value)),
                          stringsAsFactors = FALSE)
  
  # Check for duplicate across all variables
  duplicate_index <- with(tracking_table, 
                          which(from_pop == new_entry$from_pop &
                                  to_pop == new_entry$to_pop &
                                  condition == new_entry$condition)
  )
  
  if (length(duplicate_index) > 0) {
    # Replace the duplicate entry
    tracking_table[duplicate_index, ] <- new_entry
  } else {
    # Append the new entry
    tracking_table <<- rbind(tracking_table, new_entry)
  }
  
  write.csv(tracking_table, "clusterModTracker.csv", row.names = FALSE)
}

# Example Usage
# update_tracking("Cycling DC", "B cell", "CD20", 0.5, ">")
# update_tracking("Memory T", "Treg", "FoxP3", 0.3, "<")



## Move cells from one population to another one based on specifications

# Assigns cells from "target" population defined in column 'celltype_col' 
# to "newTarget" in a new column called *reassigned_celltype* based on 'value' of
# 'marker' in 'assayF' being above/below 'direction'.
moveCells <- function(speF, target, newTarget, marker, value, direction, celltype_col = "celltype", assayF = "norm_exprs") {
  # Ensure the 'reassigned_celltype' column is present
  colData(speF)$reassigned_celltype <- colData(speF)[[celltype_col]]
  
  # Calculate counts before and after reassignment
  counts_before <- data.frame(colData(speF)) %>%
    dplyr::filter(reassigned_celltype == target) %>%
    group_by(DonorID) %>%
    summarize(counts_before = n())
  
  # Determine cells to reassign based on the marker expression and direction
  if (direction == ">") {
    cellRemove <- assay(speF, assayF)[marker, ] > value
  } else if (direction == "<") {
    cellRemove <- assay(speF, assayF)[marker, ] < value
  } else {
    stop("Invalid direction. Use '>' or '<'.")
  }
  
  # Reassign cell types
  cellTarget <- colData(speF)[[celltype_col]] == target
  reassigned <- cellTarget & cellRemove
  colData(speF)$reassigned_celltype[reassigned] <- newTarget
  
  counts_after <- data.frame(colData(speF)) %>%
    dplyr::filter(reassigned_celltype == target) %>%
    group_by(DonorID) %>%
    summarize(counts_after = n())
  
  # Combine the counts
  counts_combined <- merge(counts_before, counts_after, by = "DonorID", all = TRUE)
  counts_combined[is.na(counts_combined)] <- 0
  
  # Sum of counts before and after
  total_counts <- counts_combined %>%
    summarise(
      total_counts_before = sum(counts_before, na.rm = TRUE),
      total_counts_after = sum(counts_after, na.rm = TRUE)
    )
  
  # Print both outputs to the console
  print(counts_combined)
  print(total_counts)
  
  # Update tracking table
  update_tracking(target, newTarget, marker, value, direction)
  
  return(speF)
}

# Example usage:
# moveCells(spe_subset, "Cycling DC", "B cell", "CD20", 0.5, ">")