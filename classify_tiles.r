

options(scipen = 999)
suppressPackageStartupMessages(library(strawr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbscan))



args    <- commandArgs( trailingOnly = TRUE )
sample  <- args[1]
genome  <- args[2]
bedpe   <- args[3]
out_dir <- args[4]


command <- paste(
  "~/aqua_tools/list_samples.sh | grep ",
  sample
)
result <- system(command, intern = T)
result <- strsplit(result, "\\|")

versions   <- sapply(result, function(x) as.integer(trimws(x[4])))
is_default <- sapply(result, function(x) trimws(x[5]) == "Y")
version    <- versions[is_default]

path_hic <- paste(
  "/home/ubuntu/lab-data", "/",
  genome, "/",
  sample, "/",
  sample,"_version_",version, "/",
  sample,"_version_",version,".allValidPairs.hic",
  sep = ""
)

whole_bedpe <- read.table(
  bedpe,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2"
  )
)


flag_cis   <- TRUE
flag_trans <- TRUE



calculate_clumpiness <- function(tile_matrix) {
  # Convert to numeric matrix if not already
  mat <- as.matrix(tile_matrix)
  
  # Calculate local density using a sliding window
  nrow <- nrow(mat)
  ncol <- ncol(mat)
  
  window_size <- 5 
  
  if(sort(dim(tile_matrix))[1] < window_size){
    window_size <- sort(dim(tile_matrix))[1]
  }
  
  local_densities <- matrix(0, nrow-window_size+1, ncol-window_size+1)
  max_counts <- matrix(0, nrow-window_size+1, ncol-window_size+1)
  
  for(i in 1:(nrow-window_size+1)) {
    for(j in 1:(ncol-window_size+1)) {
      window <- mat[i:(i+window_size-1), j:(j+window_size-1)]
      local_densities[i,j] <- sum(window > 0) / (window_size^2)
    }
  }
  
  for(i in 1:(nrow-window_size+1)) {
    for(j in 1:(ncol-window_size+1)) {
      window <- mat[i:(i+window_size-1), j:(j+window_size-1)]
      local_densities[i,j] <- sum(window > 0) / (window_size^2)
    }
  }
  
  max_density_index        <- which(local_densities == max(local_densities), arr.ind = T)[1,]
  max_density_window       <- mat[max_density_index[1]:(max_density_index[1]+window_size-1), max_density_index[2]:(max_density_index[2]+window_size-1)]
  max_density_max_count    <- mean(tail(sort(as.vector(max_density_window)),3))
  
  min_density_index        <- which(local_densities == min(local_densities), arr.ind = T)[1,]
  min_density_window       <- mat[min_density_index[1]:(min_density_index[1]+window_size-1), min_density_index[2]:(min_density_index[2]+window_size-1)]
  min_density_min_count    <- mean(head(sort(as.vector(max_density_window)),3))
  
  separation_score         <- max_density_max_count - min_density_min_count
  
  # Return metrics about clustering
  return(list(
    max_density = max(local_densities),
    density_variance = var(as.vector(local_densities)),
    separation_score = separation_score
  ))
}

sort_bedpe <- function(bedpe) {
  # Helper function to convert chromosome to numeric rank
  chr_to_rank <- function(chr) {
    chr <- gsub("chr", "", chr)
    if(grepl("^\\d+$", chr)) {  # If numeric
      return(as.numeric(chr))
    } else {  # Handle X, Y, M
      switch(chr,
             "X" = 23,
             "Y" = 24,
             "M" = 25,
             "MT" = 25,
             999)  # For any other cases, rank them last
    }
  }
  
  sorted_bedpe <- data.frame(t(apply(bedpe, 1, function(row) {
    # Case 1: Same chromosome - compare starts
    if(row["chr1"] == row["chr2"]) {
      if(as.numeric(row["start1"]) > as.numeric(row["start2"])) {
        # Swap everything
        temp_chr <- row["chr1"]
        row["chr1"] <- row["chr2"]
        row["chr2"] <- temp_chr
        
        temp_start <- row["start1"]
        row["start1"] <- row["start2"]
        row["start2"] <- temp_start
        
        temp_end <- row["end1"]
        row["end1"] <- row["end2"]
        row["end2"] <- temp_end
      }
    }
    # Case 2: Different chromosomes - compare using ranking system
    else {
      chr1_rank <- chr_to_rank(row["chr1"])
      chr2_rank <- chr_to_rank(row["chr2"])
      
      if(chr1_rank > chr2_rank) {
        # Swap everything
        temp_chr <- row["chr1"]
        row["chr1"] <- row["chr2"]
        row["chr2"] <- temp_chr
        
        temp_start <- row["start1"]
        row["start1"] <- row["start2"]
        row["start2"] <- temp_start
        
        temp_end <- row["end1"]
        row["end1"] <- row["end2"]
        row["end2"] <- temp_end
      }
    }
    return(row)
  })))
  
  return(sorted_bedpe)
}



whole_bedpe <- unique(sort_bedpe(whole_bedpe))


############################################################################
############################################################################
###                                                                      ###
###                          CIS TRANSLOCATIONS                          ###
###                                                                      ###
############################################################################
############################################################################

if(flag_cis){
  
  bedpe <- whole_bedpe[whole_bedpe$chr1 == whole_bedpe$chr2,]
  
  ##################################################################
  ##                        Gradient tiles                        ##
  ##################################################################
  
  
  feature_matrix <- matrix(
    data = 0,
    nrow = nrow(bedpe),
    ncol = 3
  )
  colnames(feature_matrix) <- c("max_density","density_variance", "separation_score")
  
  
  for( i in 1:nrow(bedpe)){
    
    chr1    <- bedpe[i,1]
    start1  <- bedpe[i,2]
    end1    <- bedpe[i,3]
    chr2    <- bedpe[i,4]
    start2  <- bedpe[i,5]
    end2    <- bedpe[i,6]
    
    tile_contacts <- straw(
      "NONE",
      path_hic,
      paste( chr1, start1, end1, sep = ":"  ),
      paste( chr2, start2, end2, sep = ":"  ),
      "BP",
      5000)
    
    tile_matrix <- matrix(
      data = 0,
      nrow = length(seq(start1,end1,5000)),
      ncol = length(seq(start2,end2,5000))
    )
    
    rownames(tile_matrix) <- seq(start1,end1,5000)
    colnames(tile_matrix) <- seq(start2,end2,5000)
    
    for( j in 1:nrow(tile_contacts)){
      x <- as.character(tile_contacts[j,"x"])
      y <- as.character(tile_contacts[j,"y"])
      
      if(x %in% rownames(tile_matrix) && y %in% colnames(tile_matrix)){
        tile_matrix[
          as.character(tile_contacts[j,"x"]),
          as.character(tile_contacts[j,"y"])] <- tile_contacts[j,"counts"]
      } else if(x %in% colnames(tile_matrix) && y %in% rownames(tile_matrix)){
        tile_matrix[
          as.character(tile_contacts[j,"y"]),
          as.character(tile_contacts[j,"x"])] <- tile_contacts[j,"counts"]
      }
    }
    
    feature_matrix[i,"max_density"]      <- calculate_clumpiness(tile_matrix)[[1]]
    feature_matrix[i,"density_variance"] <- calculate_clumpiness(tile_matrix)[[2]]
    feature_matrix[i,"separation_score"] <- calculate_clumpiness(tile_matrix)[[3]]
    
  }
  
  separation_scores <- feature_matrix[,3]
  
  
  # this cut-off method seemed too aggresive
  if(FALSE){
    
    separation_scores <- smooth.spline(1:length(separation_scores), sort(separation_scores), spar = 0.5)$y
    
    get_slope <- function(data_window, indices) {
      model <- lm(data_window ~ indices)
      return(coefficients(model)[2])  # Return the slope
    }
    
    # Create sliding windows and calculate slopes
    n            <- length(separation_scores)
    window_size  <- 10
    slopes       <- numeric(n - window_size + 1)  # Vector to store slopes
    bin_indices  <- 1:(n - window_size + 1)  # Starting index of each bin
    
    # Calculate slope for each window
    for(i in bin_indices) {
      window_data    <- separation_scores[i:(i + window_size - 1)]
      window_indices <- 1:window_size
      slopes[i]      <- get_slope(window_data, window_indices)
    }
    
    slopes              <- round(slopes,2)
    window_with_slope_1 <- which(slopes == min(slopes[slopes>=1]))
    window_of_interest  <- separation_scores[window_with_slope_1:(window_with_slope_1 + window_size - 1)]
    cutoff_value        <- window_of_interest[1]
    
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_classification-separation-scores_cis.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    plot(
      1:length(sort((abs(feature_matrix[,3])))),
      sort((abs(feature_matrix[,3]))),
      xlab = "sorted bedpe indices",
      ylab = "separation score",
      main = sample
    )
    abline(h=cutoff_value, col="red")
    dev.off()
  }
  
  
  if(TRUE){
    
    separation_scores <- log2(separation_scores)
    mean              <- mean(separation_scores)
    sd                <- sd(separation_scores)
    
    cutoff_value      <- mean + sd
    cutoff_value      <- 2^cutoff_value  
    
  }
  
  bedpe$sep_score <- as.numeric(feature_matrix[,"separation_score"])
  bedpe$max_dens  <- as.numeric(feature_matrix[,"max_density"])
  
  
  high_max_dens_bedpe <- bedpe[bedpe$max_dens >= 0.7,1:6]
  mid_max_dens_bedpe  <- bedpe[bedpe$max_dens >= 0.5 & 
                                 bedpe$max_dens < 0.7 &
                                 bedpe$sep_score >= cutoff_value, 1:6]
  
  gradient_bedpe <- rbind(
    high_max_dens_bedpe,
    mid_max_dens_bedpe
  )
  
  uniform_bedpe       <- setdiff(bedpe[,1:6], gradient_bedpe)
  
  if(nrow(uniform_bedpe)>0){
    uniform_bedpe$class <- "uniform"
  }
  
  pdf(
    file.path(out_dir, "plots", paste0(sample, "_classification-separation-scores_cis.pdf")),
    height = 8,
    width = 16,
    useDingbats = F
  )
  plot(
    feature_matrix[,1],
    feature_matrix[,3],
    xlab = "density of max density window",
    ylab = "separation score",
    main = sample
  )
  abline(h=(cutoff_value), col="red")
  abline(v=(0.5), col="red")
  dev.off()
  
  ##################################################################
  ##                      Floating Diagonals                      ##
  ##################################################################
  
  if(nrow(gradient_bedpe)>2){
    
    feature_matrix <- matrix(
      data = 0,
      nrow = nrow(gradient_bedpe),
      ncol = 1
    )
    colnames(feature_matrix) <- c("floating_score")
    
    
    for( i in 1:nrow(gradient_bedpe)){
      
      ##################################################################
      ##                    Extract contact matrix                    ##
      ##################################################################
      
      chr1    <- gradient_bedpe[i,1]
      start1  <- gradient_bedpe[i,2]
      end1    <- gradient_bedpe[i,3]
      chr2    <- gradient_bedpe[i,4]
      start2  <- gradient_bedpe[i,5]
      end2    <- gradient_bedpe[i,6]
      
      tile_contacts <- straw(
        "NONE",
        path_hic,
        paste( chr1, start1, end1, sep = ":"  ),
        paste( chr2, start2, end2, sep = ":"  ),
        "BP",
        5000)
      
      tile_matrix <- matrix(
        data = 0,
        nrow = length(seq(start1,end1,5000)),
        ncol = length(seq(start2,end2,5000))
      )
      
      rownames(tile_matrix) <- seq(start1,end1,5000)
      colnames(tile_matrix) <- seq(start2,end2,5000)
      
      for( j in 1:nrow(tile_contacts)){
        x <- as.character(tile_contacts[j,"x"])
        y <- as.character(tile_contacts[j,"y"])
        
        if(x %in% rownames(tile_matrix) && y %in% colnames(tile_matrix)){
          tile_matrix[
            as.character(tile_contacts[j,"x"]),
            as.character(tile_contacts[j,"y"])] <- tile_contacts[j,"counts"]
        } else if(x %in% colnames(tile_matrix) && y %in% rownames(tile_matrix)){
          tile_matrix[
            as.character(tile_contacts[j,"y"]),
            as.character(tile_contacts[j,"x"])] <- tile_contacts[j,"counts"]
        }
      }
      
      max_pos     <- which(tile_matrix == max(tile_matrix), arr.ind = TRUE)[1,]
      
      top_right_vec <- sapply(
        0:(min(max_pos[1], ncol(tile_matrix) - max_pos[2]) - 1), 
        function(i) tile_matrix[max_pos[1] - i, max_pos[2] + i])[-1]
      
      bottom_right_vec <- sapply(
        0:(min(nrow(tile_matrix) - max_pos[1], ncol(tile_matrix) - max_pos[2])), 
        function(i) tile_matrix[max_pos[1] + i, max_pos[2] + i])[-1]
      
      bottom_left_vec <- sapply(
        0:(min(nrow(tile_matrix) - max_pos[1], max_pos[2] - 1)), 
        function(i) tile_matrix[max_pos[1] + i, max_pos[2] - i])[-1]
      
      top_left_vec <- sapply(
        0:(min(max_pos[1], max_pos[2]) - 1), 
        function(i) tile_matrix[max_pos[1] - i, max_pos[2] - i])[-1]
      
      diagonal_1 <- sum(top_right_vec, bottom_left_vec)
      diagonal_2 <- sum(top_left_vec,  bottom_right_vec)
      
      diagonal_sum   <- sum(diagonal_1,diagonal_2)
      tile_sum       <- sum(tile_matrix)
      
      feature_matrix[i,"floating_score"] <- diagonal_sum/tile_sum
      
    }
    
    
    # cluster using k-means with 2 centers
    set.seed(123) 
    clusters         <- kmeans(feature_matrix[,"floating_score"], centers = 2)
    cluster_centers  <- clusters$centers
    diagonal_tiles   <- which.max(cluster_centers[,1])
    gradient_tiles   <- which.min(cluster_centers[,1])
    
    classifications  <- ifelse(
      clusters$cluster == diagonal_tiles, "floating_diagonal",
      ifelse(clusters$cluster == gradient_tiles, "gradient", "ambiguous"))
    
    gradient_bedpe$class <- classifications
    diagonal_bedpe       <- gradient_bedpe[gradient_bedpe$class == "floating_diagonal",]
    gradient_bedpe       <- gradient_bedpe[gradient_bedpe$class == "gradient",]
    
    if(nrow(diagonal_bedpe)>0){
      
      final_bedpe_cis <- rbind(
        gradient_bedpe,
        uniform_bedpe,
        diagonal_bedpe
      )
    } else {
      final_bedpe_cis <- rbind(
        gradient_bedpe,
        uniform_bedpe
      )
    }
    
  } else {
    final_bedpe_cis <- uniform_bedpe
  }
  
}




############################################################################
############################################################################
###                                                                      ###
###                         TRANS TRANSLOCATIONS                         ###
###                                                                      ###
############################################################################
############################################################################

if(flag_trans){
  
  bedpe <- whole_bedpe[whole_bedpe$chr1 != whole_bedpe$chr2,]
  
  ##################################################################
  ##                        Gradient tiles                        ##
  ##################################################################
  
  
  feature_matrix <- matrix(
    data = 0,
    nrow = nrow(bedpe),
    ncol = 3
  )
  colnames(feature_matrix) <- c("max_density","density_variance", "separation_score")
  
  
  for( i in 1:nrow(bedpe)){
    
    chr1    <- bedpe[i,1]
    start1  <- bedpe[i,2]
    end1    <- bedpe[i,3]
    chr2    <- bedpe[i,4]
    start2  <- bedpe[i,5]
    end2    <- bedpe[i,6]
    
    tile_contacts <- straw(
      "NONE",
      path_hic,
      paste( chr1, start1, end1, sep = ":"  ),
      paste( chr2, start2, end2, sep = ":"  ),
      "BP",
      5000)
    
    tile_matrix <- matrix(
      data = 0,
      nrow = length(seq(start1,end1,5000)),
      ncol = length(seq(start2,end2,5000))
    )
    
    rownames(tile_matrix) <- seq(start1,end1,5000)
    colnames(tile_matrix) <- seq(start2,end2,5000)
    
    for( j in 1:nrow(tile_contacts)){
      x <- as.character(tile_contacts[j,"x"])
      y <- as.character(tile_contacts[j,"y"])
      
      if(x %in% rownames(tile_matrix) && y %in% colnames(tile_matrix)){
        tile_matrix[
          as.character(tile_contacts[j,"x"]),
          as.character(tile_contacts[j,"y"])] <- tile_contacts[j,"counts"]
      } else if(x %in% colnames(tile_matrix) && y %in% rownames(tile_matrix)){
        tile_matrix[
          as.character(tile_contacts[j,"y"]),
          as.character(tile_contacts[j,"x"])] <- tile_contacts[j,"counts"]
      }
    }
    
    feature_matrix[i,"max_density"]      <- calculate_clumpiness(tile_matrix)[[1]]
    feature_matrix[i,"density_variance"] <- calculate_clumpiness(tile_matrix)[[2]]
    feature_matrix[i,"separation_score"] <- calculate_clumpiness(tile_matrix)[[3]]
    
  }
  
  separation_scores <- feature_matrix[,3]
  
  
  # this cut-off method seemed too aggresive
  if(FALSE){
    
    separation_scores <- smooth.spline(1:length(separation_scores), sort(separation_scores), spar = 0.5)$y
    
    get_slope <- function(data_window, indices) {
      model <- lm(data_window ~ indices)
      return(coefficients(model)[2])  # Return the slope
    }
    
    # Create sliding windows and calculate slopes
    n            <- length(separation_scores)
    window_size  <- 10
    slopes       <- numeric(n - window_size + 1)  # Vector to store slopes
    bin_indices  <- 1:(n - window_size + 1)  # Starting index of each bin
    
    # Calculate slope for each window
    for(i in bin_indices) {
      window_data    <- separation_scores[i:(i + window_size - 1)]
      window_indices <- 1:window_size
      slopes[i]      <- get_slope(window_data, window_indices)
    }
    
    slopes              <- round(slopes,2)
    window_with_slope_1 <- which(slopes == min(slopes[slopes>=1]))
    window_of_interest  <- separation_scores[window_with_slope_1:(window_with_slope_1 + window_size - 1)]
    cutoff_value        <- window_of_interest[1]
    
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_classification-separation-scores_cis_trans.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    plot(
      1:length(sort((abs(feature_matrix[,3])))),
      sort((abs(feature_matrix[,3]))),
      xlab = "sorted bedpe indices",
      ylab = "separation score",
      main = sample
    )
    abline(h=cutoff_value, col="red")
    dev.off()
  }
  
  
  if(TRUE){
    
    separation_scores <- log2(separation_scores)
    mean              <- mean(separation_scores)
    sd                <- sd(separation_scores)
    
    cutoff_value      <- mean + sd
    cutoff_value      <- 2^cutoff_value  
    
  }
  
  bedpe$sep_score <- as.numeric(feature_matrix[,"separation_score"])
  bedpe$max_dens  <- as.numeric(feature_matrix[,"max_density"])
  
  
  high_max_dens_bedpe <- bedpe[bedpe$max_dens >= 0.7,1:6]
  mid_max_dens_bedpe  <- bedpe[bedpe$max_dens >= 0.5 & 
                                 bedpe$max_dens < 0.7 &
                                 bedpe$sep_score >= cutoff_value, 1:6]
  
  gradient_bedpe <- rbind(
    high_max_dens_bedpe,
    mid_max_dens_bedpe
  )
  
  uniform_bedpe       <- setdiff(bedpe[,1:6], gradient_bedpe)
  if(nrow(uniform_bedpe)>0){
    uniform_bedpe$class <- "uniform"
  }
  
  pdf(
    file.path(out_dir, "plots", paste0(sample, "_classification-separation-scores_trans.pdf")),
    height = 8,
    width = 16,
    useDingbats = F
  )
  plot(
    feature_matrix[,1],
    feature_matrix[,3],
    xlab = "density of max density window",
    ylab = "separation score",
    main = sample
  )
  abline(h=(cutoff_value), col="red")
  abline(v=(0.5), col="red")
  dev.off()
  
  ##################################################################
  ##                      Floating Diagonals                      ##
  ##################################################################
  
  if(nrow(gradient_bedpe)>2){
    
    feature_matrix <- matrix(
      data = 0,
      nrow = nrow(gradient_bedpe),
      ncol = 1
    )
    colnames(feature_matrix) <- c("floating_score")
    
    
    for( i in 1:nrow(gradient_bedpe)){
      
      ##################################################################
      ##                    Extract contact matrix                    ##
      ##################################################################
      
      chr1    <- gradient_bedpe[i,1]
      start1  <- gradient_bedpe[i,2]
      end1    <- gradient_bedpe[i,3]
      chr2    <- gradient_bedpe[i,4]
      start2  <- gradient_bedpe[i,5]
      end2    <- gradient_bedpe[i,6]
      
      tile_contacts <- straw(
        "NONE",
        path_hic,
        paste( chr1, start1, end1, sep = ":"  ),
        paste( chr2, start2, end2, sep = ":"  ),
        "BP",
        5000)
      
      tile_matrix <- matrix(
        data = 0,
        nrow = length(seq(start1,end1,5000)),
        ncol = length(seq(start2,end2,5000))
      )
      
      rownames(tile_matrix) <- seq(start1,end1,5000)
      colnames(tile_matrix) <- seq(start2,end2,5000)
      
      for( j in 1:nrow(tile_contacts)){
        x <- as.character(tile_contacts[j,"x"])
        y <- as.character(tile_contacts[j,"y"])
        
        if(x %in% rownames(tile_matrix) && y %in% colnames(tile_matrix)){
          tile_matrix[
            as.character(tile_contacts[j,"x"]),
            as.character(tile_contacts[j,"y"])] <- tile_contacts[j,"counts"]
        } else if(x %in% colnames(tile_matrix) && y %in% rownames(tile_matrix)){
          tile_matrix[
            as.character(tile_contacts[j,"y"]),
            as.character(tile_contacts[j,"x"])] <- tile_contacts[j,"counts"]
        }
      }
      
      max_pos     <- which(tile_matrix == max(tile_matrix), arr.ind = TRUE)[1,]
      
      top_right_vec <- sapply(
        0:(min(max_pos[1], ncol(tile_matrix) - max_pos[2]) - 1), 
        function(i) tile_matrix[max_pos[1] - i, max_pos[2] + i])[-1]
      
      bottom_right_vec <- sapply(
        0:(min(nrow(tile_matrix) - max_pos[1], ncol(tile_matrix) - max_pos[2])), 
        function(i) tile_matrix[max_pos[1] + i, max_pos[2] + i])[-1]
      
      bottom_left_vec <- sapply(
        0:(min(nrow(tile_matrix) - max_pos[1], max_pos[2] - 1)), 
        function(i) tile_matrix[max_pos[1] + i, max_pos[2] - i])[-1]
      
      top_left_vec <- sapply(
        0:(min(max_pos[1], max_pos[2]) - 1), 
        function(i) tile_matrix[max_pos[1] - i, max_pos[2] - i])[-1]
      
      diagonal_1 <- sum(top_right_vec, bottom_left_vec)
      diagonal_2 <- sum(top_left_vec,  bottom_right_vec)
      
      diagonal_sum   <- sum(diagonal_1,diagonal_2)
      tile_sum       <- sum(tile_matrix)
      
      feature_matrix[i,"floating_score"] <- diagonal_sum/tile_sum
      
    }
    
    
    # cluster using k-means with 2 centers
    set.seed(123) 
    clusters         <- kmeans(feature_matrix[,"floating_score"], centers = 2)
    cluster_centers  <- clusters$centers
    diagonal_tiles   <- which.max(cluster_centers[,1])
    gradient_tiles   <- which.min(cluster_centers[,1])
    
    classifications  <- ifelse(
      clusters$cluster == diagonal_tiles, "floating_diagonal",
      ifelse(clusters$cluster == gradient_tiles, "gradient", "ambiguous"))
    
    gradient_bedpe$class <- classifications
    diagonal_bedpe       <- gradient_bedpe[gradient_bedpe$class == "floating_diagonal",]
    gradient_bedpe       <- gradient_bedpe[gradient_bedpe$class == "gradient",]
    
    if(nrow(diagonal_bedpe)>0){
      
      final_bedpe_trans <- rbind(
        gradient_bedpe,
        uniform_bedpe,
        diagonal_bedpe
      )
    } else {
      final_bedpe_trans <- rbind(
        gradient_bedpe,
        uniform_bedpe
      )
    }
    
  } else {
    final_bedpe_trans <- uniform_bedpe
  }
  
}


#################################################################
##                         Final bedpe                         ##
#################################################################

final_bedpe <- rbind(
  final_bedpe_cis,
  final_bedpe_trans
)


write.table(
  final_bedpe,
  file.path(out_dir, "results", paste0(sample, "_translocation-blocks_50KB_merged_fixed_classified.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

