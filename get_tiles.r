#!/usr/bin/env Rscript

options(scipen = 999)
suppressPackageStartupMessages(library(strawr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(parallel))
suppressPackageStartupMessages(library(doParallel))
suppressPackageStartupMessages(library(foreach))


args         <- commandArgs( trailingOnly = TRUE )
sample       <- as.character(args[1])
genome       <- as.character(args[2])
genome_size  <- as.character(args[3])
out_dir      <- as.character(args[4])


chr_sizes <- read.table(
  genome_size,
  as.is = T,
  col.names = c("chr","size")
)




#################################################################
##                          Functions                          ##
#################################################################

create_chr_bed <- function(chr, size, interval = tile_size) {
  starts <- seq(0, size - interval, by = interval)
  ends <- pmin(starts + interval, size)
  data.frame(
    chrom = rep(chr, length(starts)),
    start = starts,
    end = ends,
    stringsAsFactors = FALSE
  )
}

create_tiled_bedpe_parallel <- function(round1_tiles, tile_size) {
  # 1. Calculate total number of rows that will be generated
  n_tiles_left <- pmax(0, floor((round1_tiles$end1 - round1_tiles$start1) / tile_size))
  n_tiles_right <- pmax(0, floor((round1_tiles$end2 - round1_tiles$start2) / tile_size))
  rows_per_input <- n_tiles_left * n_tiles_right
  total_rows <- sum(rows_per_input)
  
  # 2. Set up parallel processing
  n_cores <- max(1, detectCores() - 1)
  chunks <- cut(1:nrow(round1_tiles), breaks = n_cores, labels = FALSE)
  
  # Create cluster and register it
  cl <- makeCluster(n_cores)
  registerDoParallel(cl)
  
  # 3. Process chunks in parallel
  results <- foreach(chunk_id = 1:n_cores, 
                     .combine = rbind,
                     .packages = c('data.table')) %dopar% {
                       # Get chunk indices
                       chunk_rows <- which(chunks == chunk_id)
                       chunk_data <- round1_tiles[chunk_rows, ]
                       
                       # Pre-calculate approximate size for this chunk
                       chunk_total_rows <- sum(rows_per_input[chunk_rows])
                       
                       # Initialize results for this chunk
                       chunk_results <- data.table(
                         chr.x = character(),
                         start.x = numeric(),
                         end.x = numeric(),
                         chr.y = character(),
                         start.y = numeric(),
                         end.y = numeric()
                       )
                       
                       # Process each row in the chunk
                       for(i in 1:nrow(chunk_data)) {
                         n_left <- n_tiles_left[chunk_rows[i]]
                         n_right <- n_tiles_right[chunk_rows[i]]
                         
                         if(n_left > 0 && n_right > 0) {
                           # Create sequences
                           left_starts <- seq(chunk_data$start1[i], 
                                              by = tile_size, 
                                              length.out = n_left)
                           right_starts <- seq(chunk_data$start2[i], 
                                               by = tile_size, 
                                               length.out = n_right)
                           
                           # Create all combinations using CJ (cross join) from data.table
                           grid <- CJ(left = left_starts, right = right_starts)
                           
                           # Create result rows
                           temp_result <- data.table(
                             chr.x = chunk_data$chr1[i],
                             start.x = grid$left,
                             end.x = grid$left + tile_size,
                             chr.y = chunk_data$chr2[i],
                             start.y = grid$right,
                             end.y = grid$right + tile_size
                           )
                           
                           # Combine with chunk results
                           chunk_results <- rbindlist(list(chunk_results, temp_result))
                         }
                       }
                       
                       chunk_results
                     }
  
  # 4. Clean up
  stopCluster(cl)
  
  # Convert back to data.frame if needed
  as.data.frame(results)
}



############################################################################
############################################################################
###                                                                      ###
###                          ROUND 1: 2MB TILES                          ###
###                                                                      ###
############################################################################
############################################################################

tile_size <- 2000000

cat("____Round-1: Tiling at 2Mb \n")

#################################################################
##                          Cis tiles                          ##
#################################################################

chr_beds <- list()
for (i in 1:nrow(chr_sizes)) {
  chr  <- chr_sizes$chr[i]
  size <- chr_sizes$size[i]
  
  chr_beds[[chr]] <- create_chr_bed(chr, size,tile_size)
}

cis_bedpe <- data.frame(
  chr1   = character(),
  start1 = numeric(),
  end1   = numeric(),
  stringsAsFactors = F
)

for(chr in names(chr_beds)){
  chr_df        <- chr_beds[[chr]]
  chr_cis_bedpe <- data.frame(
    chr1   = chr_df$chrom[-nrow(chr_df)],
    start1 = chr_df$start[-nrow(chr_df)],
    end1   = chr_df$end[-nrow(chr_df)],
    chr2   = chr_df$chrom[-1],
    start2 = chr_df$start[-1],
    end2   = chr_df$end[-1],
    stringsAsFactors = FALSE
  )
  
  cis_bedpe <- rbind(cis_bedpe,chr_cis_bedpe)
  rm(chr_df,chr_cis_bedpe)
}

#################################################################
##                         Trans tiles                         ##
#################################################################

trans_bedpe <- data.frame(
  chr1   = character(),
  start1 = numeric(),
  end1   = numeric(),
  chr2   = character(),
  start2 = numeric(),
  end2   = numeric(),
  stringsAsFactors = F
)

chr_names <- paste(
  "chr",
  c(1:22,"X","Y"),
  sep=""
)

for (i in 1:(length(chr_names) - 1)) {
  chr1 <- chr_names[i]
  chr1_df <- chr_beds[[chr1]]
  
  for (j in (i + 1):length(chr_names)) {
    chr2 <- chr_names[j]
    chr2_df <- chr_beds[[chr2]]
    
    trans_chr_bedpe <- merge(chr1_df, chr2_df, by = NULL)
    colnames(trans_chr_bedpe) <- c(
      "chr1", "start1", "end1",
      "chr2", "start2", "end2"
    )
    
    trans_bedpe <- rbind(trans_bedpe, trans_chr_bedpe)
    
    rm(trans_chr_bedpe)
  }
  rm(chr1_df,chr2_df)
}

# Add intra-chromosomal tiles beyond 3Mb
for(chr in chr_names) {
  chr_df <- chr_beds[[chr]]
  
  # Create all possible cis pairs
  cis_chr_bedpe <- merge(chr_df, chr_df, by = NULL)
  colnames(cis_chr_bedpe) <- c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2"
  )
  
  # Keep only pairs that are at least 3Mb apart
  cis_chr_bedpe <- cis_chr_bedpe[abs(cis_chr_bedpe$start2 - cis_chr_bedpe$start1) >= 3500000,]
  
  trans_bedpe <- rbind(trans_bedpe, cis_chr_bedpe)
  rm(cis_chr_bedpe)
}

# write cis tiles to file
write.table(
  cis_bedpe,
  file.path(out_dir, "intermediates", paste0(sample, "_cis.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

# write trans tiles to file
write.table(
  trans_bedpe,
  file.path(out_dir, "intermediates", paste0(sample, "_trans.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

# get cpm counts for cis tiles
command <- paste0(
  "bash ~/aqua_tools/annotate_loops.sh ",
  "-P ", file.path(out_dir, "intermediates", paste0(sample, "_cis.bedpe ")),
  "-A ", sample, " ",
  "-G ", genome, " ",
  "-Q cpm ",
  "--formula sum ",
  "-R 100000 > ", file.path(out_dir, "intermediates", paste0(sample, "_cis_2MB-tiles_100kb-res_cpm-counts.bedpe"))
)
system(command)

# get cpm counts for trans tiles
command <- paste0(
  "bash ~/aqua_tools/annotate_loops.sh ",
  "-P ", file.path(out_dir, "intermediates", paste0(sample, "_trans.bedpe ")),
  "-A ", sample, " ",
  "-G ", genome, " ",
  "-Q cpm ",
  "--formula sum ",
  "-R 100000 > ", file.path(out_dir, "intermediates", paste0(sample, "_trans_2MB-tiles_100kb-res_cpm-counts.bedpe"))
)
system(command)


cis_bedpe <- read.table(
  file.path(out_dir, "intermediates", paste0(sample, "_cis_2MB-tiles_100kb-res_cpm-counts.bedpe")),
  as.is = T,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "count"
  )
)

trans_bedpe <- read.table(
  file.path(out_dir, "intermediates", paste0(sample, "_trans_2MB-tiles_100kb-res_cpm-counts.bedpe")),
  as.is = T,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "count"
  )
)


##################################################################
##                         Thresholding                         ##
##################################################################

counts_cis <- cis_bedpe$count
counts_cis <- counts_cis[counts_cis != 0]
counts_cis <- log2(counts_cis)

if(sample %in% c("SJOS030605-X1_NT_H3K27ac")){
  threshold_inter  <- 0
  threshold_intra  <- 0
} else if(sample %in% c("SJOS063833-X1_NT_H3K27ac")) {
  threshold_inter  <- 0
  threshold_intra  <- 0
} else if(sample %in% c("RHB-P3F-463_NT_H3K27ac")) {
  threshold_inter  <- 2
  threshold_intra  <- 2
} else {
  threshold_inter <- mean(counts_cis) - 3*sd(counts_cis)
  threshold_intra <- mean(counts_cis) - 1*sd(counts_cis)
}


# trans tiles in trans space
trans_bedpe_inter <- trans_bedpe[trans_bedpe$chr1 != trans_bedpe$chr2,]
trans_bedpe_inter <- trans_bedpe_inter[trans_bedpe_inter$count > 0, ]
# trans tiles in cis space
trans_bedpe_intra <- trans_bedpe[trans_bedpe$chr1 == trans_bedpe$chr2,]
trans_bedpe_intra <- trans_bedpe_intra[trans_bedpe_intra$count > 0, ]

flag_plot <- TRUE
if(flag_plot){
  
  flag_intra <- TRUE
  flag_inter <- TRUE
  
  if(flag_intra){
    
    trans_bedpe <- trans_bedpe_intra
    
    set.seed(123)
    if(nrow(trans_bedpe) > nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe[sample(1:nrow(trans_bedpe),nrow(cis_bedpe),replace = F),]$count
      )
    } else if(nrow(trans_bedpe) < nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe[sample(1:nrow(cis_bedpe),nrow(trans_bedpe),replace = F),]$count,
        trans = trans_bedpe$count
      )
    } else {
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe$count
      )
    }
    
    plot_data <- plot_data[plot_data$cis   != 0,]
    plot_data <- plot_data[plot_data$trans != 0,]
    plot_data$cis   <- log2(plot_data$cis)
    plot_data$trans <- log2(plot_data$trans)
    plot_data_long  <- pivot_longer(plot_data, cols = c(cis, trans), names_to = "variable", values_to = "value")
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_2mb-distribution_cis.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    print(
      ggplot(plot_data_long, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.8) +
        geom_vline(xintercept = threshold_intra, color = "red", linetype = "dashed") +
        labs(
          x = "log2 sum (CPM) of tile",
          y = "Density") +
        theme_classic() +
        ggtitle(paste(sample,"intrachromosomal_tiles", sep = "-")) +
        scale_fill_manual(values = c("salmon", "lightblue"))
    )
    dev.off()
    
  }
  
  if(flag_inter){
    
    trans_bedpe <- trans_bedpe_inter
    
    set.seed(123)
    if(nrow(trans_bedpe) > nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe[sample(1:nrow(trans_bedpe),nrow(cis_bedpe),replace = F),]$count
      )
    } else if(nrow(trans_bedpe) < nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe[sample(1:nrow(cis_bedpe),nrow(trans_bedpe),replace = F),]$count,
        trans = trans_bedpe$count
      )
    } else {
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe$count
      )
    }
    
    plot_data       <- plot_data[plot_data$cis   != 0,]
    plot_data       <- plot_data[plot_data$trans != 0,]
    plot_data$cis   <- log2(plot_data$cis)
    plot_data$trans <- log2(plot_data$trans)
    plot_data_long  <- pivot_longer(plot_data, cols = c(cis, trans), names_to = "variable", values_to = "value")
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_2mb-distribution_trans.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    print(
      ggplot(plot_data_long, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.8) +
        geom_vline(xintercept = threshold_inter, color = "red", linetype = "dashed") +
        labs(
          x = "log2 sum (CPM) of tile",
          y = "Density") +
        theme_classic() +
        ggtitle(paste(sample,"interchromosomal_tiles", sep = "-")) +
        scale_fill_manual(values = c("salmon", "lightblue"))
    )
    dev.off()
    
  }
}


round1_tiles <- rbind(
  trans_bedpe_inter[log2(trans_bedpe_inter$count) >= threshold_inter,],
  trans_bedpe_intra[log2(trans_bedpe_intra$count) >= threshold_intra,]
)


write.table(
  round1_tiles,
  file.path(out_dir, "intermediates", paste0(sample, "_translocation-blocks_2MB.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

############################################################################
############################################################################
###                                                                      ###
###                         ROUND 2: 200KB TILES                         ###
###                                                                      ###
############################################################################
############################################################################

cat("____Round-2: Tiling at 200kb \n")

tile_size <- 200000


#################################################################
##                          Cis tiles                          ##
#################################################################

chr_beds <- list()
for (i in 1:nrow(chr_sizes)) {
  chr  <- chr_sizes$chr[i]
  size <- chr_sizes$size[i]
  
  chr_beds[[chr]] <- create_chr_bed(chr, size,tile_size)
}

cis_bedpe <- data.frame(
  chr1   = character(),
  start1 = numeric(),
  end1   = numeric(),
  stringsAsFactors = F
)

for(chr in names(chr_beds)){
  chr_df        <- chr_beds[[chr]]
  chr_cis_bedpe <- data.frame(
    chr1   = chr_df$chrom[-nrow(chr_df)],
    start1 = chr_df$start[-nrow(chr_df)],
    end1   = chr_df$end[-nrow(chr_df)],
    chr2   = chr_df$chrom[-1],
    start2 = chr_df$start[-1],
    end2   = chr_df$end[-1],
    stringsAsFactors = FALSE
  )
  
  cis_bedpe <- rbind(cis_bedpe,chr_cis_bedpe)
  rm(chr_df,chr_cis_bedpe)
}

#################################################################
##                         Trans tiles                         ##
#################################################################

trans_bedpe <- data.frame(
  chr1   = character(),
  start1 = numeric(),
  end1   = numeric(),
  chr2   = character(),
  start2 = numeric(),
  end2   = numeric(),
  stringsAsFactors = F
)

# making round-2 tiles using round-1 tiles as scaffold
trans_bedpe           <- create_tiled_bedpe_parallel(round1_tiles,tile_size)
colnames(trans_bedpe) <- c("chr1","start1","end1","chr2","start2","end2")


# write cis tiles to file
write.table(
  cis_bedpe,
  file.path(out_dir, "intermediates", paste0(sample, "_cis.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

# write trans tiles to file
write.table(
  trans_bedpe,
  file.path(out_dir, "intermediates", paste0(sample, "_trans.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

# get cpm counts for cis tiles
command <- paste0(
  "bash ~/aqua_tools/annotate_loops.sh ",
  "-P ", file.path(out_dir, "intermediates", paste0(sample, "_cis.bedpe ")),
  "-A ", sample, " ",
  "-G ", genome, " ",
  "-Q cpm ",
  "--formula sum ",
  "-R 5000 > ", file.path(out_dir, "intermediates", paste0(sample, "_cis_200KB-tiles_5kb-res_cpm-counts.bedpe"))
)
system(command)

# get cpm counts for trans tiles
command <- paste0(
  "bash ~/aqua_tools/annotate_loops.sh ",
  "-P ", file.path(out_dir, "intermediates", paste0(sample, "_trans.bedpe ")),
  "-A ", sample, " ",
  "-G ", genome, " ",
  "-Q cpm ",
  "--formula sum ",
  "-R 5000 > ", file.path(out_dir, "intermediates", paste0(sample, "_trans_200KB-tiles_5kb-res_cpm-counts.bedpe"))
)
system(command)


cis_bedpe <- read.table(
  file.path(out_dir, "intermediates", paste0(sample, "_cis_200KB-tiles_5kb-res_cpm-counts.bedpe")),
  as.is = T,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "count"
  )
)

trans_bedpe <- read.table(
  file.path(out_dir, "intermediates", paste0(sample, "_trans_200KB-tiles_5kb-res_cpm-counts.bedpe")),
  as.is = T,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "count"
  )
)


##################################################################
##                         Thresholding                         ##
##################################################################

counts_cis <- cis_bedpe$count
counts_cis <- counts_cis[counts_cis != 0]
counts_cis <- log2(counts_cis)


if(sample %in% c("SJOS030605-X1_NT_H3K27ac")){
  threshold_inter  <- -4
  threshold_intra  <- -4
} else if(sample %in% c("SJOS063833-X1_NT_H3K27ac")) {
  threshold_inter  <- -4
  threshold_intra  <- -4
} else if(sample %in% c("RHB-P3F-463_NT_H3K27ac")) {
  threshold_inter  <- -2
  threshold_intra  <- -2
} else {
  threshold_inter <- mean(counts_cis) - 3*sd(counts_cis)
  threshold_intra <- mean(counts_cis) - 1*sd(counts_cis)
}


# trans tiles in trans space
trans_bedpe_inter <- trans_bedpe[trans_bedpe$chr1 != trans_bedpe$chr2,]
trans_bedpe_inter <- trans_bedpe_inter[trans_bedpe_inter$count > 0, ]
# trans tiles in cis space
trans_bedpe_intra <- trans_bedpe[trans_bedpe$chr1 == trans_bedpe$chr2,]
trans_bedpe_intra <- trans_bedpe_intra[trans_bedpe_intra$count > 0, ]


flag_plot <- TRUE
if(flag_plot){
  
  flag_intra <- TRUE
  flag_inter <- TRUE
  
  if(flag_intra){
    
    trans_bedpe <- trans_bedpe_intra
    
    set.seed(123)
    if(nrow(trans_bedpe) > nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe[sample(1:nrow(trans_bedpe),nrow(cis_bedpe),replace = F),]$count
      )
    } else if(nrow(trans_bedpe) < nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe[sample(1:nrow(cis_bedpe),nrow(trans_bedpe),replace = F),]$count,
        trans = trans_bedpe$count
      )
    } else {
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe$count
      )
    }
    
    plot_data <- plot_data[plot_data$cis   != 0,]
    plot_data <- plot_data[plot_data$trans != 0,]
    plot_data$cis   <- log2(plot_data$cis)
    plot_data$trans <- log2(plot_data$trans)
    plot_data_long  <- pivot_longer(plot_data, cols = c(cis, trans), names_to = "variable", values_to = "value")
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_200kb-distribution_cis.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    print(
      ggplot(plot_data_long, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.8) +
        geom_vline(xintercept = threshold_intra, color = "red", linetype = "dashed") +
        labs(
          x = "log2 sum (CPM) of tile",
          y = "Density") +
        theme_classic() +
        ggtitle(paste(sample,"intrachromosomal_tiles", sep = "-")) +
        scale_fill_manual(values = c("salmon", "lightblue"))
    )
    dev.off()
    
  }
  
  if(flag_inter){
    
    trans_bedpe <- trans_bedpe_inter
    
    set.seed(123)
    if(nrow(trans_bedpe) > nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe[sample(1:nrow(trans_bedpe),nrow(cis_bedpe),replace = F),]$count
      )
    } else if(nrow(trans_bedpe) < nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe[sample(1:nrow(cis_bedpe),nrow(trans_bedpe),replace = F),]$count,
        trans = trans_bedpe$count
      )
    } else {
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe$count
      )
    }
    
    plot_data       <- plot_data[plot_data$cis   != 0,]
    plot_data       <- plot_data[plot_data$trans != 0,]
    plot_data$cis   <- log2(plot_data$cis)
    plot_data$trans <- log2(plot_data$trans)
    plot_data_long  <- pivot_longer(plot_data, cols = c(cis, trans), names_to = "variable", values_to = "value")
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_200kb-distribution_trans.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    print(
      ggplot(plot_data_long, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.8) +
        geom_vline(xintercept = threshold_inter, color = "red", linetype = "dashed") +
        labs(
          x = "log2 sum (CPM) of tile",
          y = "Density") +
        theme_classic() +
        ggtitle(paste(sample,"interchromosomal_tiles", sep = "-")) +
        scale_fill_manual(values = c("salmon", "lightblue"))
    )
    dev.off()
    
  }
}

round2_tiles <- rbind(
  trans_bedpe_inter[log2(trans_bedpe_inter$count) >= threshold_inter,],
  trans_bedpe_intra[log2(trans_bedpe_intra$count) >= threshold_intra,]
)

write.table(
  round2_tiles,
  file.path(out_dir, "intermediates", paste0(sample, "_translocation-blocks_200KB.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

############################################################################
############################################################################
###                                                                      ###
###                          ROUND 3: 50KB TILES                         ###
###                                                                      ###
############################################################################
############################################################################

cat("____Round-3: Tiling at 50kb \n")

tile_size <- 50000


#################################################################
##                          Cis tiles                          ##
#################################################################

chr_beds <- list()
for (i in 1:nrow(chr_sizes)) {
  chr  <- chr_sizes$chr[i]
  size <- chr_sizes$size[i]
  
  chr_beds[[chr]] <- create_chr_bed(chr, size,tile_size)
}

cis_bedpe <- data.frame(
  chr1   = character(),
  start1 = numeric(),
  end1   = numeric(),
  stringsAsFactors = F
)

for(chr in names(chr_beds)){
  chr_df        <- chr_beds[[chr]]
  chr_cis_bedpe <- data.frame(
    chr1   = chr_df$chrom[-nrow(chr_df)],
    start1 = chr_df$start[-nrow(chr_df)],
    end1   = chr_df$end[-nrow(chr_df)],
    chr2   = chr_df$chrom[-1],
    start2 = chr_df$start[-1],
    end2   = chr_df$end[-1],
    stringsAsFactors = FALSE
  )
  
  cis_bedpe <- rbind(cis_bedpe,chr_cis_bedpe)
  rm(chr_df,chr_cis_bedpe)
}

#################################################################
##                         Trans tiles                         ##
#################################################################

trans_bedpe <- data.frame(
  chr1   = character(),
  start1 = numeric(),
  end1   = numeric(),
  chr2   = character(),
  start2 = numeric(),
  end2   = numeric(),
  stringsAsFactors = F
)

# making round-3 tiles using round-2 tiles as scaffold
trans_bedpe           <- create_tiled_bedpe_parallel(round2_tiles,tile_size)
colnames(trans_bedpe) <- c("chr1","start1","end1","chr2","start2","end2")

# write cis tiles to file
write.table(
  cis_bedpe,
  file.path(out_dir, "intermediates", paste0(sample, "_cis.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

# write trans tiles to file
write.table(
  trans_bedpe,
  file.path(out_dir, "intermediates", paste0(sample, "_trans.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

# get cpm counts for cis tiles
command <- paste0(
  "bash ~/aqua_tools/annotate_loops.sh ",
  "-P ", file.path(out_dir, "intermediates", paste0(sample, "_cis.bedpe ")),
  "-A ", sample, " ",
  "-G ", genome, " ",
  "-Q cpm ",
  "--formula sum ",
  "-R 5000 > ", file.path(out_dir, "intermediates", paste0(sample, "_cis_50KB-tiles_5kb-res_cpm-counts.bedpe"))
)
system(command)

# get cpm counts for trans tiles
command <- paste0(
  "bash ~/aqua_tools/annotate_loops.sh ",
  "-P ", file.path(out_dir, "intermediates", paste0(sample, "_trans.bedpe ")),
  "-A ", sample, " ",
  "-G ", genome, " ",
  "-Q cpm ",
  "--formula sum ",
  "-R 5000 > ", file.path(out_dir, "intermediates", paste0(sample, "_trans_50KB-tiles_5kb-res_cpm-counts.bedpe"))
)
system(command)


cis_bedpe <- read.table(
  file.path(out_dir, "intermediates", paste0(sample, "_cis_50KB-tiles_5kb-res_cpm-counts.bedpe")),
  as.is = T,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "count"
  )
)

trans_bedpe <- read.table(
  file.path(out_dir, "intermediates", paste0(sample, "_trans_50KB-tiles_5kb-res_cpm-counts.bedpe")),
  as.is = T,
  col.names = c(
    "chr1", "start1", "end1",
    "chr2", "start2", "end2",
    "count"
  )
)


##################################################################
##                         Thresholding                         ##
##################################################################

counts_cis <- cis_bedpe$count
counts_cis <- counts_cis[counts_cis != 0]
counts_cis <- log2(counts_cis)

if(sample %in% c("SJOS030605-X1_NT_H3K27ac")){
  threshold_inter  <- -4
  threshold_intra  <- -4
} else if(sample %in% c("SJOS063833-X1_NT_H3K27ac")) {
  threshold_inter  <- -4
  threshold_intra  <- -4
} else if(sample %in% c("RHB-P3F-463_NT_H3K27ac")) {
  threshold_inter  <- -3
  threshold_intra  <- -3
}else {
  threshold_inter <- mean(counts_cis) - 3*sd(counts_cis)
  threshold_intra <- mean(counts_cis) - 1*sd(counts_cis)
}


# trans tiles in trans space
trans_bedpe_inter <- trans_bedpe[trans_bedpe$chr1 != trans_bedpe$chr2,]
trans_bedpe_inter <- trans_bedpe_inter[trans_bedpe_inter$count > 0, ]
# trans tiles in cis space
trans_bedpe_intra <- trans_bedpe[trans_bedpe$chr1 == trans_bedpe$chr2,]
trans_bedpe_intra <- trans_bedpe_intra[trans_bedpe_intra$count > 0, ]


flag_plot <- TRUE
if(flag_plot){
  
  flag_intra <- TRUE
  flag_inter <- TRUE
  
  if(flag_intra){
    
    trans_bedpe <- trans_bedpe_intra
    
    set.seed(123)
    if(nrow(trans_bedpe) > nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe[sample(1:nrow(trans_bedpe),nrow(cis_bedpe),replace = F),]$count
      )
    } else if(nrow(trans_bedpe) < nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe[sample(1:nrow(cis_bedpe),nrow(trans_bedpe),replace = F),]$count,
        trans = trans_bedpe$count
      )
    } else {
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe$count
      )
    }
    
    plot_data <- plot_data[plot_data$cis   != 0,]
    plot_data <- plot_data[plot_data$trans != 0,]
    plot_data$cis   <- log2(plot_data$cis)
    plot_data$trans <- log2(plot_data$trans)
    plot_data_long  <- pivot_longer(plot_data, cols = c(cis, trans), names_to = "variable", values_to = "value")
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_50kb-distribution_cis.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    print(
      ggplot(plot_data_long, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.8) +
        geom_vline(xintercept = threshold_intra, color = "red", linetype = "dashed") +
        labs(
          x = "log2 sum (CPM) of tile",
          y = "Density") +
        theme_classic() +
        ggtitle(paste(sample,"intrachromosomal_tiles", sep = "-")) +
        scale_fill_manual(values = c("salmon", "lightblue"))
    )
    dev.off()
    
  }
  
  if(flag_inter){
    
    trans_bedpe <- trans_bedpe_inter
    
    set.seed(123)
    if(nrow(trans_bedpe) > nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe[sample(1:nrow(trans_bedpe),nrow(cis_bedpe),replace = F),]$count
      )
    } else if(nrow(trans_bedpe) < nrow(cis_bedpe)){
      plot_data <- data.frame(
        cis   = cis_bedpe[sample(1:nrow(cis_bedpe),nrow(trans_bedpe),replace = F),]$count,
        trans = trans_bedpe$count
      )
    } else {
      plot_data <- data.frame(
        cis   = cis_bedpe$count,
        trans = trans_bedpe$count
      )
    }
    
    plot_data       <- plot_data[plot_data$cis   != 0,]
    plot_data       <- plot_data[plot_data$trans != 0,]
    plot_data$cis   <- log2(plot_data$cis)
    plot_data$trans <- log2(plot_data$trans)
    plot_data_long  <- pivot_longer(plot_data, cols = c(cis, trans), names_to = "variable", values_to = "value")
    
    pdf(
      file.path(out_dir, "plots", paste0(sample, "_50kb-distribution_trans.pdf")),
      height = 8,
      width = 8,
      useDingbats = F
    )
    print(
      ggplot(plot_data_long, aes(x = value, fill = variable)) +
        geom_density(alpha = 0.8) +
        geom_vline(xintercept = threshold_inter, color = "red", linetype = "dashed") +
        labs(
          x = "log2 sum (CPM) of tile",
          y = "Density") +
        theme_classic() +
        ggtitle(paste(sample,"interchromosomal_tiles", sep = "-")) +
        scale_fill_manual(values = c("salmon", "lightblue"))
    )
    dev.off()
    
  }
}

round3_tiles <- rbind(
  trans_bedpe_inter[log2(trans_bedpe_inter$count) >= threshold_inter,],
  trans_bedpe_intra[log2(trans_bedpe_intra$count) >= threshold_intra,]
)


write.table(
  round3_tiles,
  file.path(out_dir, "intermediates", paste0(sample, "_translocation-blocks_50KB.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)

##################################################################
##                      Merging 50kb tiles                      ##
##################################################################

flag_cis   <- TRUE
flag_trans <- TRUE


if(flag_trans){
  
  chromosomes <- paste0("chr", c(1:22, "X", "Y"))
  
  chr1 <- character()
  chr2 <- character()
  
  for (i in 1:(length(chromosomes) - 1)) {
    for (j in (i + 1):length(chromosomes)) {
      chr1 <- c(chr1, chromosomes[i])
      chr2 <- c(chr2, chromosomes[j])
    }
  }
  chr_pairs <- data.frame(chr1 = chr1, chr2 = chr2)
  
  
  trans_tiles_merged <- data.frame()
  
  for(i in 1:nrow(chr_pairs)){
    
    trans_bedpe <- round3_tiles[
      round3_tiles$chr1 == chr_pairs[i,1] & round3_tiles$chr2 == chr_pairs[i,2],]
    
    if(nrow(trans_bedpe) > 0){
      
      trans_bedpe <- trans_bedpe[order(
        trans_bedpe$chr1,trans_bedpe$start1,
        trans_bedpe$chr2,trans_bedpe$start2),]
      rownames(trans_bedpe) <- 1:nrow(trans_bedpe)
      
      chra_size <- chr_sizes[chr_sizes$chr == chr_pairs[i,1],2]
      chra_bins <- seq(0,chra_size,by=50000)
      
      chrb_size <- chr_sizes[chr_sizes$chr == chr_pairs[i,2],2]
      chrb_bins <- seq(0,chrb_size,by=50000)
      
      
      positions <- matrix(
        data=NA,
        nrow=nrow(trans_bedpe),
        ncol=2)
      colnames(positions) <- c("row","col")
      
      
      for(j in 1:nrow(trans_bedpe)){
        
        row_index <- which(trans_bedpe[j,"start1"] == chra_bins)
        col_index <- which(trans_bedpe[j,"start2"] == chrb_bins)
        
        positions[j,"row"] <- row_index
        positions[j,"col"] <- col_index
        
      }
      
      clust               <- dbscan(positions, eps = 2, minPts = 2)
      trans_bedpe$cluster <- clust$cluster
      
      trans_bedpe        <- trans_bedpe[order(trans_bedpe$cluster),]
      trans_bedpe_blocks <- trans_bedpe[trans_bedpe$cluster > 0,]
      
      unique_blocks <- unique(trans_bedpe_blocks$cluster)
      
      merged_bedpe <- data.frame()
      for(block in unique_blocks){
        
        df   <- trans_bedpe_blocks[trans_bedpe_blocks$cluster == block,]
        df_2 <- data.frame(
          chr1   = df[1,"chr1"],
          start1 = min(df[,"start1"]),
          end1   = max(df[,"end1"]),
          chr2   = df[1,"chr2"],
          start2 = min(df[,"start2"]),
          end2   = max(df[,"end2"])
        )
        
        merged_bedpe <- rbind(merged_bedpe,df_2)
        rm(df,df_2)
        
      }
      
      merged_bedpe <- rbind(
        merged_bedpe,
        trans_bedpe[trans_bedpe$cluster == 0,1:6]
      )
      
      merged_bedpe <- merged_bedpe[order(
        merged_bedpe$chr1,merged_bedpe$start1,
        merged_bedpe$chr2,merged_bedpe$start2),]
      
      trans_tiles_merged <- rbind(
        trans_tiles_merged,
        merged_bedpe
      )
      
    } else {
      next
    }
    
  }
  
}

if(flag_cis){
  
  chromosomes <- paste0("chr", c(1:22, "X", "Y"))
  
  chr_pairs   <- data.frame(chr1 = chromosomes, chr2 = chromosomes)
  
  
  cis_tiles_merged <- data.frame()
  
  for(i in 1:nrow(chr_pairs)){
    
    trans_bedpe <- round3_tiles[
      round3_tiles$chr1 == chr_pairs[i,1] & round3_tiles$chr2 == chr_pairs[i,2],]
    
    if(nrow(trans_bedpe) > 0){
      
      trans_bedpe <- trans_bedpe[order(
        trans_bedpe$chr1,trans_bedpe$start1,
        trans_bedpe$chr2,trans_bedpe$start2),]
      rownames(trans_bedpe) <- 1:nrow(trans_bedpe)
      
      chra_size <- chr_sizes[chr_sizes$chr == chr_pairs[i,1],2]
      chra_bins <- seq(0,chra_size,by=50000)
      
      chrb_size <- chr_sizes[chr_sizes$chr == chr_pairs[i,2],2]
      chrb_bins <- seq(0,chrb_size,by=50000)
      
      
      positions <- matrix(
        data=NA,
        nrow=nrow(trans_bedpe),
        ncol=2)
      colnames(positions) <- c("row","col")
      
      
      for(j in 1:nrow(trans_bedpe)){
        
        row_index <- which(trans_bedpe[j,"start1"] == chra_bins)
        col_index <- which(trans_bedpe[j,"start2"] == chrb_bins)
        
        positions[j,"row"] <- row_index
        positions[j,"col"] <- col_index
        
      }
      
      clust               <- dbscan(positions, eps = 2, minPts = 2)
      trans_bedpe$cluster <- clust$cluster
      
      trans_bedpe        <- trans_bedpe[order(trans_bedpe$cluster),]
      trans_bedpe_blocks <- trans_bedpe[trans_bedpe$cluster > 0,]
      
      unique_blocks <- unique(trans_bedpe_blocks$cluster)
      
      merged_bedpe <- data.frame()
      for(block in unique_blocks){
        
        df   <- trans_bedpe_blocks[trans_bedpe_blocks$cluster == block,]
        df_2 <- data.frame(
          chr1   = df[1,"chr1"],
          start1 = min(df[,"start1"]),
          end1   = max(df[,"end1"]),
          chr2   = df[1,"chr2"],
          start2 = min(df[,"start2"]),
          end2   = max(df[,"end2"])
        )
        
        merged_bedpe <- rbind(merged_bedpe,df_2)
        rm(df,df_2)
        
      }
      
      merged_bedpe <- rbind(
        merged_bedpe,
        trans_bedpe[trans_bedpe$cluster == 0,1:6]
      )
      
      merged_bedpe <- merged_bedpe[order(
        merged_bedpe$chr1,merged_bedpe$start1,
        merged_bedpe$chr2,merged_bedpe$start2),]
      
      cis_tiles_merged <- rbind(
        cis_tiles_merged,
        merged_bedpe
      )
      
    } else {
      next
    }
    
  }
  
}


write.table(
  rbind(
    trans_tiles_merged,
    cis_tiles_merged
  ),
  file.path(out_dir, "results", paste0(sample, "_translocation-blocks_50KB_merged.bedpe")),
  row.names = F,
  col.names = F,
  quote = F,
  sep = "\t"
)


