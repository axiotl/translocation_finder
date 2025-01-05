#!/usr/bin/env Rscript

options(scipen = 999)
suppressPackageStartupMessages(library(strawr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(dbscan))
suppressPackageStartupMessages(library(parallel))


args    <- commandArgs( trailingOnly = TRUE )
sample  <- as.character(args[1])
genome  <- as.character(args[2])
bedpe   <- as.character(args[3])
out_dir <- as.character(args[4])



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

# Set up parallel processing
num_cores <- detectCores() - 1
search_space <- 100000

# Calculate number of rows and chunk size
total_rows <- length(readLines(bedpe))
chunk_size <- 1000  

# Calculate number of chunks
num_chunks <- ceiling(total_rows / chunk_size)

# Open input and output connections
con_in   <- file(description = bedpe, open = "r")
out_file <- file.path(out_dir, "results", paste0(sample, "_translocation-blocks_50KB_merged_fixed.bedpe"))
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)


process_tile <- function(tile, path_hic, search_space) {
  
  if(tile[,"start2"] < tile[,"start1"]) {
    tile_2 <- cbind(tile[,4:6], tile[,1:3])
    tile <- tile_2
    rm(tile_2)
  }
  
  
  tile_contacts <- straw(
    "NONE",
    path_hic,
    paste(tile[1,1], tile[1,2], tile[1,3], sep = ":"),
    paste(tile[1,4], tile[1,5], tile[1,6], sep = ":"),
    "BP",
    5000
  )
  invisible(gc())
  
  ##################################################################
  ##                    Expand left coordinate                    ##
  ##################################################################
  left           <- tile 
  left$start1    <- left$start1 - search_space
  
  left_contacts  <- straw(
    "NONE",
    path_hic,
    paste( left[1,1], left[1,2], left[1,3], sep = ":"  ),
    paste( left[1,4], left[1,5], left[1,6], sep = ":"  ),
    "BP",
    5000)
  invisible(gc())
  
  left_contacts <- anti_join(
    left_contacts, 
    tile_contacts, 
    by = c("x", "y"))
  
  new_left_coordinate <- tile$start1
  if(nrow(left_contacts)>0){
    if(max(left_contacts$counts) > max(tile_contacts$counts)){
      
      id <- which(left_contacts$counts > max(tile_contacts$counts))
      df <- left_contacts[id,]
      df <- df[order(df$counts),]
      
      # because i'm not sure if the straw output coordinate of interest
      # will be in column x or column y
      x  <- abs(df[nrow(df),"x"] - new_left_coordinate)
      y  <- abs(df[nrow(df),"y"] - new_left_coordinate)
      
      new_left_coordinate <- df[nrow(df),which.min(c(x,y))] - 5000
      
    }
  }
  
  #################################################################
  ##                   Expand right coordinate                   ##
  #################################################################
  right           <- tile
  right$end1      <- right$end1 + search_space
  
  right_contacts  <- straw(
    "NONE",
    path_hic,
    paste( right[1,1], right[1,2], right[1,3], sep = ":"  ),
    paste( right[1,4], right[1,5], right[1,6], sep = ":"  ),
    "BP",
    5000)
  invisible(gc())
  
  right_contacts <- anti_join(
    right_contacts, 
    tile_contacts, 
    by = c("x", "y"))
  
  new_right_coordinate <- tile$end1
  if(nrow(right_contacts)>0){
    if(max(right_contacts$counts) > max(tile_contacts$counts)){
      
      id <- which(right_contacts$counts > max(tile_contacts$counts))
      df <- right_contacts[id,]
      df <- df[order(df$counts),]
      
      # because i'm not sure if the straw output coordinate of interest
      # will be in column x or column y
      x  <- abs(df[nrow(df),"x"] - new_right_coordinate)
      y  <- abs(df[nrow(df),"y"] - new_right_coordinate)
      
      new_right_coordinate <- df[nrow(df),which.min(c(x,y))] + 5000
      
    }
  }
  
  ##################################################################
  ##                     Expand up coordinate                     ##
  ##################################################################
  up           <- tile
  up$start2    <- up$start2 - search_space
  
  up_contacts  <- straw(
    "NONE",
    path_hic,
    paste( up[1,1], up[1,2], up[1,3], sep = ":"  ),
    paste( up[1,4], up[1,5], up[1,6], sep = ":"  ),
    "BP",
    5000)
  invisible(gc())
  
  up_contacts <- anti_join(
    up_contacts, 
    tile_contacts, 
    by = c("x", "y"))
  
  new_up_coordinate <- tile$start2
  if(nrow(up_contacts)>0){
    if(max(up_contacts$counts) > max(tile_contacts$counts)){
      
      id <- which(up_contacts$counts > max(tile_contacts$counts))
      df <- up_contacts[id,]
      df <- df[order(df$counts),]
      
      # because i'm not sure if the straw output coordinate of interest
      # will be in column x or column y
      x  <- abs(df[nrow(df),"x"] - new_up_coordinate)
      y  <- abs(df[nrow(df),"y"] - new_up_coordinate)
      
      new_up_coordinate <- df[nrow(df),which.min(c(x,y))] - 5000
      
    }
  }
  
  ##################################################################
  ##                    Expand down coordinate                    ##
  ##################################################################
  down           <- tile
  down$end2      <- down$end2 + search_space
  
  down_contacts  <- straw(
    "NONE",
    path_hic,
    paste( down[1,1], down[1,2], down[1,3], sep = ":"  ),
    paste( down[1,4], down[1,5], down[1,6], sep = ":"  ),
    "BP",
    5000)
  invisible(gc())
  
  down_contacts <- anti_join(
    down_contacts, 
    tile_contacts, 
    by = c("x", "y"))
  
  new_down_coordinate <- tile$end2
  if(nrow(down_contacts)>0){
    if(max(down_contacts$counts) > max(tile_contacts$counts)){
      
      id <- which(down_contacts$counts > max(tile_contacts$counts))
      df <- down_contacts[id,]
      df <- df[order(df$counts),]
      
      # because i'm not sure if the straw output coordinate of interest
      # will be in column x or column y
      x  <- abs(df[nrow(df),"x"] - new_down_coordinate)
      y  <- abs(df[nrow(df),"y"] - new_down_coordinate)
      
      new_down_coordinate <- df[nrow(df),which.min(c(x,y))] + 5000
      
    }
  }
  
  ##################################################################
  ##                        Fixing borders                        ##
  ##################################################################
  
  new_tile <- data.frame(
    chr1   = tile$chr1,
    start1 = new_left_coordinate,
    end1   = new_right_coordinate,
    chr2   = tile$chr2,
    start2 = new_up_coordinate,
    end2   = new_down_coordinate
  )
  
  new_tile_contacts <- straw(
    "NONE",
    path_hic,
    paste( new_tile[1,1], new_tile[1,2], new_tile[1,3], sep = ":"  ),
    paste( new_tile[1,4], new_tile[1,5], new_tile[1,6], sep = ":"  ),
    "BP",
    5000)
  invisible(gc())
  
  new_tile_matrix <- matrix(
    data = 0,
    nrow = length(seq(new_tile[1,2],new_tile[1,3],5000)),
    ncol = length(seq(new_tile[1,5],new_tile[1,6],5000))
  )
  
  rownames(new_tile_matrix) <- seq(new_tile[1,2],new_tile[1,3],5000)
  colnames(new_tile_matrix) <- seq(new_tile[1,5],new_tile[1,6],5000)
  
  for( j in 1:nrow(new_tile_contacts)){
    
    x <- as.character(new_tile_contacts[j,"x"])
    y <- as.character(new_tile_contacts[j,"y"])
    
    if(x %in% rownames(new_tile_matrix) && y %in% colnames(new_tile_matrix)){
      new_tile_matrix[
        as.character(new_tile_contacts[j,"x"]),
        as.character(new_tile_contacts[j,"y"])] <- new_tile_contacts[j,"counts"]
    } else if(x %in% colnames(new_tile_matrix) && y %in% rownames(new_tile_matrix)){
      new_tile_matrix[
        as.character(new_tile_contacts[j,"y"]),
        as.character(new_tile_contacts[j,"x"])] <- new_tile_contacts[j,"counts"]
    }
  }
  
  x <- apply(new_tile_matrix,1,sum)
  y <- apply(new_tile_matrix,2,sum)
  
  x <- x[x!=0]
  y <- y[y!=0]
  
  new_tile_matrix <- new_tile_matrix[names(x),names(y), drop = FALSE]
  
  new_tile <- data.frame(
    chr1   = tile$chr1,
    start1 = as.numeric(rownames(new_tile_matrix)[1]),
    end1   = as.numeric(rownames(new_tile_matrix)[nrow(new_tile_matrix)]) + 5000,
    chr2   = tile$chr2,
    start2 = as.numeric(colnames(new_tile_matrix)[1]),
    end2   = as.numeric(colnames(new_tile_matrix)[ncol(new_tile_matrix)]) + 5000
  )
  return(new_tile)
}


# Process chunks in a loop
index <- 0
repeat {
  index <- index + 1
  #cat(sprintf("Processing chunk %d of %d\n", index, num_chunks))
  
  # Read chunk
  current_chunk <- try(read.table(
    con_in,
    nrows = chunk_size,
    col.names = c("chr1", "start1", "end1", "chr2", "start2", "end2")
  ))
  
  # Check if we've reached the end of file
  if (inherits(current_chunk, "try-error") || nrow(current_chunk) == 0) {
    break
  }
  
  # Process chunk in parallel
  results <- mclapply(1:nrow(current_chunk), function(i) {
    tile <- current_chunk[i,]
    tryCatch({
      process_tile(tile, path_hic, search_space)
    }, error = function(e) {
      cat(sprintf("Error processing row %d: %s\n", i, e$message))
      return(NULL)
    })
  }, mc.cores = num_cores)
  
  # Combine results
  valid_results <- do.call(rbind, results[!sapply(results, is.null)])
  
  # Write results
  if (nrow(valid_results) > 0) {
    write.table(
      valid_results,
      file = out_file,
      append = TRUE,
      quote = FALSE,
      sep = "\t",
      row.names = FALSE,
      col.names = FALSE
    )
  }
  
  # Clean up memory
  rm(current_chunk, results, valid_results)
  invisible(gc())
}

# Close connections
close(con_in)


