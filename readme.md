# translocation_finder

 **Table of contents:**

 - [translocation_finder.sh](#translocation_finder.sh)
 - [get_tiles.r](#get_tiles.r)
 - [fix_tiles.r](#fix_tiles.r)
 - [classify_tiles.r](#classify_tiles.r)


 <a id="translocation_finder.sh"></a>
## translocation_finder.sh

**Summary:** Outputs bedpe files (double bedfile) containing translocations using H3K27ac HiChIP .hic files.

**Process:**
* 1-115: Arguments and help
    - `$A` = Sample name
    - `$G` = Genome build (easily checked with command `list_samples`)
    - `$S` = Chromosome sizes (2 columns: chromosome and size)
    - `$O` = Out path 
> Help/usage descriptions do not correspond with code!
* 117-122: Setup directories (results, plots, intermediates). 
* 127-134: **Call** [get_tiles.r](#get_tiles.r), which receives as **input** sample name (`$A`), genome build (`$G`), genome size (`$S`) and output dir (`$O`). Warns that **it takes a while**.
* 136-143: **Call** [fix_tiles.r](#fix_tiles.r), which receives as **input** sample name (`$A`), genome build (`$G`), output from [get_tiles.r](#get_tiles.r) (translocation-blocks_50KB_merged.bedpe) and output dir (`$O`). 
* 145-155: **Call** [classify_tiles.r](#classify_tiles.r), which receives as **input** sample name (`$A`), genome build (`$G`), output from [fix_tiles.r](#fix_tiles.r) (translocation-blocks_50KB_merged_fixed.bedpe) and output dir (`$O`). 
* Final output is translocation-blocks_50KB_merged_fixed_classified.bedpe


 <a id="get_tiles.r"></a>
## get_tiles.r

**Summary:** Tiles at 50 kb with interactions above a manually set threshold in both the cis and trans space are selected. To optimize efforts, it is done in rounds: loops identified at the 2Mb resolution are subdivided into 200 kb tiles; and significant loops among those tiles are further subdivided into 50 kb tiles. Significant loops among 50 kb tiles are clusterized by proximity.

**Libraries:** 
* strawr: API for fast data extraction for .hic files that provides programmatic access to the matrices. It doesn't store the pointer data for all the matrices, only the one queried, and currently they are only supporting matrices (not vectors).
* ggplot2
* tidyr
* dplyr
* dbscan: A fast reimplementation of several density-based algorithms of the DBSCAN family. Includes the clustering algorithms DBSCAN (density-based spatial clustering of applications with noise) and HDBSCAN (hierarchical DBSCAN), the ordering algorithm OPTICS (ordering points to identify the clustering structure), shared nearest neighbor clustering, and the outlier detection algorithms LOF (local outlier factor) and GLOSH (global-local outlier score from hierarchies). The implementations use the kd-tree data structure (from library ANN) for faster k-nearest neighbor search. 
* parallel
* doParallel
* foreach: looping that supports parallel execution.

**Functions:** 
* `create_chr_bed`: Divide chromosome into chunks of a given `tile_size`, returns a bedfile format dataframe. 
* `create_tiled_bedpe_parallel`: Takes a data frame (`round1_tiles`) containing genomic interval pairs in bedpe-like format and generates smaller tile intervals from each pair. It is parallelized.

**Process:**
* 14-18: Obtain and open inputs (from  [get_tiles.r](#get_tiles.r), lines 127-134, `$A`, `$G`, `$S`, `$O`).
    - Sample name = `sample`
    - Genome build = `genome` 
    - Genoem size = `genome_size` &#8594; `chr_sizes`
    - Output directory = `out_dir`
* 125-452: Round 1: 2 Mb tiles
    Cis tiles
    - 141-147: **Call** `create_chr_bed` to generate a per-chromosome list of dataframes `chr_beds` where chromosomes are divided into `tile_size` (2 Mb) tiles. 
    - 149-170: Cis tiles (`cis_bedpe`): a bedpe-like dataframe where every tile in the chromosome is paired with the next adjacent tile.
    Trans tiles
    - 173-212: Interchromosomal trans tiles (`trans_bedpe`): Make all possible trans pairs (every chromosome with every other chromosome, and every coordinate in chr A with every coordinate in chr B). 
    - 213-230: Intrachromosomal trans tiles (`cis_chr_bedpe` &#8594;appended to&#8594;`trans_bedpe`): Make all possible cis pairs (every coordinate with every coordinate in the same chromosome) and keep only pairs that are at least 3.5 Mb apart. 
    - 231-249: Write `cis_bedpe` and `trans_bedpe` to files.
    - 251-273: **Call** `aqua_tools/annotate_loops.sh` on cis and trans bedpe files, with options `-Q cpm` (), `--formula sum` and `-R 100000`.
    > find out about aqua tools 
    - 275-296: Load results from  aqua tools into `cis_bedpe` and `trans_bedpe`.
    Thresholding
    - 301-303: Make log2 of non-zero cis count values (`counts_cis`). 
    - 305-323: Specify thresholds for inter and intra chromosomal trans tiles (manually depending on sample name; for unknown samples, the threshold is from 3 SDs below the mean to 1 SD below the mean).
    - 326-331: Select non-zero count values for interchromosomal trans tiles  (`trans_bedpe_inter`) and non-zero count values for intrachromosomal trans tiles (`trans_bedpe_intra`).
    - 333-437: Generate optional plots. The same number of points for cis and trans tiles from `cis_bedpe` and `trans_bedpe` (`_inter` and `_intra`) are selected (depending on shortest vector). After sampling, select non-zero values and transform into log2 space. Then plot.
    - 438-452: Select and save final trans tiles `round1_tiles`: all those inter and intra chromosomal trans tiles with log2(non-zero counts) above their respective thresholds, specified in l. 305-323. 
* 454-743: Round 2: 200 kb tiles. It is essentially the same as round 1, with some differences:
    - 468: Cis tiles
    - 503: Trans tiles: the final trans tiles in round 1 (`round1_tiles`) are recycled **calling** `create_tiled_bedpe_parallel`. When calling `aqua_tools/annotate_loops.sh` on cis and trans bedpe files, option R is different: `-R 5000`.
    - 588: Thresholding: manually-set thresholds are different. Final trans tiles saved to file are `round2_tiles`.
* 744-: Round 3: 50 kb tiles. It is essentially the same as round 1, with some differences:
    - 758: Cis tiles
    - 503: Trans tiles: the final trans tiles in round 2 (`round2_tiles`) are recycled **calling** `create_tiled_bedpe_parallel`. When calling `aqua_tools/annotate_loops.sh` on cis and trans bedpe files, option R is different: `-R 5000`.
    - 588: Thresholding: manually-set thresholds are different. Final trans tiles saved to file are `round2_tiles`.
    - 1033: Merging 50 kb tiles. Nearby tiles are clustered into larger blocks based on spatial proximity using DBSCAN (a density-based clustering algorithm). 

 <a id="fix_tiles.r"></a>
## fix_tiles.r

**Summary:** 
**Libraries:** 
* 
*
*

**Functions:** 
* _`get_counts`: Defined but not used_

**Process:**
* 54-84: Obtain and open inputs (from [call_loops.sh](#call_loops.sh), lines 332-355: `$file`, `$R`, path to copy of `$P` in the disk directory).
    - Scaffolds path = `path_scaffold` &#8594; `scaffold`. Chr prefix gets checked to match .hic file. 
    - Resolution = `bin_size`
    - .hic file = `path_hic`. Chr prefix gets checked.
* 103-130: Use an Alternative method instead of function `get_counts`. It iterates over each chromosome. **Call** straw, with **options** NONE normalization, path to .hic file, all the range of coordinates within this chromosome in unit = BP, at specified resolution. It generates a square matrix in `bin_resolution` bins with the counts from the .hic file. The diagonal counts are used to populate the empty `scaffold` bins.
* 138-140: Print the matrix. The result is saved as an annotated bed when this R function is called from [call_loops.sh](#call_loops.sh).

 <a id="classify_tiles.r"></a>
## classify_tiles.r

**Summary:** 

**Libraries:** 
* 
*
*
*

**Functions:** 
* `add_tad`: Takes x (active and inactive bins) and loop over TADs, assigning TAD IDs to those bins that overlap with any, and removing those bins that do not overlap with TADs.
* _`zero_diag`: Overwritten later_
* _`get_tad_matrix`: Defined but not used_
* _`get_tuple`: Defined but not used_
* _`zero_diag`: Defined but not used_

**Process:**

* 132-180: Obtain and open inputs (from [call_loops.sh](#call_loops.sh), lines 401-428:             `$tads`, `$output_dir`, `$N`, `$annotated_scaffold`, `$P`, `$R`, `$pivots_dir`).
    - TAD bedfile = `tad_path` &#8594; `tads` at line 328.
    - Output directory is set as working directory.
    - Output name = `sample_name` &#8594; `idx`.
    - Annotated scaffold bedfile = `scaffold_path` &#8594; `scaffold`. Column "unit" for unit type (diagonal, pro_off, enh_off, etc.) and counts column is named as `idx` (sample name). Later, chromosome names are checked and some constants regarding the unit types are defined. Specifically, "enh_off" and "pro_on" are stored as `inh_units`.
    - .hic file ($P) = `path_hic`.
    - Resolution = `bin_size`.
    - Tmp dir for storing pivots = `pivots_dir`.
    - Optional flags.

> To see what is happening, run the following code at line 213.
<pre>ggplot()+
    geom_point( aes(x =  log2(X[[unit]] + 2^(-10)), y = 0) )+ 
    geom_line(aes(x = dens$x, y = dens$y))+
    geom_line(aes(x = dens$x, y = y), color = "red", linetype = "dashed")+ 
    geom_vline(xintercept = pivot, color = "red")+
    ggtitle(unit)
</pre>

 


