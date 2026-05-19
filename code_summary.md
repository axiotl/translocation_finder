# translocation_finder

 **Table of contents:**

 - [translocation_finder.sh](#translocation_finder.sh)
 - [get_tiles.r](#get_tiles.r)
 - [fix_tiles.r](#fix_tiles.r)
 - [classify_tiles.r](#classify_tiles.r)


 <a id="translocation_finder.sh"></a>
## translocation_finder.sh

**Summary:** Outputs bedpe files containing chromothripsis events classified as diagonals, gradients and uniforms from H3K27ac HiChIP data.

**Process:**
* 1-115: Arguments and help
    - `$A` = Sample name
    - `$G` = Genome build
    - `$S` = Chromosome sizes (2 columns: chromosome and size)
    - `$O` = Out path 

* 117-122: Setup directories (results, plots, intermediates). 
* 127-134: **Call** [get_tiles.r](#get_tiles.r), which receives as **input** sample name (`$A`), genome build (`$G`), genome size (`$S`) and output dir (`$O`). Warns that **it takes a while**.
* 136-143: **Call** [fix_tiles.r](#fix_tiles.r), which receives as **input** sample name (`$A`), genome build (`$G`), output from [get_tiles.r](#get_tiles.r) (translocation-blocks_50KB_merged.bedpe) and output dir (`$O`). 
* 145-155: **Call** [classify_tiles.r](#classify_tiles.r), which receives as **input** sample name (`$A`), genome build (`$G`), output from [fix_tiles.r](#fix_tiles.r) (translocation-blocks_50KB_merged_fixed.bedpe) and output dir (`$O`). 
* Final output is translocation-blocks_50KB_merged_fixed_classified.bedpe


<a id="get_tiles.r"></a>
## get_tiles.r

**Summary:** Select tiles at 50 kb with interactions with CPM values above a threshold based on diagonal values. To optimize efforts, it is done in rounds: loops identified at the 2Mb resolution are subdivided into 200 kb tiles; and significant loops among those tiles are further subdivided into 50 kb tiles. Significant loops among 50 kb tiles are clusterized by proximity.

**Functions:**

- `create_chr_bed`: Divide chromosome into chunks of a given `tile_size`, returns a bedfile format dataframe.
- `create_tiled_bedpe_parallel`: Takes a data frame (`round1_tiles`) containing genomic interval pairs in bedpe-like format and generates smaller tile intervals from each pair. It is parallelized.

**Process:**

- 14-18: Obtain and open inputs (from [translocation_finder.sh](#translocation_finder.sh), lines 127-134, `$A`, `$G`, `$S`, `$O`).
- Sample name = `sample`
- Genome build = `genome`
- Genome size = `genome_size` → `chr_sizes`
- Output directory = `out_dir`
- **125-452: Round 1: 2 Mb tiles**
    - **Cis tiles**
    - 141-147: **Call** `create_chr_bed` to generate a per-chromosome list of dataframes `chr_beds` where chromosomes are divided into `tile_size` (2 Mb) tiles.
    - 149-170: Cis tiles (`cis_bedpe`): a bedpe-like dataframe where every tile in the chromosome is paired with the next adjacent tile.
    - **Trans tiles**
    - 173-212: Interchromosomal trans tiles (`trans_bedpe`): Make all possible trans pairs (every chromosome with every other chromosome, and every coordinate in chr A with every coordinate in chr B).
    - 213-230: Intrachromosomal trans tiles (`cis_chr_bedpe` →appended to→`trans_bedpe`): Make all possible cis pairs (every coordinate with every coordinate in the same chromosome) and keep only pairs that are at least 3.5 Mb apart.
    - 231-249: Write `cis_bedpe` and `trans_bedpe` to files.
    - 251-273: **Call** `aqua_tools/annotate_loops.sh` on cis and trans bedpe files, with options `-Q cpm`, `--formula sum` and `-r 100000`.
        - 275-296: Load results from aqua tools into `cis_bedpe` and `trans_bedpe`.
    - **Thresholding**
    - 301-303: Make log2 of non-zero cis count values (`counts_cis`).
    - 305-323: Specify thresholds for inter and intra chromosomal trans tiles (manually depending on sample name; for unknown samples, the threshold is from 3 SDs below the mean to 1 SD below the mean).
    - 326-331: Select non-zero count values for interchromosomal trans tiles (`trans_bedpe_inter`) and non-zero count values for intrachromosomal trans tiles (`trans_bedpe_intra`).
    - 333-437: Generate optional plots. The same number of points for cis and trans tiles from `cis_bedpe` and `trans_bedpe` (`_inter` and `_intra`) are selected (depending on shortest vector). After sampling, select non-zero values and transform into log2 space. Then plot.
    - 438-452: Select and save final trans tiles `round1_tiles`: all those inter and intra chromosomal trans tiles with log2(non-zero counts) above their respective thresholds, specified in l. 305-323.
- **454-743: Round 2: 200 kb tiles.** It is essentially the same as round 1, with some differences:
    - **Cis tiles**
    - **Trans tiles**: the final trans tiles in round 1 (`round1_tiles`) are recycled **calling** `create_tiled_bedpe_parallel`. When calling `aqua_tools/annotate_loops.sh` on cis and trans bedpe files, option r is different: `-r 5000`.
    - **Thresholding**: manually-set thresholds are different. Final trans tiles saved to file are `round2_tiles`.
- **744-: Round 3: 50 kb tiles.** It is essentially the same as round 1, with some differences:
    - **Cis tiles**
    - **Trans tiles**: the final trans tiles in round 2 (`round2_tiles`) are recycled **calling** `create_tiled_bedpe_parallel`. When calling `aqua_tools/annotate_loops.sh` on cis and trans bedpe files, option r is different: `-r 5000`.
    - **Thresholding**: manually-set thresholds are different. Final trans tiles saved to file are `round3_tiles`.
- 1033: **Merging** 50 kb tiles. Nearby tiles are clustered into larger blocks based on spatial proximity using DBSCAN (a density-based clustering algorithm).

<a id="fix_tiles.r"></a>
## fix_tiles.r

**Summary:** Each loop (tile) from [get_tiles.r](#get_tiles.r) is refined by exploring the nearby coordinates at 100 kb around (x and y axes) and taking the maximum counts at each direction within the area. This expanded region is then trimmed to include only rows and columns with count data, determining the final coordinates. 

**Functions:**

- `process_tile` : *Macro function called in chunks and parallelized within a loop*. *Explained in Process.*

**Process:**

- 12-54: Obtain and open inputs (from [translocation_finder.sh](#translocation_finder.sh), lines 136-143, `$A`, `$G`, `$O` and 1 more), and additional setup actions.
    - Sample name = `sample`
    - Genome build = `genome`
    - Output from [get_tiles.r](#get_tiles.r) (translocation-blocks_50KB_merged.bedpe) = `bedpe`
    - Output directory = `out_dir`
- 294-346: Process `con_in` (connection to `bedpe` file path) in chunks of `chunk_size` tiles (set to 1000), calling `process_tile` and checking multiple times for errors. Valid results are stored in `out_file`.
- **Calls to process_tile:** receive as input the `tile`, `path_hic` which is the local, automatically-built path to the .hic file from the sample name, `search_space`, set to 100,000.
    - 59-63: Sort coordinates to have coord.1 before coord.2 within the tile (coord pair).
    - 66-74: Extract Hi-C contact counts with straw in 5kb bins to `tile_contacts`.
    - 76-150: **Expand coordinates.** Obtain Hi-C counts with straw for the tile plus 100kb to the left (`search_space`), removing those contacts already in `tile_contacts` → `left_contacts`. If there is a bin with higher counts than the maximum counts in `tile_contacts`, adjust the left coordinate to the bin with highest counts in `left_contacts`, plus one bin. Do the same to the right, up and down.
    - 228-291: **Fix expanded tile borders.** Define the new tile with expanded coordinates. Obtain Hi-C counts for that tile with straw, at 5 kb bins. Generate an empty matrix representing the bins within the tile. Fill the matrix with the counts, handling cases where the x and y coordinates may be swapped due to coordinate symmetry. Remove rows and columns without counts, and use the remaining matrix to define final tile coordinates.

<a id="classify_tiles.r"></a>
## classify_tiles.r

**Summary:** Classifies each tile from [fix_tiles.r](#fix_tiles.r) into one of three categories: gradient, uniform, or floating diagonal, based on Hi-C contact patterns within the tile.

**Functions:**

- `calculate_clumpiness`: Quantify how "clumpy" or spatially clustered is the signal (non-zero values) within a Hi-C counts matrix. While examining the matrix in sliding windows of size *n x n*:
    1. Calculate the proportion of non-zero values out of *n²* tiles (density).
    2. Find sliding windows with the maximum and minimum observed density.
    3. Take the mean of the 3 largest values in the maximum density window and the mean of the 3 smallest values in the **minimum** density window.
    4. Returns a list of three metrics:
        - **`max_density`**: The maximum density value found among *n x n* regions.
        - **`density_variance`**: Variance in density values (higher = more heterogeneous).
        - **`separation_score`**: Contrast between max and min intensities; the difference between the two quantities in step 3.
- `sort_bedpe`: Sort BEDPE dataframe, with chr1 < chr2 for all trans rows, and start1 < start2 for all cis rows.

**Process:**

- 12-18: Obtain and open inputs (from [translocation_finder.sh](#translocation_finder.sh), lines 145-155, `$A`, `$G`, `$O` and 1 more), and additional setup actions.
    - Sample name = `sample`
    - Genome build = `genome`
    - Output from [fix_tiles.r](#fix_tiles.r) (translocation-blocks_50KB_merged_fixed.bedpe) = `bedpe`
    - Output directory = `out_dir`
- 48-58: Generate flags to indicate whether there are cis and/or trans relationships within the input `whole_bedpe`.
- 173: Get `whole_bedpe` sorted **calling** `sort_bedpe`.
- 176-471: **CIS translocations** (if flag present).
    - 186: Select CIS relationships (chr1 == chr2).
    - **Gradient tiles**
        - 193-198: Generate empty matrix `feature_matrix` to store max density, density variance and separation score.
        - 201-246: For each CIS tile, obtain the contact counts with `straw` and store them into a matrix. Then, **call** `calculate_clumpiness` to fill `feature_matrix`.
        - *252-297: Deprecated method.*
        - 298-307: Set a cutoff based on the separation scores obtained for each CIS tile using `calculate_clumpiness`: 2 ^ (mean + 1SD in log2 space).
        - 309-328: From the bedpe with all CIS relationships, select those with ≥ 0.7 maximum density; and those with density between 0.5 and 0.7 with a separation above the cutoff value. These relationships form the `gradient_bedpe`; and all the other relationships represent the `uniform_bedpe`.
        - 329-344: Print plots.
    - **Floating diagonals.** For each tile in `gradient_bedpe`:
        - 350-404: Obtain Hi-C counts using `straw`.
        - 405-431: Find the maximum Hi-C contact value and obtain diagonal vectors from that position using `sapply`. Find out the `floating_score`: proportion of counts in the diagonal from the total counts in the tile; how diagonally enriched it is.
        - 434-469: Classify tiles into two groups according to the floating score using `kmeans`. The group with high floating score is classified as `floating_diagonal` and the other group as `gradient` (there's also an ambiguous option).
- 472-764: **TRANS translocations** (if flag present).
    - 482: Select TRANS relationships (chr1 != chr2).
    - **Gradient tiles — exactly the same as cis relationships.**
        - 489-494: Generate empty matrix `feature_matrix` to store max density, density variance and separation score.
        - 497-542: For each TRANS tile, obtain the contact counts with `straw` and store them into a matrix. Then, **call** `calculate_clumpiness` to fill `feature_matrix`.
        - *548-594: Deprecated method.*
        - 594-603: Set a cutoff based on the separation scores obtained for each TRANS tile using `calculate_clumpiness`: 2 ^ (mean + 1SD in log2 space).
        - 605-622: From the bedpe with all TRANS relationships, select those with ≥ 0.7 maximum density; and those with density between 0.5 and 0.7 with a separation above the cutoff value. These relationships form the `gradient_bedpe`; and all the other relationships represent the `uniform_bedpe`.
        - 623-639: Print plots.
    - **Floating diagonals — exactly the same as cis relationships.** For each tile in `gradient_bedpe`:
        - 645-698: Obtain Hi-C counts using `straw`.
        - 700-724: Find the maximum Hi-C contact value and obtain diagonal vectors from that position using `sapply`. Find out the `floating_score`: proportion of counts in the diagonal from the total counts in the tile; how diagonally enriched it is.
        - 729-762: Classify tiles into two groups according to the floating score using `kmeans`. The group with high floating score is classified as `floating_diagonal` and the other group as `gradient` (there's also an ambiguous option).
- 769-784: Join dataframes for classified CIS and TRANS tiles and write to file.