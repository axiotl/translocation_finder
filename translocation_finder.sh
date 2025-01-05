#!/bin/bash 

data_dir=~/lab-data
script_dir=~/aqua_scripts

sample_sheet="/home/ubuntu/shared-files/setup/sample_sheet.txt"

if [ ! -f "$sample_sheet" ]; then
    python3 /home/ubuntu/aqua_tools/restore_sample_sheet.py > /dev/null 2>&1 

    if [ ! -f "$sample_sheet" ]; then
        echo "Failed to restore sample sheet. Exiting."
        exit 1
    fi
fi

function usage {
    echo -e "usage: "
    echo -e "  translocation_finder.sh \\"
    echo -e "    -A SAMPLE_NAME \\"
    echo -e "    -G GENOME_BUILD \\"
    echo -e "    -g GENOME_SIZE \\"
    echo -e "    -O OUT_DIR \\"
    echo -e "   [-h]"
    echo -e "Use option -h|--help for more information"
}


function help {
    echo
    echo "Outputs bedpe files containing translocations using H3K27ac HiChIP .hic files"
    echo
    echo "---------------"
    echo "OPTIONS"
    echo
    echo "    -A|--sample   : Name of sample as it appears on Tinkebox"
    echo "    -G|--genome   : Genome build"
    echo "    -S|--size     : Full path to 2-col .txt file containing chr sizes"
    echo "    -O|--out      : Full path to where the outputs should be made"
    echo "  [ -h|--help   ] : Help message."
    exit;
}

if [ $# -lt 1 ]
    then
    usage
    exit
fi

# Transform long options to short ones
for arg in "$@"; do
  shift
  case "$arg" in
      "--sample")   set -- "$@" "-A" ;;
      "--genome")   set -- "$@" "-G" ;;
      "--size")     set -- "$@" "-S" ;;
      "--out")      set -- "$@" "-O" ;;
      "--help")     set -- "$@" "-h" ;;
       *)           set -- "$@" "$arg"
  esac
done



while getopts ":A:G:S:O:h" OPT
do
    case $OPT in
  A) A=$OPTARG;;
  G) G=$OPTARG;;
  S) S=$OPTARG;;
  O) O=$OPTARG;;
  h) help ;;
  \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      exit 1
      ;;
  :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      exit 1
      ;;
    esac
done


# Do all necessary parameter checks
#----------------------------------

if [[ -z $A ]];
then
    usage
    exit
fi

if [[ -z $G ]];
then
    usage
    exit
fi

if [[ -z $S ]];
then
    usage
    exit
fi


if [[ -z $O ]];
then
    usage
    exit
fi

#----------------------------------


mkdir -p "$O"
mkdir -p "$O/results"
mkdir -p "$O/plots"
mkdir -p "$O/intermediates"


echo "Translocation Finder: $A"
echo

# tile trans space and get 50kb tiles that contain potential translocations of interest
echo "Fetching tiles across the entire cis and trans space... (this takes a while)"
Rscript $script_dir/get_tiles.r \
 $A \
 $G \
 $S \
 $O
echo

# expand obtained tiles to accomodate missing diagonal pixels
echo "Fixing borders of obtained tiles"
Rscript $script_dir/fix_tiles.r \
 $A \
 $G \
 "${O}/results/${A}_translocation-blocks_50KB_merged.bedpe" \
 $O
echo

# classify tiles into: a) gradient b) uniform c) floating diagonals
echo "Classifying tiles into: a) gradient b) uniform c) floating diagonals"
Rscript $script_dir/classify_tiles.r \
 $A \
 $G \
 "${O}/results/${A}_translocation-blocks_50KB_merged_fixed.bedpe" \
 $O 
echo 
echo
echo "Done! Final .bedpe is ${O}/results/${A}_translocation-blocks_50KB_merged_fixed_classified.bedpe"




