#!/bin/bash

# Default values for optional arguments
ref=""                           # Mandatory reference file
ref_name="ref"                   # ref_name (default-"ref")
query=""                         # Mandatory query file and name pairs (mandatory)
query_name="query"               # query_name (default-"query")
output_dir="."                   # Optional output directory (default to current directory)
python_path="./SNPsyn_python.py" # Optional path to the python file (default- file exist in current directory)
kingdom="Bacteria"               # Optional kingdom (default "Bacteria")
genus="Escherichia"              # Optional genus (default "Escherichia")
threads=8                        # Optional - number of threads (default=8)



# Help message function
usage() {
    echo "Usage: $0 -r <reference_file> -q <query_file> [--ref_name <ref_name>] [--query_name <query_name>] [-o <output_directory>]  [--python_path <python_path>] [-k <kingdom>]  [-g <genus>]  [-h]"
    echo
    echo "Options:"
    echo "  -r, --ref            Reference file path (required)"
    echo "  --ref_name           Reference name, string (optional)"
    echo "  -q, --query          Query file and query name (required; can be repeated for multiple pairs)"
    echo "  --query_name         Query name, string (optional)"
    echo "  -o, --output_dir     Output directory for results (optional, defaults to current directory)"
    echo "  --python_path        path to the python file (optional,default- file exist in current directory)"
    echo "  -k, --kingdom        Optional kingdom (default Bacteria)"
    echo "  -g, --genus          Optional genus (default Escherichia)"
    echo "  -t, --threads            Optional - number of threads (default=8)"
    echo "  -h, --help           Show this help message"
    exit 1
}

SHORT=r:,q:,o:,k:,g:,t:,h
LONG=ref:,query:,ref_name:,query_name:,output_dir:,python_path:,kingdom:,genus:,threads:,help
OPTS=$(getopt -a -n SNPsyn --options $SHORT --longoptions $LONG -- "$@")


eval set -- "$OPTS"

while :
do
  case "$1" in
    -r | --ref )
      ref="$2"
      shift 2
      ;;
    -q | --query )
      query="$2"
      shift 2
      ;;
    --ref_name )
      ref_name="$2"
      shift 2
      ;;
    --query_name )
      query_name="$2"
      shift 2
      ;;
    -o | --output_dir )
      output_dir="$2"
      shift 2
      ;;
    --python_path )
      python_path="$2"
      shift 2
      ;;
    -k | --kingdom )
      kingdom="$2"
      shift 2
      ;;
    -g | --genus )
      genus="$2"
      shift 2
      ;;  
    -t | --threads )
      threads="$2"
      shift 2
      ;;
    -h | --help)
      usage
      exit 2
      ;;
    --)
      shift;
      break
      ;;
    *)
      echo "Unexpected option: $1"
      ;;
  esac
done

# Check if mandatory arguments are provided
if [ -z "$ref" ] || [ -z "$query" ]; then
    echo "Error: Both reference file (-r) and query file (-q) are required."
    usage
fi


echo "Reference: $ref"
echo "Query: $query"
echo "Reference Name: $ref_name"
echo "Query Name: $query_name"
echo "Output Directory: $output_dir"
echo "Python path: $python_path"
echo "kingdom: $kingdom"
echo "genus: $genus"


temp=$output_dir/temp
MUMMER_OUT=$temp/mummer_output
mummer_delta=$MUMMER_OUT/delta
mummer_results=$MUMMER_OUT/output
PROKKA_OUT=$temp/prokka_output
SNPsyn_output=$output_dir/$ref_name

mkdir $output_dir $temp $MUMMER_OUT $mummer_delta $mummer_results $PROKKA_OUT $SNPsyn_output

# MUMmer run 
echo "start MUMMer run"
cd $mummer_delta 
nucmer $ref $query --maxmatch --prefix=$ref_name"_"$query_name
show-snps -ClrT $mummer_delta/$ref_name"_"$query_name.delta > $mummer_results/$ref_name"_"$query_name.txt 

# Prokka run - customize the command for your genome
echo "start prokka run"
 prokka $ref --outdir $PROKKA_OUT --prefix $ref_name --locustag $ref_name --kingdom $kingdom --genus $genus --usegenus --force  --cpus $threads


# Call the Python script with arguments
# Running the Python script to analyze outputs
#three inputs: 1. mummer output of reference vs strain 2. prokka gff 3. fasta file of the refernce 
#then 2 ouputs - csv with snps and summary file

echo "calculates synonym vs nonsynonym snps"
python $python_path $mummer_results/$ref_name"_"$query_name.txt  $PROKKA_OUT/$ref_name.gff $ref $SNPsyn_output/$ref_name"_"$query_name"_SNPsyn".csv $SNPsyn_output/$ref_name"_"$query_name"_SNPsyn_summary".csv