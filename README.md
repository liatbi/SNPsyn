# SNPsyn
SNPsyn is a tool for analyzing SNPs (Single Nucleotide Polymorphisms) between a reference genome and query genomes. It provides annotations for each SNP, indicating whether it is synonymous, non-synonymous, an indel, or located in non-coding DNA (ncDNA).


**Prerequisites**
Bash: Ensure bash is installed as the script runs in a bash shell.
Python: Required to execute the Python 3.11.0 of SNPsyn (SNPsyn_python.py).
Required Python Libraries: pandas 1.5.1, Bio 1.6.2, sys
MUMmer3: MUMmer3.23, gcc-6.2.0
prokka: prokka-1.13.3, blast-2.6.0, perl/perl-5.26

**Installation**
Download the two necessary scripts to your working directory:
1. SNPsyn.sh
2. SNPsyn_python.py

**Usage**
Run the SNPsyn.sh script, specifying the reference and query files:

**Input**
Reference: A FASTA file containing the reference genome.
Query: A FASTA file containing the query genome(s).

**Output**
The output files are stored in the specified output directory and contain annotated SNP files with the following details:
1. POS- SNP Positions: Position information of each identified SNP compared to the reference genome
2. snp_pos_in_gene: The position of the SNP in the ORF- (first nucleotide in the orf will be in position "1")
3. SNP Type:
   * Synonym: SNPs that do not change the encoded amino acid.
   * Non-synonym: SNPs that result in an amino acid change.
   * Indel: Insertions or deletions.
   * ncDNA: SNPs located in non-coding DNA regions.
4. ORF SNPs: Number of SNPs identified in each open reading frame (ORF).
5. protein_covpident- represents the percent identity coverage: the percentage of identical nucleotides in the query ORF relative to the entire reference ORF.
6. SNPs in the ORF Context - the same annotation will be given to all the SNPs in the ORF
   * frame_shift - one or two indels in a row that made a frame shift
   * premature_stop_codon - the SNPs within the ORF lead to stop codons 
   * Amino acid changes (annotation in case there is no stop codon)
   * Synonymous mutations
   * Within non-coding DNA.
7. protein_changed - TRUE/FALSE indicating whether the protein has changed (not all SNPs are synonymous).
8. frame_shift - TRUE/FALSE indicating whether there is a frameshift in the ORF.
9. premature_stop_codon - TRUE/FALSE indicating whether a premature stop codon was generated due to SNPs in the ORF. 

**Examples**

Run SNPsyn with reference and query genome:
 
bash ./SNPsyn.sh SNPsyn -r /path/to/reference.fasta -q /path/to/query.fasta --ref_name ref_name --query_name query_name -o /path/to/output

Run SNPsyn with reference and query genome:

bash ./SNPsyn.sh SNPsyn -r /path/to/reference.fasta -q /path/to/query2.fasta --ref_name ref_name --query_name query_name --python_path ./SNPsyn_scripts/SNPsyn_python.py -o /path/to/output -k Bacteria -g Salmonella --threads 4


-r, --ref            Reference file path (required)   
--ref_name           Reference name, string (optional)  
-q, --query          Query file and query name (required; can be repeated for multiple pairs)  
--query_name         Query name, string (optional)  
-o, --output_dir     Output directory for results (optional, defaults to the current directory)  
--python_path        path to the python file (optional, default- file exists in the current directory)  
-k, --kingdom        Optional kingdom (default Bacteria)  
-g, --genus          Optional genus (default Escherichia)  
-t, --threads        Optional - number of threads (default=8)  
