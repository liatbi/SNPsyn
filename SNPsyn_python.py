import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import sys

mummer_snps_path = sys.argv[1]
prokka_path = sys.argv[2]
fasta_path = sys.argv[3]
output_file = sys.argv[4]  # Output filename
summary_file = sys.argv[5] #summary file

mummer_snps=pd.read_csv(mummer_snps_path, header=2, sep="\t")
orf_df=pd.read_csv(prokka_path, sep='\t', header=2)

#df for ref vs strain - mummer output into df
#order the column names correctly
cols= mummer_snps.columns[:].tolist()
cols=[item.replace('[','').replace(']','') for item in cols] #drop []
cols.extend(['ref_contig','strain_contig']) #add columns name contig
mummer_snps.reset_index(inplace=True)
mummer_snps.columns=cols
mummer_snps.rename(columns={'P1': 'POS', 'SUB': 'REF', 'SUB.1': 'ALT'}, inplace=True)
#find main contig - and make a df only with the main (without plasmids)
main_contig=mummer_snps['ref_contig'].value_counts().idxmax() #find the main contig of the strain ( in case of  multiple contigs)
mummer_snps=mummer_snps[mummer_snps['ref_contig'] == main_contig] #use only main contig (drop snps from plasmids)
print("number of snps positions (unique per position): ", len(mummer_snps['POS'].unique()))
print("number of snps: ", len(mummer_snps['POS']))
mummer_snps


#prokka output to get orfs
#shift the column row to first row and set column names
orf_df=orf_df.shift(1)
orf_df.iloc[0]=orf_df.columns
orf_df.columns=['ref_contig','source','feature_type','start','end','score','strand','phase', 'attributes']
orf_df.dropna(subset='feature_type', inplace=True)
#convert start and stop to int
orf_df['start'] = orf_df['start'].astype(int)
orf_df['end'] = orf_df['end'].astype(int)
#split the attributes col
#attribute_split=orf_df['attributes'].str.split(';', expand=True)
orf_df['ID']=orf_df['attributes'].str.split('ID=').str[1].str.split(';').str[0].str.split('_').str[-1]
orf_df['Name']=orf_df['attributes'].str.split('Name=').str[1].str.split(';').str[0]
orf_df['gene']=orf_df['attributes'].str.split('gene=').str[1].str.split(';').str[0]
orf_df['product']=orf_df['attributes'].str.split('product=').str[1].str.split(';').str[0]
orf_df

#add the ID from the prokka df that suitable for the position of the snip
for pos in (mummer_snps['POS'].unique()):
    ID_df=orf_df[(orf_df['start']<pos) & (orf_df['end']>pos)]['attributes']
    if ID_df.empty: #snp not in orf
            pass
    elif len(ID_df)==1: #one orf
        ID=ID_df.item()
        mummer_snps.loc[mummer_snps['POS']==pos,'attributes']=ID
    else: #more than one orf
        orf_count=len(ID_df) #numbers of orfs
        IDs=ID_df.tolist()
        rows_to_duplicate=mummer_snps[mummer_snps['POS'] == pos]
        #add first ID
        ID=IDs[0]
        print(IDs)
        mummer_snps.loc[mummer_snps['POS']==pos,'attributes']=ID
        for ID in IDs[1:]:
            new_rows=rows_to_duplicate.copy()
            new_rows['attributes']=ID
            mummer_snps=pd.concat([mummer_snps,new_rows], ignore_index=True)
            print(pos)
mummer_snps.sort_values(by=['POS','attributes'], inplace=True, ascending=True)
mummer_snps


#merge with orf df to get merged snps with its gene orf and info
snps_orf_df=mummer_snps.merge(orf_df, on=['attributes','ref_contig'], how='left')
snps_orf_df['feature_type'].fillna('ncDNA', inplace=True)
snps_orf_df


#add number of snps in orf
num_of_snps_in_orf=snps_orf_df.groupby(['ID'])['POS'].count().rename('num_of_snps_in_orf')
snps_orf_df=snps_orf_df.merge(num_of_snps_in_orf, on='ID', how='left')
snps_orf_df


#add syn or non syn for orf for single snp
ref_sequence=SeqIO.read(fasta_path, 'fasta')
syn=0; nonsyn=0; wrong=0; indel=0
#fasta.seq(start:end)
for index,row in snps_orf_df[snps_orf_df['feature_type']!='ncDNA'].iterrows(): #only coding
    ref=row['REF']
    alt=row['ALT']
    if (ref!='.') & (alt!='.'):
        start=int(row['start'])
        end=int(row['end'])
        snp_pos=row['POS']
        alt=row['ALT']
      #  print(start, end, snp_pos, row['REF'], alt, row['strand'])
        #generate sequence
        alt_ref_sequence = ref_sequence[:snp_pos-1] + alt + ref_sequence[snp_pos:] # whole ref genome with one snp
        if row['strand']=='+': #pos strand
            dna_sequence= ref_sequence.seq[start - 1:end]
            alt_dna_sequence= alt_ref_sequence.seq[start - 1:end]
        elif row['strand']=='-': #neg strand
            dna_sequence= ref_sequence.seq[start - 1:end].reverse_complement()
            alt_dna_sequence= alt_ref_sequence.seq[start - 1:end]
        #give notification if it not start ot stop in start or stop codon
        if (dna_sequence[0:3] != 'ATG') & (dna_sequence[0:3] != 'GTG')& (dna_sequence[0:3] != 'TTG'):
                print('not start codon')
                print(dna_sequence)
        if (dna_sequence[-3:] != 'TAA') & (dna_sequence[-3:] != 'TAG') & (dna_sequence[-3:] != 'TGA'):
                print('not stop codon')
                print(dna_sequence)
        rna_seq = dna_sequence.transcribe()
        protein_seq = rna_seq.translate()
        alt_rna_seq = alt_dna_sequence.transcribe()
        alt_protein_seq = alt_rna_seq.translate()
        if (dna_sequence!=alt_dna_sequence): #dna should not be the same this should always be true
             if (protein_seq==alt_protein_seq):
                syn+=1
                snps_orf_df.at[index, 'syn_nonsyn'] = 'syn'
             else:
                nonsyn+=1
                snps_orf_df.at[index, 'syn_nonsyn'] = 'nonsyn'
            #    snps_orf_df.at[index, 'nonsyn_type'] = 'AA_change'
        else:
            wrong+=1
    else:
        indel+=1
        snps_orf_df.at[index, 'syn_nonsyn'] = 'indel'
        #snps_orf_df.at[index, 'syn_nonsyn'] = 'nonsyn
        #snps_orf_df.at[index, 'nonsyn_type'] = 'indel'
dna_sequence
print(syn,nonsyn,indel, wrong)
snps_orf_df


#fill in the syn_nonsyn column if it non coding dna etc...
snps_orf_df['syn_nonsyn'].fillna(snps_orf_df['feature_type'], inplace=True)
snps_orf_df


#add syn or non syn for whole orf
ref_sequence=SeqIO.read(fasta_path, 'fasta')
snps_orf_df['orf_syn_non_syn'] = None
#snps_orf_df['orf_syn_non_syn_type'] = None
for orf in snps_orf_df['ID'].dropna().unique(): #go orf by orf
    alt_ref_sequence=ref_sequence #initializing the sequence
    start=int(snps_orf_df[snps_orf_df['ID']==orf]['start'].iloc[0]) #strat and and should be similar to all snps in the orf
    end=int(snps_orf_df[snps_orf_df['ID']==orf]['end'].iloc[0]) #strat and and should be similar to all snps in the orf
    addtional=0
    # print(start, end, snp_pos, row['REF'], alt, row['strand'])
    for index,row in snps_orf_df[snps_orf_df['ID']==orf].iterrows(): # all the rows in this orf
        ref=row['REF']
        alt=row['ALT']
        snp_pos=row['POS']+addtional
        alt=row['ALT']
       # print(start, end, snp_pos, row['REF'], alt, row['strand'])
        #generate sequence
        if (ref!='.') & (alt!='.'):
            alt_ref_sequence = alt_ref_sequence[:snp_pos-1] + alt + alt_ref_sequence[snp_pos:] # adding snp
        elif ref=='.':
            alt_ref_sequence=alt_ref_sequence[:snp_pos-1] + alt + alt_ref_sequence[snp_pos-1:]
            #if row['strand']=='+': #pos strand
            addtional+=1
          #  print('ref=.', snp_pos, alt, start, end)
            #elif row['strand']=='-': #neg strand
        elif alt=='.':
            alt_ref_sequence=alt_ref_sequence[:snp_pos-1] + alt_ref_sequence[snp_pos:]
            addtional-=1
           # print('alt=.')
    #sequence contains all snps of the orf
    #   print(orf)
    if row['strand']=='+': #pos strand
        dna_sequence= ref_sequence.seq[start - 1:end]
        alt_dna_sequence= alt_ref_sequence.seq[start - 1:end+addtional]
    elif row['strand']=='-': #neg strand
        dna_sequence= ref_sequence.seq[start - 1:end].reverse_complement()
        alt_dna_sequence= alt_ref_sequence.seq[start - 1:end+addtional]
    #give notification if it not start ot stop in start or stop codon
    if (dna_sequence[0:3] != 'ATG') & (dna_sequence[0:3] != 'GTG')& (dna_sequence[0:3] != 'TTG'):
        print('not start codon')
        print(dna_sequence)
    if (dna_sequence[-3:] != 'TAA') & (dna_sequence[-3:] != 'TAG') & (dna_sequence[-3:] != 'TGA'):
        print('not stop codon')
        print(dna_sequence)
    rna_seq = dna_sequence.transcribe()
    protein_seq = rna_seq.translate()
    alt_rna_seq = alt_dna_sequence.transcribe()
    alt_protein_seq = alt_rna_seq.translate()
    #adding meaning of syn non syn in case of all the orf contect
    if (dna_sequence!=alt_dna_sequence): #dna should not be the same this should always be true
        if (protein_seq==alt_protein_seq): #whole protein is synonym
            syn+=1
            snps_orf_df.loc[snps_orf_df['ID'] == orf, 'orf_syn_non_syn'] = 'syn'
        else: #non synonym
            nonsyn+=1
          #  snps_orf_df.loc[snps_orf_df['ID'] == orf, 'orf_syn_non_syn'] = 'nonsyn'
            stop_codon_position = alt_protein_seq.find('*')
            if (stop_codon_position>-1) & (stop_codon_position+1!=len(protein_seq)):
             #   snps_orf_df.loc[snps_orf_df['ID'] == orf, 'orf_syn_non_syn_type'] = 'orf_contains_snp_of_stop_codon'
                snps_orf_df.loc[snps_orf_df['ID'] == orf, 'orf_syn_non_syn'] = 'orf_contains_snp_of_stop_codon'
            else:
              #  snps_orf_df.loc[snps_orf_df['ID'] == orf, 'orf_syn_non_syn_type'] = 'AA_changes'
                snps_orf_df.loc[snps_orf_df['ID'] == orf, 'orf_syn_non_syn'] = 'AA_changes'

#fill in the syn_nonsyn column if it non coding dna etc...
snps_orf_df['orf_syn_non_syn'].fillna(snps_orf_df['feature_type'], inplace=True)
snps_orf_df

SNPsyn_output = snps_orf_df[['ID','POS','REF','ALT','strand','feature_type', 'num_of_snps_in_orf','syn_nonsyn','orf_syn_non_syn', 'start', 'end', 'product','Name','gene','ref_contig', 'strain_contig','P2','attributes', 'source','phase','BUFF', 'DIST', 'LEN R', 'LEN Q', 'FRM', 'TAGS' ,'score']] #'nonsyn_type'
SNPsyn_output


summary_SNPsyn=snps_orf_df['syn_nonsyn'].value_counts()
summary_SNPsyn

#export
SNPsyn_output.to_csv(output_file, index=False)
summary_SNPsyn.to_csv(summary_file)