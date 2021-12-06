# My normal script is very complicated to read so I have made this to keep it simple!
# This script won't run without actual input files, but it details how each step works and what and why each parament is used
# my examples use sample_p as +OHT smaple_m as -OHT





############################# Files #############################

# The mutseq_primers.csv file is just the amplicon name, the forward primer sequence and the reverse primer seqeunce as a csv



# the bed file of the coordinates is the chromosome, the start cordinate, end coordinate, the DSB position and the amplicon name for each amplicon with the primer sequences removed (explained more in the mutation calling section)
# but because it is my custom genome of the specific amplified regions with 50bp added either side, the chromosome name is the same as the amplicon name and the coordinates are all weird becasue the "chromosome" starts only 50bp upstream of the primer
# You can create your own genome or you can align to the whole genome and use a bed file of the genomic coordinates instead. Whole genome mapping considerably increases mapping time and mutlimapping, but should work overall.





############################# Translocation mapping #############################

# Translocapture is part of the mprofile package I made for the analysis, it can be isntalled via pip 'pip install mprofile-mut' or via Github https://github.com/aldob/mProfile 

# The t1/t2 parameters are translocated fastq outputs, the fq1/fq2 parameters are the non-translocated fastq outputs for further analysis

# the primer file should be a csv with: amplicon name, forward primer sequence, reverse primer sequence. The amplicon names will be used to create the column/row names in the table

TransloCapture -s 3 -1 sample_R1.fastq -2 sample_R2.fastq -o sample_tmap.csv -t1 sample_tmap_R1.fq -t2 sample_tmap_R2.fq -fq1 sample_norm_R1.fq -fq2 sample_norm_R2.fq -p mutseq_primers.csv

# the -pp parameter is for pre-processed, it then only inputs the translocation tables and will delta the -c input from the -i input and output a new table -o

TransloCapture -pp -i sample_p_tmap.csv -c sample_m_tmap.csv -o sample_delt_tmap.csv





############################# Qaulity filtering #############################

# I have found that stringent processing is best here to reduce the noise in the downstream results
# You lose a lot of reads, but the results end up cleaner as the lower quality reads have lots of errors 
# The -q and -u arguments may need altering if you have too few reads passing the filters

fastp -w 1 -A -l 30 -q 30 -u 20 -h sample_fastp.html -i sample_norm_R1.fq -I sample_norm_R2.fq --out1 sample_qc1.fq --out2 sample_qc2.fq





############################# Alignment #############################

# I have optimised the alignment parameters for Bowtie2 to get the highest alignment efficiency without getting many multimappers 
# I use a custom genome file based on hg38 that only has the amplified regions padded by 50bp either side, this improves alignment accuracy and greatly reduces processing time
# the custom genome was made by using getfasta to make a fasta file of regions specified by a bed file

# I pipe the bowtie2 results straight into samtools sort so the alignment is sorted while it is generated, but you can specify an output file and sort it after if you prefer 

bowtie2 -t -p 60 --fr --maxins 400 --no-discordant --no-mixed --ignore-quals --no-1mm-upfront -D 100 -R 50 -L 28 -N 1 --np 0 --dpad 49 --gbar 2 --mp 3.2,0.35 --rdg 1,1 --rfg 5,2 --score-min L,-1.0,-0.5 -x genome \
-1 sample_qc1.fq -2 sample_qc2.fq \
| samtools sort -o sample_sort.bam -O bam -l 0 -@ 60 -m 1G

samtools index sample_sort.bam





############################# Mutation calling #############################

# mpileup calculates the mutation for every read at every nucleotide
# the output file is very annoying to use which is why I wrote callMUT to convert it
# -aa ensures that even bases with 0 reads are still reported with a readcount of 0, otherwise some samples can have different numbers of nucleotides which will ruin the callMUT analysis
# --no-BAQ is important as normally BAQ filtering will remove reads with multiple errors thinking they are low quality, but in this case they are simply mutated
# -d limits the read depth, setting 0 should make it unlimited but this has caused problems for me in the past so I just stick to a really high number now
# The -l bedfile is the regions you want to analyse (one row for each amplicon: chromosome, first coordinate, second coordinate), remove the primers! You don't want to call mutations at the primers, so make the coordinates the first nucleotides inside each primer

samtools mpileup -aa --no-BAQ -d 1000000000 sample_sort.bam -o sample.mpileup -f genome.fa -l mutseq.bed



# callMUT will convert the mpileup to an mprofile which is much easier to use for downstream analysis
# -ic is a percentage of read cutoff for reporting individual indels. Sometimes if there are lots the files can get quite big, but most of the time just set NA unless you aren't interested in specific indel sequences

callMUT -i sample.mpileup -o sample.mprofile -ic NA

# again -pp makes callMUT accept pre-processed mprofiles and does a delta of the -c file from the -i file

callMUT -i sample_p.mprofile -c sample_m.mprofile -o sample_delt.mprofile -ic NA

# Then you have an mprofile which I generally used directly for plotting in R. 



# The imut_analysis.py python script is used to calculated the average mutation rate per DSB for different conditions
