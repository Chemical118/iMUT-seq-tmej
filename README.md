# iMUT-seq
Analytical scripts for the iMUT-seq manuscript that is currently in preparation.

<br>

The shell script provides a full breakdown of the analytical pipeline starting with raw fastq files going through to fully process mprofile mutation mapping files. 
This script is not designed to be run as is, but instead provides a simplified list of commands with the exact parameters necessary for iMUT-seq analysis.

<br>

The mutseq.bed and mutseq_primers.csv files are required for this shell script pipeline, their functions are explained within the script.

<br>

The python script takes mprofile files and calculates MMEJ rates and base substitution signatures, then outputs tables of average mutation rates per loci and per amplicon type as well as tables of total deletion lengths, MMEJ deletion lengths and MMEJ homology lengths.

<br>
The R script contains the ggplot commands used to create the mutation line plots presented in the iMUT-seq manuscript.
