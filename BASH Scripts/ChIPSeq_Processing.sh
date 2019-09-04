#!/bin/bash
#----------------------------------------------------------------------------------------------------------------------------------
# 					           ChIP-Seq Analysis
#----------------------------------------------------------------------------------------------------------------------------------
# This script takes as input FASTQ files from ChIP-Seq experiments and carries out all the required steps in order to
# obtain the desired output files in BedGraph and Wig formats. Its main argument is the folder containing the FASTQ files.
# The names of the processed files must be changed according to which information they store (Input or Histone H1) 
# because the 6th step consists on subtracting the Input signal (Control) from the histones H1 signal (Treatment).

#----------------------------------------------------------------------------------------------------------------------------------
# Cheking and storage of the Arguments
USAGE="USAGE= $0 -p path"
echo "$@"
if [ "$#" -ne "2" ]; then
	echo "Incorrect arguments: $USAGE"
	exit 1
fi
if [ "$1" == "-h" ] || [ "$1" == "--help" ]; then
	echo "$USAGE"
	exit 1
fi
case $1 in
	-p|--path)
	WorkingDir="$2"
esac

#----------------------------------------------------------------------------------------------------------------------------------
# Iterate over all the files included in the working directory with the .gz extension (compressed FASTQ files)
cd $WorkingDir
for file in *gz
	do echo "----------- PROCESSING SAMPLE: $file -------------"

# 1| Bowtie2 alignment
	echo "-------------------- Bowtie2 ---------------------"
	bowtie2 -x /home/ajvlab/Softwares/hg19/hg19 -p 10 -U $file -S "$(basename "$file").sam"
	rm $file

# 2| SAM to BAM Conversion
	echo "------------- SAM to BAM Conversion --------------"
	samtools view -b -S $(basename "$file").sam > "$(basename "$file").bam"

# 3| BAM Sorting
	echo "------------------ BAM Sorting -------------------"
	samtools sort $(basename "$file").bam > "sorted_$(basename "$file").bam"
	rm $(basename "$file").bam

# 4| BAM Filtering
	echo "----------------- BAM Filtering ------------------"
	samtools view -c sorted_$(basename "$file").bam
	samtools view -b -F 3844  sorted_$(basename "$file").bam > "filtered_sorted_$(basename "$file").bam"

# 5| BAM to BedGraph Conversion
	echo "---------- BAM to BedGraph Conversion ------------"
	N=$(samtools view -c filtered_sorted_$(basename "$file").bam)
	echo "$N"
	RPM=$(echo "1 / $N * 1000000" | bc -l)
	echo "$RPM"
	bedtools genomecov -ibam filtered_sorted_$(basename "$file").bam -bga -scale $RPM > "$(basename "$file").bdg"
done

#----------------------------------------------------------------------------------------------------------------------------------
# 6| Input Subtraction: Building Signal Tracks
for file in H1*bdg
	do echo "--------------- Input Subtraction: $file ----------------"
	macs2 bdgcmp -t $file -c Input*bdg -m subtract -o "$(basename "$file")_InputSubtracted.bdg"
	rm $file

# 7| BedGraph Sorting
	echo "------------------ BedGraph Sorting ---------------------"
	sort -k1,1 -k2,2n $(basename "$file")_InputSubtracted.bdg > "sorted_$(basename "$file")_InputSubtracted.bdg"
	rm $(basename "$file")_InputSubtracted.bdg

# 8| BedGraph Mapping against CytoBands
	for band in /home/ajvlab/Nuria/New_Cytobands/*bed
		do echo "-------- BedGraph Mapping: $band --------"
		bedtools map -a $band -b sorted_$(basename "$file")_InputSubtracted.bdg -c 4 -o mean > "$(basename "$file")_$(basename "$band")_overlap.bed"
	done

# 8| BedGraph to Wig Conversion
	echo "-------------- BedGraph to Wig Conversion ---------------"
	perl /home/ajvlab/Softwares/bedgraph_to_wig.pl --bedgraph sorted_$(basename "$file")_InputSubtracted.bdg --wig "$(basename "$file")_InputSubtracted.wig" --step 50
done

#----------------------------------------------------------------------------------------------------------------------------------


