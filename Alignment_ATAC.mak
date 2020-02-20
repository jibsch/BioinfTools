SHELL=/bin/bash
DATA_DIR=../data/
ALIGN_DIR=../align/
MACS_DIR=../macs/
THREADS=8
BWA_MM10=~/references/mouse/mm38_BWA_index/Mus_musculus.GRCm38.dna.primary_assembly
BWA_HG19=~/references/human/hg19_bwa_index
BL_MM10=~/references/mouse/mm10.blacklist.bed
BL_HG19=~/references/human/hg19-blacklist.v2.bed

ifdef path
TMP1 := ${ALIGN_DIR}$(path)/
ALIGN_DIR=$(TMP1)
TMP2 := ${DATA_DIR}$(path)/
DATA_DIR=$(TMP2)
endif


#Set index according to user paramter
ifeq ($(genome),hg19)
INDEX=$(BWA_HG19)
BL=$(BL_HG19)
MACS_ORG=hs
endif
ifeq ($(genome),mm)
INDEX=$(BWA_MM10)
BL=$(BL_MM10)
MACS_ORG=mm
endif

BAMS=$(subst R1_001.fastq.gz,q10_bl_filtered.bam, $(subst ${DATA_DIR},${ALIGN_DIR}, $(wildcard ${DATA_DIR}*R1_001.fastq.gz)))
PEAKS=$(subst _R1_001.fastq.gz,, $(subst ${DATA_DIR},${MACS_DIR}, $(wildcard ${DATA_DIR}*R1_001.fastq.gz)))
BIGWIGS=$(subst R1_001.fastq.gz,q10_bl_filtered.bam.bw, $(subst ${DATA_DIR},${ALIGN_DIR}, $(wildcard ${DATA_DIR}*R1_001.fastq.gz)))

#ls ../rna-seq/*R1* | sed -e 's/\.\.\/rna-seq\//${ALIGN_DIR}/' -e 's/_R1_001\.fastq\.gz/-starAligned.sortedByCoord.out.bam/'

all: ${BAMS} ${PEAKS}

bigwigs: ${BAMS} ${BIGWIGS}

#Ctrl-2_r1.fastq.gz  Ctrl-2_r2.fastq.gz dd
${ALIGN_DIR}%_s.bam: ${DATA_DIR}%_R1_001.fastq.gz ${DATA_DIR}%_R2_001.fastq.gz
	bwa mem -t ${THREADS} ${INDEX} $< $(word 2,$^) | samtools view -bS - > ${ALIGN_DIR}$*.bam
	samtools sort -@ 15 --output-fmt BAM -o $@ ${ALIGN_DIR}$*.bam
	samtools index $@


${ALIGN_DIR}%_q10_bl_filtered.bam: ${ALIGN_DIR}%_s.bam
	intersectBed -v -a $^  -b ${BL} -ubam | samtools view -h - | awk 'NR<69{print $$0} NR>=69{if($$5>=10) print $$0}' | samtools view -bS - > $@
	samtools index $@

${MACS_DIR}%: ${ALIGN_DIR}%_q10_bl_filtered.bam
	macs2 callpeak -t $^ -f BAMPE -g ${MACS_ORG} --keep-dup 1 --outdir $@  -n $* --nomodel

${MACS_DIR}%_joined_reps: ${ALIGN_DIR}%-1_q10_bl_filtered.bam ${ALIGN_DIR}%-2_q10_bl_filtered.bam
	macs2 callpeak -t $^ -f BAMPE -g mm --keep-dup 1 --outdir $@  -n $* --nomodel

${ALIGN_DIR}%.bam.bw: ${ALIGN_DIR}%.bam
	~/tools/bamToBigwig.sh $^ ~/references/human/Homo_sapiens.GRCh37.dna.primary_assembly.sizes

${MACS_DIR}merged_peaks.bed: ${PEAKS}
	mergePeaks
	awk 'NR>1{OFS="\t"; print $2, ".", "exon", $3, $4, ".", "+",".", "peak_id \""$1"\"" }' mergedPeaks.bed > merged_peaks.gtf


closestBed -t first -D b -a <(awk 'NR>1{OFS="\t"; print $2, $3, $4, $1}' mergedPeaks.bed | sort -k1,1 -k2n,2) -b ~/references/human/Homo_sapiens.GRCh37.87_TSS.bed | awk '{if($NF>-5000 && $NF<5000) print $0}' > mergedPeaks_toGene5k.bed

	awk 'BEGIN{o=""} {if($8!=o) {c=1; o=$8}; OFS="\t"; print $1, ".", "exon", $2, $3, ".", "+",  ".", "gene_id "$8"-"c; c+=1}' mergedPeaks_toGene5k.bed > mergedPeaks_toGene5k.gtf
