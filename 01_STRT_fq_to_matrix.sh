#!/usr/bin/bash

# ===============================================================
# File: STRT_fq_to_matrix.sh
# File Created: Th Sep 2023
# Author: Zhang Bin 
# -----
# Last Modified: Thu Sep 07 2023
# Modified By: Zhang Bin 
# ===============================================================

#
fastp=
samtools=
#
STAR=
index=
#
featureCounts=
gtf=
#
perl=
STRT_fq_format=./steps/STRT_fq_format.pl
STRT_umi_count=./steps/STRT_umi_count.pl
STRT_cell_mapping_ratio=./steps/STRT_cell_mapping_ratio.pl
#
Rscript=
matrix_to_seurat_rds=./steps/matrix_to_seurat_rds.R

#######################################################################################

#
barcode=
data_dir=
out_dir=

#
thread=20
sample_list=

for sample in ${sample_list[@]}; do

    sample_dir=${out_dir}/${sample}
    [[ -e ${sample_dir} ]] || mkdir ${sample_dir}

    extracted_dir=${sample_dir}/fastq_extracted
    [[ -e ${extracted_dir} ]] || mkdir ${extracted_dir}

    qc_dir=${sample_dir}/qc_fastp
    [[ -e ${qc_dir} ]] || mkdir ${qc_dir}

    mapping_dir=${sample_dir}/mapping_star
    [[ -e ${mapping_dir} ]] || mkdir ${mapping_dir}

    count_dir=${sample_dir}/count_featurecounts
    [[ -e ${count_dir} ]] || mkdir ${count_dir}

    matrix_dir=${sample_dir}/count_matrix
    [[ -e ${matrix_dir} ]] || mkdir ${matrix_dir}

    pbs=${sample_dir}/${sample}.pbs

cat >${pbs} <<-EOF
#!/bin/bash

#PBS -N ${sample}
#PBS -o ${sample_dir}/${sample}.out
#PBS -e ${sample_dir}/${sample}.err
#PBS -d ${sample_dir}
#PBS -l nodes=1:ppn=${thread}
#PBS -q long
#PBS -V

{

    echo "01 format fastq \`date\` !"

    time ${perl} ${STRT_fq_format} \\
        -w ${thread} \\
        -r1 ${data_dir}/${sample}_R1.fq.gz \\
        -r2 ${data_dir}/${sample}_R2.fq.gz \\
        -b ${barcode} \\
        -o ${extracted_dir}/${sample}_extracted_R1.fq.gz \\
    2>&1 |
    tee ${extracted_dir}/${sample}.log

    echo "02 QC \`date\` !"

    time ${fastp} \\
        -w ${thread} \\
        -i ${extracted_dir}/${sample}_extracted_R1.fq.gz \\
        -o ${qc_dir}/${sample}_R1_trimmed.fq.gz \\
        --trim_poly_x \\
        --adapter_sequence TGGTATCAACGCAGAGTACAT \\
        --report_title ${sample} \\
        --html ${qc_dir}/${sample}.html \\
        --json ${qc_dir}/${sample}.json \\
    2>&1  |
    tee ${qc_dir}/${sample}.log

    echo "03 mapping \`date\` !"

    time ${STAR} \\
        --runMode alignReads \\
        --runThreadN ${thread} \\
        --genomeDir ${index} \\
        --readFilesIn ${qc_dir}/${sample}_R1_trimmed.fq.gz \\
        --readFilesCommand zcat \\
        --outFilterMultimapNmax 1 \\
        --outSAMtype BAM SortedByCoordinate \\
         --outFileNamePrefix ${mapping_dir}/${sample}_ \\
        2>&1 |
    tee ${mapping_dir}/${sample}.log

    time ${perl} ${cell_mapping_ratio} \\
        -w ${thread} \\
        -r1 ${qc_dir}/${sample}_R1_trimmed.fq.gz \\
        -b ${mapping_dir}/${sample}_Aligned.sortedByCoord.out.bam \\
        -o ${mapping_dir}/${sample}_cell_mapping_ratio.tsv \\
        2>&1 |
    tee ${mapping_dir}/${sample}_cell_mapping_ratio.log

    echo "04 count \`date\` !"

    time ${featureCounts} \\
        -T ${thread} \\
        -a ${gtf} \\
        -o ${count_dir}/${sample}_count.tsv \\
        -R BAM \\
        ${mapping_dir}/${sample}_Aligned.sortedByCoord.out.bam \\
        2>&1 |
    tee ${count_dir}/${sample}.log

    time ${samtools} sort \\
        -@ ${thread} \\
        ${count_dir}/${sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam \\
        -o ${count_dir}/${sample}_assigned_sorted.bam
    rm ${count_dir}/${sample}_Aligned.sortedByCoord.out.bam.featureCounts.bam

    echo "05 build matrix \`date\` !"

    time ${perl} ${STRT_umi_count} \\
        -w ${thread} \\
        -b ${count_dir}/${sample}_assigned_sorted.bam \\
        -o ${matrix_dir}/${sample}_count_matrix.tsv \\
        2>&1 |
    tee ${matrix_dir}/${sample}.log

    echo "06 matrix to seurat rds \`date\` !"

    time ${Rscript} ${STRT_cell_mapping_ratio} \\
        ${sample} \\
        ${matrix_dir}/${sample}_count_matrix.tsv \\
        ${barcode} \\
        ${matrix_dir}/${sample}.rds \\ 
        ${mapping_dir}/${sample}_cell_mapping_ratio.tsv \\
        2>&1 |
    tee ${matrix_dir}/${sample}_seurart.log

    echo "end \`date\` !"

} 2>&1 |
tee ${sample_dir}/${sample}_\`date '+%Y%m%d_%H%M%S'\`.log

EOF

    echo "write ${pbs} !"
    echo `date`
    qsub ${pbs}
    echo ""

done 2>&1 | 
tee ${out_dir}/`date '+%Y%m%d_%H%M%S'`.log

