#!/bin/bash
export PATH=~/software_database/rnaseq:$PATH
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334276_1.fastq -S ./dried_bam/dried_gills_1d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334277_1.fastq -S ./dried_bam/dried_gills_3d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334278_1.fastq -S ./dried_bam/dried_gills_5d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334279_1.fastq -S ./dried_bam/dried_gills_7d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334280_1.fastq -S ./dried_bam/dried_gills_9d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334281_1.fastq -S ./dried_bam/dried_gills_10d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334282_1.fastq -S ./dried_bam/dried_gills_11d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334283.fastq -S ./dried_bam/dried_muscle_0d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334284.fastq -S ./dried_bam/dried_muscle_1d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334285.fastq -S ./dried_bam/dried_muscle_3d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334286.fastq -S ./dried_bam/dried_muscle_5d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334287.fastq -S ./dried_bam/dried_muscle_7d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334288.fastq -S ./dried_bam/dried_muscle_9d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334289.fastq -S ./dried_bam/dried_muscle_10d.sam
hisat2 -p 32 -x ./genome/oyster_genome -U ./SRR/SRR334290.fastq -S ./dried_bam/dried_muscle_11d.sam
echo "hisat2 complete"
cd dried_bam
samtools view -bS dried_gills_1d.sam >unsorted_dried_gills_1d.bam
samtools view -bS dried_gills_3d.sam >unsorted_dried_gills_3d.bam
samtools view -bS dried_gills_5d.sam >unsorted_dried_gills_5d.bam
samtools view -bS dried_gills_7d.sam >unsorted_dried_gills_7d.bam
samtools view -bS dried_gills_9d.sam >unsorted_dried_gills_9d.bam
samtools view -bS dried_gills_10d.sam >unsorted_dried_gills_10d.bam
samtools view -bS dried_gills_11d.sam >unsorted_dried_gills_11d.bam
samtools view -bS dried_muscle_0d.sam >unsorted_dried_muscle_0d.bam
samtools view -bS dried_muscle_1d.sam >unsorted_dried_muscle_1d.bam
samtools view -bS dried_muscle_3d.sam >unsorted_dried_muscle_3d.bam
samtools view -bS dried_muscle_5d.sam >unsorted_dried_muscle_5d.bam
samtools view -bS dried_muscle_7d.sam >unsorted_dried_muscle_7d.bam
samtools view -bS dried_muscle_9d.sam >unsorted_dried_muscle_9d.bam
samtools view -bS dried_muscle_10d.sam >unsorted_dried_muscle_10d.bam
samtools view -bS dried_muscle_11d.sam >unsorted_dried_muscle_11d.bam
echo "samtools view complete"
samtools sort -@ 32 unsorted_dried_gills_1d.bam dried_gills_1d
samtools sort -@ 32 unsorted_dried_gills_3d.bam dried_gills_3d
samtools sort -@ 32 unsorted_dried_gills_5d.bam dried_gills_5d
samtools sort -@ 32 unsorted_dried_gills_7d.bam dried_gills_7d
samtools sort -@ 32 unsorted_dried_gills_9d.bam dried_gills_9d
samtools sort -@ 32 unsorted_dried_gills_10d.bam dried_gills_10d
samtools sort -@ 32 unsorted_dried_gills_11d.bam dried_gills_11d
samtools sort -@ 32 unsorted_dried_muscle_0d.bam dried_muscle_0d
samtools sort -@ 32 unsorted_dried_muscle_1d.bam dried_muscle_1d
samtools sort -@ 32 unsorted_dried_muscle_3d.bam dried_muscle_3d
samtools sort -@ 32 unsorted_dried_muscle_5d.bam dried_muscle_5d
samtools sort -@ 32 unsorted_dried_muscle_7d.bam dried_muscle_7d
samtools sort -@ 32 unsorted_dried_muscle_9d.bam dried_muscle_9d
samtools sort -@ 32 unsorted_dried_muscle_10d.bam dried_muscle_10d
samtools sort -@ 32 unsorted_dried_muscle_11d.bam dried_muscle_11d
echo "samtools sort complete"