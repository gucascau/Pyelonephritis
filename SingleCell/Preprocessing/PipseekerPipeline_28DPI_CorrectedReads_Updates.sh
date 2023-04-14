#!/bin/sh
#SBATCH --cpus-per-task=16
#SBATCH --mail-type=FAIL,REQUEUE
#SBATCH --mail-user=xin.wang@nationwidechildrens.org
#SBATCH --job-name=Pipseeker
#SBATCH --time=40:00:00
#SBATCH --ntasks=1


set -e

ml purge

ml STAR/2.7.9a

id=28DPI

pipseeker=/home/gdbecknelllab/xxw004/Software/pipseeker-v1.1.7-linux-X86
starindex=/home/gdbecknelllab/xxw004/Datasets/pipseeker-gex-reference-GRCm39-2022.04
script=/home/gdbecknelllab/xxw004/Projects/Ruizrosado/Pyelonephritis/PIPseqScripts

in=/home/gdbecknelllab/xxw004/Projects/Ruizrosado/Pyelonephritis/Datasets/221221_Ruiz-Rosado_GSL-JDDRR-3176/${id}
out=/home/gdbecknelllab/xxw004/Projects/Ruizrosado/Pyelonephritis/PIPseqResults/${id}



cd ${in}

#gunzip *.gz

perl ${script}/ExtractSequenceCorrectedBaseOnLength.pl -i 51 -f ${id}_S4_L001_R1_001.fastq -r ${id}_S4_L001_R2_001.fastq -o ${id}_Cfiltered_L001
#
perl ${script}/ExtractSequenceCorrectedBaseOnLength.pl  -i 51 -f ${id}_S4_L002_R1_001.fastq -r ${id}_S4_L002_R2_001.fastq -o ${id}_Cfiltered_L002
#
#
# gzip *.fastq

mkdir -p ${out}
cd ${out}

${pipseeker} full --fastq ${in}/${id}_Cfiltered  --output-path ${out} --star-index-path ${starindex}
