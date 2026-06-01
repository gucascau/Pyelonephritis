#!/bin/bash
#SBATCH --cpus-per-task=32
#SBATCH --time=20:00:00
#SBATCH --mail-user=xin.wang@nationwidechildrens.org
#SBATCH --mail-type=FAIL

set -e

ml purge
ml load spaceranger

Probe=/home/gdbecknelllab/xxw004/Projects/Ruizrosado/Pyelonephritis/Datasets/Spatial_Probe/Visium_Mouse_Transcriptome_Probe_Set_v1.0_mm10-2020-A.csv
DataFold=/home/gdbecknelllab/xxw004/Projects/Ruizrosado/Pyelonephritis/Datasets/221215_Ruiz-Rosado_GSL-JDDRR-3171

cd /home/gdbecknelllab/xxw004/Projects/Ruizrosado/Pyelonephritis/Results/Pyelonephritis_31

spaceranger count --id=Pyelonephritis_31 --transcriptome=/reference/spaceranger/refdata-gex-mm10-2020-A --probe-set=$Probe --fastqs=$DataFold/31 --sample=31 --image=$DataFold/VisiumSlides/Slide3/31_Kidney1dpi.tif  --slide=V12U01-134 --area=A1 --reorient-images --localcores=32 --localmem=240
