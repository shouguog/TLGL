#! /bin/bash
# this file is cellranger_batch.sh
#SBATCH -J AGGREGATION
ulimit -u 10240
ulimit -n 16384
#sbatch --cpus-per-task=40  --time=120:00:00 --mem=56g VDJ_GRCh38.sh
module load cellranger/2.1.1 || exit 1
#cellranger count --help
for ii in {32..32}
do
	echo "start sample"
	echo ${ii}
	cellranger vdj --id=SNOLANE1_${ii}_VDJT --fastqs=../180427_J00139_0307_AHV3YCBBXX/HV3YCBBXXNOLANE8/outs/fastq_path/LBL_VDJ/S${ii},../180309_J00139_0286_AHTCMGBBXX/HTCMGBBXXNEWNOLANE1_SIGA_FINAL/outs/fastq_path/LBL_VDJ/S${ii} --localmem=48 --reference=/data/gaos2/10xgenomics/refdata-cellranger-vdj-GRCh38-alts-ensembl-2.0.0 --sample=S${ii}
done


