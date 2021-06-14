#! /bin/bash
# this file is cellranger_batch.sh
#SBATCH -J AGGREGATION
ulimit -u 10240
ulimit -n 16384
#sbatch --cpus-per-task=20  --time=60:00:00 --mem=56g count_LBL_5prime.sh
module load cellranger/2.1.1 || exit 1
cellranger count --help
for ii in {4..4}
do
	echo "start sample"
	echo ${ii}
	cellranger count --id=SfivePrimeR2GRCh38_${ii} --fastqs=\
../180227_J00139_0282_AHMWWJBBXX/HMWWJBBXXv211wolfgang/outs/fastq_path/LBL/S${ii},\
../180323_J00139_0293_AHMYV3BBXX/HMYV3BBXX/outs/fastq_path/LBL/S${ii},\
../180326_J00139_0294_AHNGW3BBXX/HNGW3BBXX/outs/fastq_path/LBL/S${ii},\
../180328_J00139_0295_AHNGLHBBXX/HNGLHBBXX/outs/fastq_path/LBL/S${ii},\
../180402_J00139_0297_AHT373BBXX/HT373BBXX/outs/fastq_path/LBL/S${ii},\
../180409_J00139_0299_AHV52NBBXX/HV52NBBXX/outs/fastq_path/LBL/S${ii},\
../180419_J00139_0305_AHV3TCBBXX/HV3TCBBXX/outs/fastq_path/LBL/S${ii} \
--localmem=48 --transcriptome=/data/gaos2/10xgenomics/refdata-cellranger-GRCh38-1.2.0 --chemistry=SC5P-R2
done


