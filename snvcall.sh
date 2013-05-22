#!/bin/bash
set -v
nice samtools index HSS41_GATK.bam &
nice samtools index HSS43_GATK.bam &
wait
nice bam-somaticsniper -q 1 -Q 3 -J -s 0.0001 -f hs37d5.fa HSS41_GATK.bam HSS43_GATK.bam ssp1diagraw &
nice samtools mpileup -q 1 -f hs37d5.fa HSS41_GATK.bam > HSS41_GATK.pileup &
nice samtools mpileup -q 1 -f hs37d5.fa HSS43_GATK.bam > HSS43_GATK.pileup &
nice jsm.py train joint_snv_mix_two hs37d5.fa HSS43_GATK.bam HSS41_GATK.bam /media/Data/JointSNVMix-0.7.5/config/joint_priors.cfg /media/Data/JointSNVMix-0.7.5/config/joint_params.cfg jsm2p1diagparams &
wait
nice jsm.py classify joint_snv_mix_two hs37d5.fa HSS43_GATK.bam HSS41_GATK.bam jsm2p1diagparams jsm2p1rawcalls &
nice java -jar /media/Data/VarScan.v2.2.11.jar somatic HSS43_GATK.pileup HSS41_GATK.pileup vsp1diag --somatic-p-value 0.5 --strand-filter 1 &
wait
nice java -jar /media/Data/VarScan.v2.2.11.jar processSomatic vsp1diag.snp &
wait
nice awk '$9<=0.5' jsm2p1rawcalls > jsm2p1soi.txt &
nice /media/Data/strelka_workflow/configureStrelkaWorkflow.pl --normal=HSS43_GATK.bam --tumor=HSS41_GATK.bam --ref=hs37d5.fa --config=config.ini --output-dir=./StrelkaP1 &
wait
