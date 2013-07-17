#!/bin/sh

gene=IGF2BP2
model=Dvar0.01

dir=p1simms/t2dgenes_t2dnew

for model in Dvar0.01 Pvar0.01 
do
	#for gene in ADCY5 CDKAL1 CDKN2B FTO HHEX JAZF1 SLC30A8 TCF7L2 WFS1
	#for gene in ADCY5 CDKAL1 CDKN2B FTO HHEX IGF2BP2 JAZF1 SLC30A8 TCF7L2 WFS1
	for gene in ADCY5 TCF7L2 
	do

		geno=$gene\_$model
		vcf=$geno.vcf
		tfam=$geno.tfam

		echo GENOTYPE $geno

		#---
		# Create TFAM
		#---
		cases=$dir/$gene/$model/P1_$gene\_$model\_1.cases.vcf.gz
		controls=$dir/$gene/$model/P1_$gene\_$model\_1.controls.vcf.gz

		# Mark cases as '2' nd controls as '1'
		echo TFAM $tfam
		zcat $controls | grep "^#CHROM" | cut -f 10- | tr "\t" "\n" | grep -v "^$" | sed "s/^/0 /" | sed "s/$/ 0 0 0 1/" > $tfam
		zcat $cases | grep "^#CHROM" | cut -f 10- | tr "\t" "\n" | grep -v "^$" | sed "s/^/0 /" | sed "s/$/ 0 0 0 2/" >> $tfam

		#---
		# Calculate block and subblock size
		#---
		block=`grep -v "^#" $vcf | wc -l`
		subblock=`zcat $cases | grep -v "^#" | wc -l`
		echo "$gene $model:"
		echo "    Block size     : $block"
		echo "    Sub-Block size : $subblock"

		#---
		# Run ASKAT
		#---
		# ( java -Xmx10G -jar Askat.jar -v -p 1 -b $block -kin all -sb $subblock -maxMaf 1.0 $geno 2>&1 | tee $geno.out ) &
		java -Xmx10G -jar Askat.jar -v -p 16 -b $block -kin all -sb $subblock -maxMaf 1.0 $geno 2>&1 | tee $geno.out

	done
	wait
done

