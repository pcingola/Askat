#!/bin/sh

model="Dvar0.01"		# Deleterious (risk)
model="Pvar0.01"		# Protective
model="Mvar0.01"		# Mixed
model="Dvar0.01"		# Deleterious (risk)

dir="t2dgenes_t2dnew"
zcat="gunzip -c"

for model in Dvar0.01 Pvar0.01 Mvar0.01
do
	for gene in ADCY5 CDKAL1 CDKN2B FTO HHEX IGF2BP2 JAZF1 SLC30A8 TCF7L2 WFS1
	do

		#---
		# Output file name
		#---
		vcf=$gene\_$model.vcf
		echo Output file: $vcf

		#---
		# Create header
		#---
		i=1

		case=$dir/$gene/$model/P1_$gene\_$model\_$i.cases.vcf.gz
		control=$dir/$gene/$model/P1_$gene\_$model\_$i.controls.vcf.gz

		$zcat $case | grep "^#" | grep -v "^#CHROM" > $vcf
		$zcat $case | grep "^#CHROM" > header_chr_case
		$zcat $control | grep "^#CHROM" > header_chr_control
		( cat header_chr_control ; cat header_chr_case | cut -f 10- ) | tr -d "\n" >> $vcf
		echo >> $vcf

		#---
		# Add all case / control files
		#---
		for(( i=1 ; i <= 100 ; i++ )); do
			case=$dir/$gene/$model/P1_$gene\_$model\_$i.cases.vcf.gz
			control=$dir/$gene/$model/P1_$gene\_$model\_$i.controls.vcf.gz

			echo "    $i    $case    $control"

			$zcat $control | grep -v "^#" > control.vcf
			$zcat $case | grep -v "^#" | cut -f 10- > case.vcf

			paste control.vcf case.vcf >> $vcf
		done

		#---
		# Add all NULL distribution files
		#---
		for(( i=1 ; i <= 1000 ; i++ )); do
			case=$dir/t2d/$gene/varNULL/P1_$gene\_varNULL\_$i.cases.vcf.gz
			control=$dir/t2d/$gene/varNULL/P1_$gene\_varNULL\_$i.controls.vcf.gz

			echo "    $i    $case    $control"

			$zcat $control | grep -v "^#" > control.vcf
			$zcat $case | grep -v "^#" | cut -f 10- > case.vcf

			paste control.vcf case.vcf >> $vcf
		done

		rm -vf case.vcf control.vcf header_chr_case header_chr_control

		echo Convert positions and IDs
		cat $vcf | ./changeRowNames.pl > new.$vcf
		mv new.$vcf $vcf
	done
done
