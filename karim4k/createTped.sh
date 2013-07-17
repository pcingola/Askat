#!/bin/sh

# Create PED
cat SimDATA4000SNPs.dat | cut -f 6- -d " " > snps.txt
paste SimulDATA_pheno.dat snps.txt | tr "\t" " " | tr -s " " | ./dat2ped.pl  > karim4k.ped
./createMap.pl > karim4k.map


# Convert to TPED
$HOME/tools/plink/plink --noweb --recode --transpose --file karim4k
mv plink.tfam karim4k.tfam
mv plink.tped karim4k.tped
