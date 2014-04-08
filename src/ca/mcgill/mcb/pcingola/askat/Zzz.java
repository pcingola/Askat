package ca.mcgill.mcb.pcingola.askat;

import ca.mcgill.mcb.pcingola.askat.Askat.UseMissing;

public class Zzz {

	public static void main(String[] args) {
		Askat askat = new Askat(args);
		String vcfFile = "/home/pcingola/askat/lof.private.vcf";
		String tpedFile = "/home/pcingola/askat/lof.private.tped";
		String tfamFile = "/home/pcingola/askat/lof.private.tfam";

		askat.setVerbose(true);
		askat.setDebug(true);
		askat.setUseMissing(UseMissing.REFERENCE);
		askat.setTfamFile(tfamFile);

		askat.vcf2Tped(vcfFile, tpedFile);
	}

}
