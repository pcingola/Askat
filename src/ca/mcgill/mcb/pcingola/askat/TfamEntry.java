package ca.mcgill.mcb.pcingola.askat;

import ca.mcgill.mcb.pcingola.ped.Sex;
import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * An entry in a TFAM file
 * @author pablocingolani
 */
public class TfamEntry {

	// All members are public, since the format is not likely to change 
	public String id, familyId, paternalId, maternalId;
	public Sex sex;
	public double phenotype;

	public TfamEntry() {
		id = familyId = paternalId = maternalId = "";
		sex = Sex.Unknown;
		phenotype = 0;
	}

	public TfamEntry(String tfamLine) {
		parse(tfamLine);
	}

	public boolean isPhenotypeAffected() {
		return phenotype == 2;
	}

	/**
	 * Is phenotype missing?
	 * @return
	 */
	public boolean isPhenotypeMissing() {
		return phenotype <= 0;
	}

	public boolean isPhenotypeUnaffected() {
		return phenotype == 1;
	}

	void parse(String tfamLine) {
		String fields[] = tfamLine.split("\\s", -1); // Split, but do not discard trailing empty fields
		familyId = fields[0];
		id = fields[1];
		paternalId = fields[2];
		maternalId = fields[3];

		// Coding sex: 	1=male; 2=female; other=unknown
		String sexStr = fields[4];
		if (sexStr.equals("1")) sex = Sex.Male;
		else if (sexStr.equals("2")) sex = Sex.Female;
		else sex = Sex.Unknown;

		// Phenotype
		phenotype = Gpr.parseDoubleSafe(fields[5]);
	}

	@Override
	public String toString() {
		String sexStr = "";
		if (sex == Sex.Male) sexStr = "1";
		else if (sex == Sex.Female) sexStr = "2";
		else sexStr = "0";

		return familyId + " " + id + " " + paternalId + " " + maternalId + " " + sexStr + " " + phenotype;
	}
}
