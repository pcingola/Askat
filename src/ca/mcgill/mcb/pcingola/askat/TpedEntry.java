package ca.mcgill.mcb.pcingola.askat;

import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.GprSeq;

/**
 * Entry in a TPED file (a line)
 * 
 * @author pcingola
 */
public class TpedEntry extends Marker {

	private static final long serialVersionUID = 3407992013631475007L;

	Genome genome;
	String line;
	String recs[];
	int count[];
	char genotypes[];
	char maxBase;

	public TpedEntry(Genome genome, String line) {
		this.genome = genome;
		this.line = line;
		parse();
	}

	/**
	 * Calculate Minimum allele frequency
	 * @param line
	 * @return
	 */
	public double maf() {
		double maf = 1.0;

		// Total
		int tot = 0;
		for (int i = 0; i < count.length; i++)
			tot += count[i];

		// Minor allele frequency
		for (int i = 0; i < count.length; i++) {
			if (count[i] > 0) {
				double af = ((double) count[i]) / ((double) tot);
				maf = Math.min(maf, af);
			}
		}

		return maf;
	}

	/**
	 * Parse a TPED line. Fill array of genotypes and counts
	 * @param line
	 * @param genotypes
	 * @param count
	 * @return Base having max frequency
	 */
	void parse() {
		recs = line.split("\\s+");
		if (recs.length % 2 != 0) throw new RuntimeException("Odd number of records. This should never happen!\n\t'" + line + "'");

		// Chromosome
		parent = genome.getOrCreateChromosome(recs[0]);

		// Start and End position
		start = end = Gpr.parseIntSafe(recs[3]);

		// Count bases
		count = new int[4];
		genotypes = new char[recs.length - 4];

		// Convert to char and count
		for (int i = 4, j = 0; i < recs.length; i++, j++) {
			char base = Character.toUpperCase(recs[i].charAt(0));
			genotypes[j] = base;

			switch (base) {
			case 'A':
				count[0]++;
				break;
			case 'C':
				count[1]++;
				break;
			case 'G':
				count[2]++;
				break;
			case 'T':
				count[3]++;
				break;
			default: // Nothing to do (missing data)
			}
		}

		// Get base having maximum frequency
		int max = 0, maxIdx = 0;
		for (int i = 0; i < count.length; i++) {
			if (max < count[i]) {
				max = count[i];
				maxIdx = i;
			}
		}

		// Major allele
		maxBase = GprSeq.BASES[maxIdx]; // This is assumed to be the reference
	}

	/**
	 * Transform a TPED line into a data line usable by askat
	 * @param line
	 * @return
	 */
	public String tped2askatDat() {
		// Create a numeric line
		StringBuilder datLine = new StringBuilder();
		datLine.append(recs[0] + " ");
		datLine.append(recs[1] + " ");
		datLine.append(recs[2] + " ");
		datLine.append(recs[3] + " ");
		for (int i = 0; i < genotypes.length; i += 2) {
			int num0 = (genotypes[i] == maxBase ? 0 : 1);
			int num1 = (genotypes[i + 1] == maxBase ? 0 : 1);

			int num = num0 + num1;
			datLine.append(num + " ");
		}

		return datLine.toString();
	}

}
