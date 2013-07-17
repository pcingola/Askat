package ca.mcgill.mcb.pcingola.askat;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import ca.mcgill.mcb.pcingola.util.Gpr;

/**
 * A TFAM file: 
 * 
 * References: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
 * 
 * @author pablocingolani
 */
public class Tfam implements Iterable<TfamEntry> {

	ArrayList<TfamEntry> tfamEntries;

	public Tfam() {
		tfamEntries = new ArrayList<TfamEntry>();
	}

	public Tfam(String tfamFileName) {
		load(tfamFileName);
	}

	public void add(TfamEntry te) {
		tfamEntries.add(te);
	}

	/**
	 * Find a TfamEntry by id
	 * @param id
	 * @return
	 */
	public TfamEntry find(String id) {
		for (TfamEntry te : this)
			if (te.id.equals(id)) return te;
		return null;
	}

	/**
	 * Get TfamEntry number 'idx'
	 * @param idx
	 * @return
	 */
	public TfamEntry getEntry(int idx) {
		return tfamEntries.get(idx);
	}

	/**
	 * Get a list of sample names
	 * @return
	 */
	public List<String> getSampleNames() {
		ArrayList<String> sampleNames = new ArrayList<String>();
		for (TfamEntry te : this)
			sampleNames.add(te.id);
		return sampleNames;
	}

	@Override
	public Iterator<TfamEntry> iterator() {
		return tfamEntries.iterator();
	}

	/**
	 * Load: Read and parse files
	 * @param tfamFileName
	 */
	void load(String tfamFileName) {
		// Read file
		if (!Gpr.canRead(tfamFileName)) throw new RuntimeException("Cannot read TFAM file '" + tfamFileName + "'");
		String lines[] = Gpr.readFile(tfamFileName).split("\n");

		// Add all entries
		tfamEntries = new ArrayList<TfamEntry>();
		for (String line : lines) {
			TfamEntry te = new TfamEntry(line);
			tfamEntries.add(te);
		}
	}

	public void save(String fileName) {
		StringBuilder sb = new StringBuilder();
		for (TfamEntry te : this)
			sb.append(te + "\n");

		Gpr.toFile(fileName, sb.toString());
	}

	public int size() {
		return tfamEntries.size();
	}

	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (TfamEntry te : this)
			sb.append(te + "\n");
		return sb.toString();
	}
}
