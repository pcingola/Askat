package ca.mcgill.mcb.pcingola.askat;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.collections.MultivalueHashMap;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.interval.Marker;
import ca.mcgill.mcb.pcingola.interval.Markers;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.interval.tree.IntervalForest;
import ca.mcgill.mcb.pcingola.osCmd.LineFilter;
import ca.mcgill.mcb.pcingola.osCmd.OsCmdQueue;
import ca.mcgill.mcb.pcingola.osCmd.OsCmdRunner;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;

/**
 * A class representing a "block" of SNPs to calculate kinship matrix and then call ASKAT's R script
 * 
 * @author pcingola
 */
public class KinshipBlock {

	public static boolean debug = false; // Debug mode
	public static boolean debugOnlyOnce = false; // Debug mode, just do one iteration

	String blockFile, blockName, genabelGenFile, genabelPhenFile, kinshipFile, simFile, phenoFile;
	HashSet<String> toDelete;
	List<SeqChange> intervals;
	Askat askat;

	public KinshipBlock(Askat askat, String blockFile) {
		this.askat = askat;
		this.blockFile = blockFile;
		blockName = Gpr.removeExt(blockFile);
		genabelGenFile = blockName + ".genabel.gen";
		genabelPhenFile = blockName + ".genabel.phen";
		kinshipFile = blockName + ".kinship.RData";
		simFile = blockName + ".sim";
		phenoFile = blockName + ".pheno.txt";
		toDelete = new HashSet<String>();
	}

	/**
	 * Run ASKAT on each sub-block
	 */
	public void askat() {
		if (askat.isVerbose()) Timer.showStdErr("Starting block: " + blockName);

		// Split workload into batches
		List<String> batchFiles = batchFiles();

		// Create command and add them to a queue
		OsCmdQueue queue = createQueue(batchFiles);

		// Run queue
		queue.run();

		// Delete all tmp files & directories
		if (!debug) deleteFiles();

		// OK, we are done
		if (askat.isVerbose()) Timer.showStdErr("Finished block: " + blockName);
	}

	/**
	 * Split block into similar sized "batch files"
	 * @return
	 */
	List<String> batchFiles() {
		if (intervals == null) return batchFilesFixedSize(); // Use a fixed number of SNPs in each batch
		return batchFilesIntervals(); // Use a intervals to split files		
	}

	/**
	 * Split block into similar sized "batch files"
	 * @return
	 */
	List<String> batchFilesFixedSize() {
		List<String> batchFiles = new ArrayList<String>();

		// Estimate batch size
		if (askat.isVerbose()) Timer.showStdErr("Counting lines in file '" + blockFile + "'");
		int numLines = Gpr.countLines(blockFile);
		int subBlock = askat.getSubBlockSize();
		int numLinesPerCpu = (numLines / (subBlock * askat.getNumWorkers())) * subBlock;
		int batchLines = Math.max(numLinesPerCpu, subBlock);
		if (askat.isVerbose()) Timer.showStdErr("Create batches.\n\t\t\tFile '" + blockFile + "' has " + numLines + " lines.\n\t\t\tSplit up to " + batchLines + " lines per batch.");

		// Create batches 
		try {
			String batchFileName = null;
			int batchNum = 1, lineNum = 0;
			BufferedWriter outFile = null;
			LineFileIterator lfi = new LineFileIterator(blockFile);

			for (String line : lfi) {
				// Transform from TPED to ASKAT
				TpedEntry tpedEntry = new TpedEntry(askat.getGenome(), line);
				line = tpedEntry.tped2askatDat();

				if ((outFile == null) || (lineNum >= batchLines)) {
					if (outFile != null) outFile.close();

					batchFileName = blockName + "." + batchNum + ".askat";
					batchFiles.add(batchFileName);
					if (askat.isVerbose()) Timer.showStdErr("Batch " + batchNum + ". Line " + lfi.getLineNum() + ". Creating batch : " + batchFileName);
					outFile = new BufferedWriter(new FileWriter(batchFileName));
					batchNum++;
					lineNum = 0;
				}

				outFile.write(line + "\n");
				lineNum++;
			}

			if (outFile != null) outFile.close();

		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		return batchFiles;
	}

	/**
	 * Split block into according to 'intervals'
	 * @return
	 */
	List<String> batchFilesIntervals() {
		// Build forest
		if (askat.isVerbose()) Timer.showStdErr("Building interval forest.");
		IntervalForest intervalForest = new IntervalForest();
		intervalForest.add(intervals);
		intervalForest.build();
		if (askat.isVerbose()) Timer.showStdErr("done.");

		List<String> batchFiles = new ArrayList<String>();

		// Map lines to intervals
		LineFileIterator lfi = new LineFileIterator(blockFile);
		MultivalueHashMap<Marker, String> interval2tped = new MultivalueHashMap<Marker, String>();
		for (String line : lfi) {
			// Transform from TPED to ASKAT
			TpedEntry tpedEntry = new TpedEntry(askat.getGenome(), line);
			line = tpedEntry.tped2askatDat();

			// See if entry hits ANY interval in intervalForest
			Markers results = intervalForest.query(tpedEntry);

			// Show warning if no interval is hit
			if (askat.isDebug() && results.isEmpty()) {
				if (debug) System.err.println("WARNING: TPED entry '" + tpedEntry.getChromosomeName() + ":" + tpedEntry.getStart() + "' did not hit any interval. Ignored.");
			}

			// Map entry to intervals it hits (can be more then one)
			HashSet<Marker> markersNotHit = new HashSet<Marker>();
			markersNotHit.addAll(interval2tped.keySet());
			for (Marker m : results) {
				interval2tped.add(m, tpedEntry.tped2askatDat());
				markersNotHit.remove(m);
			}

			// Data is assumed to be sorted by position. So, if a marker is 
			// no longer hit, we can assume that it will no be hit again.
			// So we can save the data for the markers that were NOT hit.
			for (Marker m : markersNotHit) {
				String batchFile = saveFileMarker(m, interval2tped); // Save file
				if (batchFile != null) batchFiles.add(batchFile); // Add file and increment number
				interval2tped.remove(m); // Remove marker (we are done)
			}
		}

		// Save all files that have not been saved so far
		for (Marker m : interval2tped.keySet()) {
			String batchFile = saveFileMarker(m, interval2tped); // Save file
			if (batchFile != null) batchFiles.add(batchFile); // Add file and increment number
			interval2tped.remove(m); // Remove marker (we are done)
		}

		return batchFiles;
	}

	/**
	 * Save all entries that march an interval (marker)
	 * @param m
	 * @param interval2tped
	 * @return
	 */
	String saveFileMarker(Marker m, MultivalueHashMap<Marker, String> interval2tped) {
		// Save file
		String mid = m.getId().replaceAll("[^a-zA-Z0-9\\-\\.]+", "_");
		String batchFile = blockName + "." //
				+ m.getChromosomeName() //
				+ ":" + (m.getStart() + 1) //
				+ "-" + (m.getEnd() + 1) //
				+ "_" + mid //
				+ ".askat";

		// Get lines
		List<String> lines = interval2tped.get(m);

		// Save file
		if (intervalsCreateFile(lines, batchFile, m)) return batchFile;
		return null;
	}

	/**
	 * Create a queue of commands to be execeuted by the OS
	 * @param batchFiles
	 * @return
	 */
	OsCmdQueue createQueue(List<String> batchFiles) {
		OsCmdQueue queue = new OsCmdQueue();
		queue.setNumThreads(askat.getNumWorkers());
		queue.setVerbose(askat.isVerbose());

		if (intervals == null) {
			// Block & sub-block method
			for (String batchFile : batchFiles) {
				toDelete.add(batchFile);
				createQueueAdd(queue, batchFile);
				if (debugOnlyOnce) break;
			}
		} else {
			// Intervals method
			int filesPerWorker = batchFiles.size() / askat.getNumWorkers() + 1;
			if (askat.isVerbose()) Timer.showStdErr("Using " + filesPerWorker + " files per worker.");

			// Create a list of batchFile separated by comma
			int count = 1, countJob = 0;
			StringBuilder batchFilesPerWorker = new StringBuilder();
			for (String batchFile : batchFiles) {
				batchFilesPerWorker.append((batchFilesPerWorker.length() > 0 ? "," : "") + batchFile);
				countJob++;

				if (count % filesPerWorker == 0) {
					// Create a job
					createQueueAdd(queue, batchFilesPerWorker.toString());
					if (askat.isVerbose()) Timer.showStdErr("\tAdded job " + queue.size() + ". Number of files : " + countJob);
					batchFilesPerWorker = new StringBuilder();
					if (debugOnlyOnce) break;
					countJob = 0;
				}

				count++;
			}

			// Submit last job
			if (batchFilesPerWorker.length() > 0) {
				createQueueAdd(queue, batchFilesPerWorker.toString());
				if (askat.isVerbose()) Timer.showStdErr("\tAdded job " + queue.size() + ". Number of files : " + countJob);
			}
		}

		return queue;
	}

	/**
	 * Create a command and add it to the queue
	 * @param queue
	 * @param batchFile
	 */
	void createQueueAdd(OsCmdQueue queue, String batchFile) {

		// Command line arguments to execute ASKAT R script 
		String args[] = { askat.getPath(Askat.CMD_RSCRIPT) //
				, askat.getrPath() + Askat.R_SCRIPT_ASKAT //
				, batchFile //
				, askat.tfamFile //
				, kinshipFile //
				, askat.getSubBlockSize() + "" //
				, askat.getpACC() + "" //UPD new command line argument for the p-value accuracy
				, Boolean.toString(debugOnlyOnce).toUpperCase() //
		};

		// Create a line filter
		LineFilter lineFilter = new LineFilter() {

			@Override
			public String filter(String line) {
				if (line.startsWith(Askat.ASKAT_RESULTS)) return line;
				if (line.startsWith(Askat.ASKAT_WARNING)) return line; //Warning in case p-value is not converged
				return null;
			}
		};

		// Create command
		String rScriptName = args[1];
		OsCmdRunner rScriptCmd = new OsCmdRunner("R_Script_" + rScriptName, args);
		//rScriptCmd.getOsCmd().setQuiet(!debug, !debug);
		rScriptCmd.getOsCmd().setQuiet(false, !askat.isDebug());
		rScriptCmd.getOsCmd().setSaveStd(true); // We want to save and parse STDOUT
		if (!askat.isDebug()) rScriptCmd.getOsCmd().setStdOutFilter(lineFilter);
		queue.add(rScriptCmd);
	}

	/**
	 * Delete tmp files and dirs
	 */
	void deleteFiles() {
		if (debug) {
			Timer.showStdErr("Debug mode: No files deleted!");
			return;
		}

		for (String fname : toDelete) {
			File f = new File(fname);
			if (askat.isVerbose()) Timer.showStdErr("Deleting '" + fname + "'");
			if (!f.delete()) System.err.println("Cannot delete '" + f.getAbsolutePath() + "'");
		}
	}

	/**
	 * Save file for batchFilesIntervals
	 * @param lines
	 * @param batchFile
	 * @param m
	 * @return
	 */
	boolean intervalsCreateFile(List<String> lines, String batchFile, Marker m) {
		if (lines != null) {
			// Don't create file if it's less than MinVariants
			if (askat.getMinVariants() >= lines.size()) {
				if (askat.isVerbose()) Timer.showStdErr("Interval " + m + "only has " + lines.size() + " variants. Skipping.");
				return false;
			}

			StringBuilder sb = new StringBuilder();
			for (String line : lines)
				sb.append(line + "\n");

			// Save file
			if (sb.length() > 0) {
				if (askat.isVerbose()) Timer.showStdErr("Saving " + lines.size() + " lines to file '" + batchFile + "' corresponding to interval " + m);
				Gpr.toFile(batchFile, sb);
				return true;
			}
		} else if (askat.isVerbose()) Timer.showStdErr("Interval " + m + " has no variants: Skipped.");
		return false;
	}

	/**
	 * Execute an R script to calculate kinship matrix
	 * @param cmd
	 */
	public void kinship() {
		if (askat.isVerbose()) Timer.showStdErr("Calculating kinship matrix for block: " + blockName);

		// Kinship file already exists? Use it!
		if (Gpr.canRead(kinshipFile)) {
			if (askat.isVerbose()) Timer.showStdErr("Kinship file '" + kinshipFile + "' alrady exists. Nothing done.");
			return;
		}

		// We should delete all these files after we are done
		toDelete.add(genabelGenFile);
		toDelete.add(genabelPhenFile);
		toDelete.add(simFile);
		toDelete.add(phenoFile);
		// These two files are expensive to calculate. We keep them in case we want to re-run
		//		toDelete.add(blockFile);
		//		toDelete.add(kinshipFile);

		// Create command line and call kinship R script
		String cmd[] = { askat.getPath(Askat.CMD_RSCRIPT) //
				, askat.getrPath() + Askat.R_SCRIPT_KINSHIP //
				, blockFile //
				, askat.tfamFile //
				, genabelGenFile //
				, genabelPhenFile //
				, kinshipFile //
				, simFile //
				, phenoFile //
				, askat.getPath(Askat.CMD_FASTLMM) //
		};

		String rScriptName = cmd[1];
		OsCmdRunner rScriptCmd = new OsCmdRunner("R_Script_" + rScriptName, cmd);
		rScriptCmd.getOsCmd().setQuiet(!debug, !debug);
		if (debug) Timer.showStdErr("\tExecuting command: " + rScriptCmd.getOsCmd());
		rScriptCmd.run();
		if (rScriptCmd.getExitValue() != 0) askat.fatalError("Execution of R script '" + rScriptName + "' failed.\n\tCommand line: " + rScriptCmd);
	}

	public void setIntervals(List<SeqChange> intervals) {
		this.intervals = intervals;
	}
}
