package ca.mcgill.mcb.pcingola.askat;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import ca.mcgill.mcb.pcingola.Pcingola;
import ca.mcgill.mcb.pcingola.fileIterator.BedFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.LineFileIterator;
import ca.mcgill.mcb.pcingola.fileIterator.VcfFileIterator;
import ca.mcgill.mcb.pcingola.interval.Genome;
import ca.mcgill.mcb.pcingola.interval.SeqChange;
import ca.mcgill.mcb.pcingola.osCmd.OsCmdRunner;
import ca.mcgill.mcb.pcingola.snpEffect.commandLine.CommandLine;
import ca.mcgill.mcb.pcingola.stats.CountByType;
import ca.mcgill.mcb.pcingola.util.Gpr;
import ca.mcgill.mcb.pcingola.util.Timer;
import ca.mcgill.mcb.pcingola.vcf.VcfEntry;
import ca.mcgill.mcb.pcingola.vcf.VcfGenotype;





/**
 * Command line wrapper to ASKAT method
 * 
 * References: "Adjusted Sequence Kernel Association Test for Rare Variants Controlling for Cryptic and Family Relatedness"
 * 				Karim Oualkacha & Celia Greenwood
 * 
 * What does the wrapper do:
 *   1- Make sure there are two files "genotype.tped" and "genotype.tfam" (or even a genotype.vcf)
 *   2- Make sure that all dependencies are installed
 *   3- Split input file in "blocks" of SNPs (kinship calculation)
 *   	For each chunk:
 *   		3.1- Calculate the kinship matrix
 *   		3.2- Create sub-blocks (typically consisting of 20 SNPs)
 *   		3.3- Create a file in an intermediate 'askat' format
 *   		3.4- Call call ASKAT R function (in parallel) using the pre-calculated kinship matrix  			
 *   		3.5- Collect the results 
 *   4- Show summary
 *
 * References about PED, MAP, TPED and TFAM formats: http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml
 *
 * @author pablocingolani
 */
public class Askat implements CommandLine {

	// Method available to calculate kinship matrix
	public enum KinshipMethod {
		BLOCK, CHROMOSOME, CHROMOSOME_AVG, ALL
	}

	// What to do with mssing genotypes
	public enum UseMissing {
		DO_NOT_USE // Do not use: Lines are ommited
		, MISSING // Mark as missing in TPED file
		, REFERENCE // Mark as reference in TPED file
	};

	public static final int VERY_LARGE_BLOCK_SIZE = 1000 * 1000 * 1000;

	// Version
	public static final String BUILD = "2013-09-05";
	public static final String VERSION_MAJOR = "1.2";
	public static final String REVISION = "d";
	public static final String VERSION_SHORT = VERSION_MAJOR + REVISION;
	public static final String VERSION = VERSION_SHORT + " (build " + BUILD + "), by " + Pcingola.BY;

	// Programs
	public static final String CMD_RSCRIPT = "Rscript";
	public static final String CMD_R = "R";
	public static final String CMD_FASTLMM = "fastlmmc";

	// R scripts
	public static final String R_SCRIPT_KINSHIP = "kinship.r";
	public static final String R_SCRIPT_ASKAT = "askat.r";
	public static final String R_SCRIPT_ERROR_TEST = "errorTest.r";
	public static final String DEPENDECY_SCRIPTS[] = { R_SCRIPT_KINSHIP, R_SCRIPT_ASKAT };

	// R Libraries
	public static final String DEPENDECY_RLIBS[] = { "GenABEL", "CompQuadForm", "nFactors", "MASS" };

	// Askat results identifiers (askat.r)
	public static final String ASKAT_RESULTS = "ASKAT_RESULTS:";
	public static final String ASKAT_WARNING = "WARNING:";

	protected int numWorkers = Gpr.NUM_CORES; // Max number of threads (if multi-threaded version is available)
	boolean debug = false; // Debug mode
	protected boolean verbose = false; // Be verbose
	protected boolean dependencyCheck = true; // Perform a dependency check
	boolean onlySnp = false; // Only use SNPs in VCF files
	protected String[] args;
	protected String genotypeName = "";
	protected String tpedFile;
	protected String tfamFile;
	protected String rPath = "./r/"; // Path to ASKAT R scripts. Note: It must end with '/'
	protected String binPath = "./"; // Path to binary programs. Note: It must end with '/'
	protected String bedFile = null; // BED file with intervals 
	protected int blockSize = VERY_LARGE_BLOCK_SIZE; // Block size: Number of SNPs used to calculate the kinship matrix (default: large number)
	protected int subBlockSize = 20; // Sub-block: Number of SNPs used in each call to ASKAT
	protected int minVariants = 3; // Don't use less than this number of variants
	protected UseMissing useMissing = UseMissing.DO_NOT_USE; // Do not use genotypes having missing values
	protected Genome genome;
	protected HashMap<String, String> pathToBin;
	protected KinshipMethod kinshipMethod = KinshipMethod.CHROMOSOME;
	protected double maxMaf = 1.0; // Maximum 'MAF' allowed for the analysis (filter out other SNPs).
	protected double pACC = 1e-9; // accuracy parameter for the r-method 'davies' computing p-value
	protected List<SeqChange> intervals;

	public static void main(String[] args) {
		Askat askat = new Askat(args);
		askat.parseArgs(args);
		askat.run();
	}

	public Askat(String[] args) {
		genome = new Genome();
		pathToBin = new HashMap<String, String>();
		this.args = args;
	}

	void askat(String blockFile, String sub) {
		String cmdSubArgs[] = { pathToBin.get("Rscript"), rPath + R_SCRIPT_KINSHIP };
		OsCmdRunner cmdSub = new OsCmdRunner("SubBlock:" + sub, cmdSubArgs);
		cmdSub.getOsCmd().setQuiet(!debug, !debug); // Show stdout and stderr only in debug mode
		cmdSub.getOsCmd().setShowExceptions(debug); // Show exceptions only in debug mode
	}

	/**
	 * Check that all dependencies are correctly installed
	 */
	void checkDependencies() {
		if (verbose) Timer.showStdErr("Checking dependencies.");

		// Check R
		String optsR[] = { CMD_R, "--vanilla", "-e", "q()" };
		checkProgramPath(optsR, "You can install it from this link: http://www.r-project.org/.\nIf already installed, you can set the PATH to the program using '-pathBin' command line option.");

		// Check Rscript
		String optsRscript[] = { CMD_RSCRIPT, "-e", "q()" };
		checkProgramPath(optsRscript, "You can install it from this link: http://www.r-project.org/.\nIf already installed, you can set the PATH to the program using '-pathBin' command line option.");

		// Check fastLmm (this one doesn't work on OSX)
		String optsFastlmmc[] = { CMD_FASTLMM };
		checkProgramPath(optsFastlmmc, "You can install it from this link: http://fastlmm.codeplex.com/  or from here  http://research.microsoft.com/en-us/um/redmond/projects/MSCompBio/Fastlmm/ .\nIf already installed, you can set the PATH to the program using '-pathBin' command line option.");

		// Are all R libraries installed?
		for (String rlib : DEPENDECY_RLIBS)
			checkRLibrary(rlib);

		// Are all R scripts available?
		if (dependencyCheck) {
			for (String script : DEPENDECY_SCRIPTS) {
				String scriptPath = rPath + script;
				if (!Gpr.canRead(scriptPath)) missingDependency(scriptPath, "You can set the path to ASKAT's R scripts using option '-s' in the command line.");
			}
		}

		if (verbose && dependencyCheck) Timer.showStdErr("All dependencies found.\n");
	}

	/**
	 * Check if input files exists (TPED or TFAM input files)
	 */
	void checkOrCreateInputFile() {
		if (!Gpr.canRead(tpedFile)) {
			// No TPED file? Try to create one from a VCF file
			String vcfFile = genotypeName + ".vcf";
			if (Gpr.canRead(vcfFile)) vcf2Tped(vcfFile, tpedFile);
			else fatalError("Cannot read file '" + tpedFile + "'");
		}

		if (!Gpr.canRead(tfamFile)) fatalError("Cannot read file '" + tfamFile + "'");



    // Check that TPED and TFAM files have the same number of genotypes
		int tfamNumSamples = Gpr.countLines(tfamFile);
		int tpedNumSamples = (Gpr.countColumns(tpedFile) - 4) / 2;
		if (tfamNumSamples != tpedNumSamples) fatalError("Number of samples in TPED and TFAM files do not match:\n\t" + tfamNumSamples + "\tsamples in " + tfamFile + "\n\t" + tpedNumSamples + "\tsamples in " + tpedFile + "\n\tNote: TFAM samples are counted as number of lines. TPED samples are counted as (number_of_columns - 4)/2");

  }


	/**
	 * Check that a program is installed
	 * @param cmdArgs
	 * @param suggest
	 * 
	 * @return True if OK (i.e. runs the program an returns using exit status 0)
	 */
	boolean checkProgram(String cmdArgs[], String suggest) {
		String progName = cmdArgs[0];

		if (!dependencyCheck) {
			String path = findExecutableOnPath(progName); // Is it in PATH?
			return path != null;
		}

		if (verbose) Timer.showStdErr("Checking dependency: Program '" + progName + "'");
		OsCmdRunner cmd = new OsCmdRunner(progName, cmdArgs);
		cmd.getOsCmd().setQuiet(!debug, !debug); // Show stdout and stderr only in debug mode
		cmd.getOsCmd().setShowExceptions(debug); // Show exceptions only in debug mode
		if (debug) Timer.showStdErr("\tExecuting command: " + cmd.getOsCmd());
		cmd.run();
		boolean ok = (cmd.getExitValue() == 0) || !cmd.getHead().isEmpty() || !cmd.getHeadStderr().isEmpty(); // Exit status OK or any output to stdout or stderr? => OK
		if (verbose) Timer.showStdErr(ok ? "OK" : "Not found");
		return ok;
	}

	/**
	 * Check that a program is installed, using binPth
	 * @param cmdArgs
	 * @param suggest
	 */
	void checkProgramPath(String cmdArgs[], String suggest) {
		String progName = cmdArgs[0];

		if (checkProgram(cmdArgs, suggest)) {
			pathToBin.put(progName, progName);
			return; // Try using no path
		}

		for (String path : binPath.split(":")) {
			cmdArgs[0] = path + progName;
			if (checkProgram(cmdArgs, suggest)) {
				pathToBin.put(progName, cmdArgs[0]); // Store path to binary
				return; // OK
			}
		}

		// Error executing? Missing program
		missingDependency(progName, suggest);
	}

	/**
	 * Check that an R library is installed
	 * @param libraryName
	 */
	void checkRLibrary(String libraryName) {
		if (!dependencyCheck) return;

		if (verbose) Timer.showStdErr("Checking dependency: R library '" + libraryName + "'");
		String opts[] = { "R", "--vanilla", "-e", "library(" + libraryName + ")" };
		OsCmdRunner checkRlib = new OsCmdRunner("R_library_" + libraryName, opts);
		checkRlib.getOsCmd().setQuiet(!debug, !debug);
		if (debug) Timer.showStdErr("\tExecuting command: " + checkRlib.getOsCmd());
		checkRlib.run();
		if (checkRlib.getExitValue() != 0) missingDependency("R library " + libraryName, "You can install it by running the following command from R:\n\tinstall.packages(\"" + libraryName + "\")");
	}

	/**
	 * Count number of lines by chromosome
	 * @return
	 */
	CountByType countByChr() {
		if (verbose) Timer.showStdErr("Analyzing TPED file");

		LineFileIterator lfi = new LineFileIterator(tpedFile);
		CountByType countByChr = new CountByType();
		int posPrev = -1;
		String chrPrev = "";
		for (String line : lfi) {
			String fields[] = line.split("\\s", 4);
			String chr = fields[0];
			int pos = Gpr.parseIntSafe(fields[3]);
			countByChr.inc(chr);

			// Is the file sorted?
			if (chr.equals(chrPrev) && (pos < posPrev)) fatalError("TPED file out of order.\n\t" + line);

			chrPrev = chr;
			posPrev = pos;
		}

		if (verbose) Timer.showStdErr("Done:\n" + countByChr);
		return countByChr;
	}

	/**
	 * Fatal error
	 * @param msg
	 */
	void fatalError(String msg) {
		System.err.println("Error: " + msg);
		System.exit(1);
	}

	/**
	 * Find the full path to executable
	 * @param executableName
	 * @return
	 */
	private String findExecutableOnPath(String executableName) {
		// Is this a full path
		if (executableName.startsWith("/") || executableName.startsWith(".")) {
			File file = new File(executableName);
			if (file.isFile()) return executableName;
			else return null;
		}

		// Prepend PATH and test
		String systemPath = System.getenv("PATH");
		String[] pathDirs = systemPath.split(File.pathSeparator);

		File fullyQualifiedExecutable = null;
		for (String pathDir : pathDirs) {
			File file = new File(pathDir, executableName);
			if (file.isFile()) {
				fullyQualifiedExecutable = file;
				break;
			}
		}
		return fullyQualifiedExecutable == null ? null : fullyQualifiedExecutable.getAbsolutePath();
	}

	
	public Genome getGenome() {
		return genome;
	}

	public int getMinVariants() {
		return minVariants;
	}

	public int getNumWorkers() {
		return numWorkers;
	}

	public String getPath(String cmd) {
		return pathToBin.get(cmd);
	}

	public String getrPath() {
		return rPath;
	}

	public int getSubBlockSize() {
		return subBlockSize;
	}

	public String getTfamFile() {
		return tfamFile;
	}

	public boolean isDebug() {
		return debug;
	}

	public boolean isVerbose() {
		return verbose;
	}
	
	public double getpACC() {
		return pACC;
	}

	/**
	 * Show missing dependency error and exit
	 * @param dependency
	 * @param suggest
	 */
	void missingDependency(String dependency, String suggest) {
		if (!dependencyCheck) return;

		System.err.println("Missing dependency: " + dependency);
		if (suggest != null) System.err.println(suggest);
		System.exit(10);
	}

	@Override
	public void parseArgs(String[] args) {

		if (args.length == 0) usage(null);

		for (int i = 0; i < args.length; i++) {

			// Argument starts with '-'?
			if (args[i].startsWith("-")) {

				if (args[i].equals("-v") || args[i].equalsIgnoreCase("-verbose")) {
					verbose = true;
				} else if (args[i].equals("-h") || args[i].equalsIgnoreCase("-help")) {
					usage(null);
					System.exit(0);
				} else if (args[i].equals("-d")) {
					verbose = true;
					debug = true;
					KinshipBlock.debug = true;
				} else if (args[i].equals("-d1")) {
					verbose = true;
					debug = true;
					KinshipBlock.debug = true;
					KinshipBlock.debugOnlyOnce = true;
				} else if (args[i].equalsIgnoreCase("-noDep")) {
					dependencyCheck = false;
				} else if (args[i].equals("-i")) {
					blockSize = VERY_LARGE_BLOCK_SIZE; // We don't divide into blocks
					if ((i + 1) < args.length) bedFile = args[++i];
					else usage("Missing BED file.");
				} else if (args[i].equals("-p")) {
					if ((i + 1) < args.length) {
						numWorkers = Gpr.parseIntSafe(args[++i]);
						if (numWorkers <= 0) usage("Number of processes should be a positive number.");
					} else usage("Missing number of processes.");
				} else if (args[i].equalsIgnoreCase("-minvar")) {
					if ((i + 1) < args.length) {
						minVariants = Gpr.parseIntSafe(args[++i]);
						if (minVariants <= 0) usage("Minimum number of variants should be a positive number.");
					} else usage("Missing minimum number of variants.");
				} else if (args[i].equalsIgnoreCase("-pathR")) {
					if ((i + 1) < args.length) {
						rPath = args[++i];
						if (!rPath.endsWith("/")) rPath += "/"; // Make sure the path has a trailing slash
					} else usage("Missing path to ASKAT R scripts.");
				} else if (args[i].equalsIgnoreCase("-pathBin")) {
					if ((i + 1) < args.length) {
						binPath = args[++i];
						if (!binPath.endsWith("/")) binPath += "/"; // Make sure the path has a trailing slash
					} else usage("Missing path to ASKAT R scripts.");
				} else if (args[i].equals("-b")) {
					if ((i + 1) < args.length) {
						blockSize = Gpr.parseIntSafe(args[++i]);
						if (blockSize <= 0) usage("Block size must be a positive number.");
					} else usage("Missing block size.");
				} else if (args[i].equals("-sb")) {
					if ((i + 1) < args.length) {
						subBlockSize = Gpr.parseIntSafe(args[++i]);
						if (subBlockSize <= 0) usage("Sub-Block size must be a positive number.");
					} else usage("Missing sub-block size.");
				} else if (args[i].equals("-maxMaf")) {
					if ((i + 1) < args.length) {
						maxMaf = Gpr.parseDoubleSafe(args[++i]);
						if (maxMaf <= 0) usage("Maximum MAF should be a positive number.");
					} else usage("Missing sub-block size.");
				} else if (args[i].equals("-useMissingRef")) {
					useMissing = UseMissing.REFERENCE;
				} else if (args[i].equals("-useMissing")) {
					useMissing = UseMissing.MISSING;
				} else if (args[i].equals("-kin")) {
					if ((i + 1) < args.length) {
						String kinStr = args[++i].toLowerCase();
						if (kinStr.equals("all")) kinshipMethod = KinshipMethod.ALL;
						else if (kinStr.equals("chr")) kinshipMethod = KinshipMethod.CHROMOSOME;
						else if (kinStr.equals("avg")) kinshipMethod = KinshipMethod.CHROMOSOME_AVG;
						else if (kinStr.equals("block")) kinshipMethod = KinshipMethod.BLOCK;
					} else usage("Missing kinship type.");
				} else if (args[i].equalsIgnoreCase("-onlySnp")) {
					onlySnp = true;
				} else if (args[i].equalsIgnoreCase("-pACC")) { // UPD: add p-value accuracy option to improve with R-function "davies" numerical precision
					if ((i + 1) < args.length) {
						pACC = Gpr.parseDoubleSafe(args[++i]);
						if (pACC <= 0) usage("Accuracy must be a positive number.");
					} else usage("Missing accuracy value.");
				} else usage("Unknow option '" + args[i] + "'");
			} else if (genotypeName.isEmpty()) genotypeName = args[i];
			else usage("Unknow parameter '" + args[i] + "'");
		}

		// Sanity checks
		if (genotypeName.isEmpty()) usage("Missing genotypeName parameter");
		if ((blockSize < subBlockSize) || (blockSize % subBlockSize != 0)) usage("Block size (" + blockSize + ") must be a multiple of sub-block size (" + subBlockSize + ")");
	}

	@Override
	public boolean run() {
		// Salutation message
		if (verbose) {
			Timer.showStdErr("ASKAT algorithm by Karim Oualkacha");
/*			Timer.showStdErr(this.getClass().getSimpleName() + " wrapper version " + VERSION + "\n"); */
		}

		// Read intervals from BED file
		if (bedFile != null) {
			if (verbose) Timer.showStdErr("Loading intervals form file '" + bedFile + "'");
			BedFileIterator bed = new BedFileIterator(bedFile);
			intervals = bed.load(); // Read all intervals
			if (verbose) Timer.showStdErr("Total intervals : " + intervals.size());
			if (intervals.size() <= 0) throw new RuntimeException("No intervals found in file '" + bedFile + "'"); // Sanity check
		}

		// Initialize
		tpedFile = genotypeName + ".tped";
		tfamFile = genotypeName + ".tfam";

		// Check dependencies
		checkDependencies();

		// Create TPED file if it doesn't exist
		checkOrCreateInputFile();

		// Run algorithm 
		switch (kinshipMethod) {

		case CHROMOSOME_AVG:
			runChrAvg();
			break;

		case BLOCK:
		case CHROMOSOME:
		case ALL:
			runByBlock();
			break;

		default:
			throw new RuntimeException("Unimplemented algorithm for kinship method " + kinshipMethod);

		}

		return true;
	}

	/**
	 * Split main file into blocks and Run algorithm on each block
	 */
	void runByBlock() {
		int filtered = 0, remaining = 0;
		int countBlock = 0, pos = 0;
		String blockFileName = null;
		String chr = "", chrPrev = null, chrBlock = null;
		BufferedWriter blockFile = null;

		if (verbose) Timer.showStdErr("Creating blocks & Running algorithm on each block.");

		try {
			boolean forceLastBlockRun = false;

			// Iterate over input file
			LineFileIterator lfi = new LineFileIterator(tpedFile);
			for (String line : lfi) {
				// Parse
				TpedEntry tpedEntry = new TpedEntry(genome, line);
				double maf = tpedEntry.maf();

				// MAF within limits?
				if (maf <= maxMaf) {
					chr = tpedEntry.getChromosomeName();
					pos = tpedEntry.getStart();

					// Change of chromosome? 
					if ((blockFile == null) // Need to open file?
							|| ((kinshipMethod != KinshipMethod.ALL) && (chrPrev != null) && !chr.equals(chrPrev)) // Change of chromosome? (except for KinshipMethod.ALL)
							|| ((kinshipMethod == KinshipMethod.BLOCK) && (countBlock >= blockSize)) // Reached block size? (only for KinshipMethod.BLOCK)
					) {
						// Close file
						if (blockFile != null) {
							blockFile.close();
							if (verbose) Timer.showStdErr("Finished block file " + blockFileName + "'. Number of entries: " + countBlock);

							// Run commands
							runByBlock(blockFileName);
						}

						// Open a new file
						blockFileName = genotypeName + "." + "block." + chr + "_" + pos + ".tped";
						if (verbose) Timer.showStdErr("Creating block '" + blockFileName + "'");

						// If the file already exists, we can skip file creating process.
						// This only works for KinshipMethod.ALL method, because it only has only one TPED file
						if ((kinshipMethod == KinshipMethod.ALL) && Gpr.canRead(blockFileName)) {
							if (verbose) Timer.showStdErr("Block file " + blockFileName + "' already exists. Nothing done.");
							forceLastBlockRun = true;// This is set just to force the 'last block' processing
							break;

						}
						blockFile = new BufferedWriter(new FileWriter(new File(blockFileName)));

						// Prepare for next iteration
						chrBlock = null;
						countBlock = 0;
						chrPrev = chr;
					}

					if (chrBlock == null) chrBlock = chr;

					blockFile.write(line + "\n");
					countBlock++;
					remaining++;
				} else filtered++; // Filter out this line (not a rare variant)
			}

			// Last block (close file and run algorithm)
			if (forceLastBlockRun || ((blockFile != null) && (countBlock > 0))) {
				if (blockFile != null) blockFile.close();
				runByBlock(blockFileName); // Run commands
			}
		} catch (Exception e) {
			throw new RuntimeException(e);
		}

		if (verbose) Timer.showStdErr("Done. Filtered out (MAF) : " + filtered + " lines. Remaining: " + remaining + " lines.");
	}

	/**
	 * Run all commands for this block / sub-block combination
	 * @param blockFile
	 * @param subBlockFiles
	 */
	void runByBlock(String blockFile) {
		if (verbose) Timer.showStdErr("Running block: " + blockFile);

		KinshipBlock block = new KinshipBlock(this, blockFile);
		if (intervals != null) block.setIntervals(intervals);
		block.kinship();
		block.askat();
	}

	/**
	 * Run algorithm performing average kinship by chromosome
	 */
	void runChrAvg() {
		throw new RuntimeException("Unimplemented algorithm for kinship method " + kinshipMethod);
	}

	/**
	 * Return a list of sample names from a TFAM file
	 * @param tfamFile
	 * @return
	 */
	ArrayList<String> sampleNames(String tfamFile) {
		return null;
	}

	public void setBedFile(String bedFile) {
		this.bedFile = bedFile;
	}

	public void setBinPath(String binPath) {
		this.binPath = binPath;
	}

	public void setBlockSize(int blockSize) {
		this.blockSize = blockSize;
	}

	public void setDebug(boolean debug) {
		this.debug = debug;
	}

	public void setTfamFile(String tfamFile) {
		this.tfamFile = tfamFile;
	}

	public void setTpedFile(String tpedFile) {
		this.tpedFile = tpedFile;
	}

	public void setUseMissing(UseMissing useMissing) {
		this.useMissing = useMissing;
	}

	public void setVerbose(boolean verbose) {
		this.verbose = verbose;
	}

	/**
	 * Return REF anf ALT values as if they were a SNP
	 * 
	 * Important: If the variant is NOT a SNP, we create a 'fake' snp ( A -> T ).
	 * 			  This is done in order to be able to MAP InDels into PED files and keep compatibility with downstream programs (GenAble).
	 * 			  Yes, it's an awful hack. YOu've been warned!
	 */
	String snpGenotype(VcfEntry ve, VcfGenotype gen, int genoNum) {
		if (ve.isSnp()) {
			if (genoNum < 0) return ve.getRef(); // Reference
			return gen.getGenotype(genoNum);
		}

		// Create fake SNP "A -> T" and map InDel values to it
		if (genoNum < 0) return "A"; // Reference
		if (gen.getGenotype(genoNum).equals(ve.getRef())) return "A"; // ALT[genoNum] == REF 
		return "T"; // ALT[genoNum] != REF
	}

	@Override
	public void usage(String message) {
		if (message != null) System.err.println("Error: " + message + "\n");
		System.err.println("ASKAT algorithm by Karim Oualkacha");
		/* System.err.println(this.getClass().getSimpleName() + " wrapper version " + VERSION + "\n"); */
		System.err.println("Usage: java -jar " + this.getClass().getSimpleName() + ".jar [options] genotype");
		System.err.println("Options:");
		System.err.println("\t-b <num>       : Number of SNPs used for calculating the kinship matrix. Default: " + blockSize);
		System.err.println("\t-d             : Debug mode (implies verbose)");
		System.err.println("\t-d1            : Debug mode. Perform only one sub-block calculation and stop");
		System.err.println("\t-i <bed>       : BED file containing intervals to group SNPs. Default: none");
		System.err.println("\t-maxMaf        : Maximum MAF (minor allelel frequency). Default: " + maxMaf);
		System.err.println("\t-minVar num    : Minimum number of variants per group when using '-i' option. Default: " + minVariants);
		System.err.println("\t-noDep         : Do not perform dependency check.");
		System.err.println("\t-h             : Show this help and exit.");
		System.err.println("\t-kin <type>    : Kinship estimation type. Options {chr, avg, all, block}. Default: " + kinshipMethod);
		System.err.println("\t-onlySnp       : Use only SNPs when converting VCF to TPED. Default: " + onlySnp);
		System.err.println("\t-pathBin <dir> : Path to binary programs (e.g. FastLmm). Default: '" + binPath + "'.");
		System.err.println("\t-pathR <dir>   : Path to R scripts (ASKAT scripts). Default '" + rPath + "'.");
		System.err.println("\t-sb <num>      : Number of SNPs used for calculating the ASKAT algorithm. Default: " + subBlockSize);
		System.err.println("\t-useMissing    : Use entries with missing genotypes (otherwise they are filtered out). ");
		System.err.println("\t-useMissingRef : Use entries with missing genotypes marking them as 'reference' instead of 'missing'. ");
		System.err.println("\t-pACC <double> : Accuracy parameter for the p-value computation, default is 1e-9.");
		System.err.println("\t-v             : Be verbose.");
		System.exit(-1);
	}

	/**
	 * Convert a VCF to a TPED file
	 * @param vcfFile
	 * @param tpedFile
	 * 
	 * Important: If the variant is NOT a SNP, we create a 'fake' snp ( A -> T ).
	 * 			  This is done in order to be able to MAP InDels into PED files and keep compatibility with downstream programs (GenAble).
	 * 			  Yes, it's an awful hack. YOu've been warned!
	 * 
	 */
	public void vcf2Tped(String vcfFile, String tpedFile) {
		if (verbose) Timer.showStdErr("Converting file '" + vcfFile + "' to TPED format: '" + tpedFile + "'");

		int countVcf = 1, countTped = 0;
		int skipMissing = 0, skipNotSnp = 0, skipNonBiAllelic = 0;
		boolean useSample[] = null; // Which samples should be used
		try {
			// Open files
			VcfFileIterator vcf = new VcfFileIterator(vcfFile);
			BufferedWriter tped = new BufferedWriter(new FileWriter(tpedFile));

			// Convert VCF to TPED
			boolean isHeader = true;
			for (VcfEntry ve : vcf) {
				// Process header information
				if (isHeader) {
					useSample = vcfAndTfamSamples(vcf); // Consolidate TFAM and VCF samples
					isHeader = false;
				}

				// Warning: More than one ALT is not currently supported
				// Warning: Only SNPs are supported
				try {
					if (ve.getAlts().length != 1) { // No bi-allelic? => We skip it
						skipNonBiAllelic++;
						if (debug) System.err.println("Skipping line " + vcf.getLineNum() + ": Not bi-allelic");
					} else if (onlySnp && !ve.isSnp()) { // Not a SNP? skip it if 'onlySnp' is true
						skipNotSnp++;
						if (debug) System.err.println("Skipping line " + vcf.getLineNum() + ": Not a SNP");
					} else {
						boolean missingValues = false; // Any missing values in this line?

						// Prepare TPED line
						StringBuilder tpedLine = new StringBuilder();

						int pos = ve.getStart() + 1;
						String chr = ve.getChromosomeName();
						String id = "id_" + vcf.getLineNum(); // Create a unique ID
						if ((id == null) || id.isEmpty()) id = chr + ":" + pos;

						tpedLine.append(chr + " "); // Chromosome
						tpedLine.append(id + " "); // Identifier
						tpedLine.append("0 "); // Genetic distance in Moragans
						tpedLine.append(pos + " "); // Base pair position

						// Add all genotypes
						int i = 0;
						for (VcfGenotype gen : ve) {
							// Should we use this sample?
							if (useSample[i++]) {
								if (gen.getGenotypeCode() < 0) { // Missing genotype?
									missingValues = true;
									if (useMissing == UseMissing.REFERENCE) {
										String ref = snpGenotype(ve, gen, -1);
										tpedLine.append(ref + " " + ref + " "); // Mark both of them as reference
									} else tpedLine.append("0 0 "); // Mark both as missing
								} else {
									String gen0 = snpGenotype(ve, gen, 0);
									String gen1 = snpGenotype(ve, gen, 1);
									if (gen.getGenotype().length == 2) tpedLine.append(gen0 + " " + gen1 + " ");
									else {
										if (useMissing == UseMissing.REFERENCE) {
											String ref = ve.getRef();
											tpedLine.append(ref + " " + ref + " "); // Mark both of them as reference
										} else tpedLine.append("0 0 "); // Mark both as missing
									}
								}
							}
						}

						// Remove last space
						int lastChar = tpedLine.length() - 1;
						if (tpedLine.charAt(lastChar) == ' ') tpedLine.deleteCharAt(lastChar);

						tpedLine.append('\n');

						// Write to TPED file
						if ((useMissing != UseMissing.DO_NOT_USE) || !missingValues) {
							tped.write(tpedLine.toString());
							countTped++;
						} else {
							// Skipped because of misisng values?
							skipMissing++;
							if (debug) System.err.println("Skipping line " + vcf.getLineNum() + ": Missing values");
						}
					}

					countVcf++;
					if (verbose && (countVcf % 1000 == 0)) Timer.showStdErr("\tVCF to TPED:\tLine " + countVcf + "\t" + ve.getChromosomeName() + ":" + (ve.getStart() + 1));
				} catch (Exception e) {
					Gpr.debug("Exception processing VCF entry : " + ve);
					e.printStackTrace();
				}

			}

			// Close
			tped.close();
		} catch (IOException e) {
			throw new RuntimeException(e);
		}

		// Show some info
		if (verbose) Timer.showStdErr("Done: " + countVcf + " VCF entries converted to " + countTped + " TPED entries.\nSkipped entries:\n\tNon Biallelic: " + skipNonBiAllelic + "\n\tNon SNPs: " + skipNotSnp + "\n\tMissing genotypes: " + skipMissing);
	}

	/**
	 * Consolidate VCF and TFAM samples
	 * @param vcf
	 */
	boolean[] vcfAndTfamSamples(VcfFileIterator vcf) {
		// Open TFAM file
		Tfam tfam = new Tfam(tfamFile);

		// Get VCF samples
		List<String> sampleNamesVcf = vcf.getSampleNames();

		// Sanity check
		if (sampleNamesVcf.size() != tfam.size()) System.err.println("WARNING: Number of samples in TFAM file and VCF file do not match\n\tSamples in VCF file: " + sampleNamesVcf.size() + "\n\tSamples in TFAM file: " + tfam.size());

		// Create a 'common' list of samples
		HashSet<String> stfam = new HashSet<String>();
		stfam.addAll(tfam.getSampleNames());

		// Create a boolean array
		boolean use[] = new boolean[sampleNamesVcf.size()];
		int i = 0;
		for (String sampleNameVcf : vcf.getSampleNames())
			use[i++] = stfam.contains(sampleNameVcf);

		// Now we have to create a new TFAM file containing ONLY the samples in both VCF and TFAM files 
		// Note: The new file is sorted in the same order as the VCF.
		Tfam newTfam = new Tfam();
		for (String sampleNameVcf : vcf.getSampleNames()) {
			TfamEntry te = tfam.find(sampleNameVcf);
			if (te != null) newTfam.add(te);
		}

		// TFAM Files differ? Save new file
		if (!tfam.toString().equals(newTfam.toString())) {
			// Sanity check
			if (verbose) Timer.showStdErr("New TFAM file has " + newTfam.size() + " entries.");
			if (newTfam.size() <= 0) throw new RuntimeException("New TFAM file has no entries!");

			// We create a backup of the old file first
			String tfamFileBu = tfamFile + "." + (new Date()).getTime();
			Timer.showStdErr("Moving TFAM file from '" + tfamFile + "' to backup file '" + tfamFileBu + "'");
			(new File(tfamFile)).renameTo(new File(tfamFileBu));

			// Save new file
			Timer.showStdErr("Saving new TFAM file '" + tfamFile + "'");
			newTfam.save(tfamFile);
		}

		return use;
	}

	//UPD: Override for consistency between old and new version of CommandLine class
	@Override
	public String[] getArgs() {
		
		return null;
	}
}
