
												ASKAT
												-----
		
-----------------------------------------------------------------------------------------------------------------------

To do
-----
		
All the user should ever type is

  askat genotype

The command should:

  1- Make sure there are two files "genotype.tped" and "genotype.tfam" (or even a genotype.vcf)
  2- Make sure that all dependencies are installed
  
  3- Split input file in "chunks" of 2000 SNPs
  		For each chunk:
  			3.1- Calculate the kinship matrix
  			3.2- Sub-divide into parts cosisting of 10 to 50 SNPs
  			3.3- Create a file in 'askat' format (see example 'Ped_EX_ASKAT_NOMissdata.dat')
  			3.4- Call call ASKAT R function (in parallel) using that part and the pre-calculated kinship matrix  			
  			3.5- Collect the results
  			
  4- Show the results in a nice summary

-----------------------------------------------------------------------------------------------------------------------
												
Intall
------

	# Create a directory for ASKAT
	
		mkdir askat
		cd askat
		

	# Install R packages
	
		# Run R
		R
		
		# Run install command (within R)
		install.packages("GenABEL")
		install.packages("CompQuadForm")
		install.packages("nFactors")
		install.packages("MASS")
		
		q()
			
	# Install FastLmm: http://fastlmm.codeplex.com/
	
		IMPORTANT: Version 1.09 seems to be broken. You should use version 1.02

		# Download and uncompress (if the link doesn't work, go to http://fastlmm.codeplex.com/ and click on "Download") 	
		wget -O FastLmm.linux.zip "http://download.codeplex.com/Download?ProjectName=fastlmm&DownloadId=365784&FileTime=129784733225770000&Build=18696"
		
		# Uncompress the package
		unzip FastLmm.linux.zip
		
		# Copy the binary to current directory and make it executable
		cp FaSTLMM.109.Linux/Linux_MKL/fastlmmc .
		chmod a+x fastlmmc
		
		
	# Run demo
	
		$ R
		
		library(nFactors)
		library(CompQuadForm)

		path = "./r"
		path = "."
		source(file = paste(path, "/ASKAT.R", sep=""))
		source(file = paste(path, "/Geno.FaST.LMM.R", sep=""))
		source(file = paste(path, "/Get_Lambda.R", sep=""))
		source(file = paste(path, "/Get_PValue.Modif.R", sep=""))
		source(file = paste(path, "/pheno.geno.kin.Without.NA.R", sep=""))
		source(file = paste(path, "/VC.FaST.LMM.R", sep=""))

		load(paste(path, "/kin1.Rdata", sep=""))
		Ped = read.table(file = paste(path, "/Ped_EX_ASKAT_NOMissdata.dat", sep=""), header=FALSE)
		Missing = FALSE
		ASKAT(Ped, kin1, Missing)
		
-----------------------------------------------------------------------------------------------------------------------

R-Java integration
------------------

References
	http://stackoverflow.com/questions/7451716/java-r-integration





Options
	rcaller	http://code.google.com/p/rcaller/
			http://stdioe.blogspot.ca/2011/07/rcaller-20-calling-r-from-java.html

			Dependencies:
				install.packages("Runiversal")			Installs OK

	JRI:	http://www.rforge.net/JRI/
			Now JRI is part of rJava

			Doesn't work:
				install.packages("rJava")
				ERROR: configuration failed for package ‘rJava’
				installation of package ‘rJava’ had non-zero exit status

	Rserve	Rserve is a TCP/IP server which allows other programs to use facilities of R (see www.r-project.org) from various languages without the need to initialize R or link against R library

			Problems: Doesn't look multi-threading

	REngine


	Rcli

	renjin	http://code.google.com/p/renjin/
			Renjin is a new JVM-based interpreter for the R language. It is intended to be 100% compatible with the original interpreter
