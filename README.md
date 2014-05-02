Askat 
=====

You can <b>download</b> the program <a href="http://sourceforge.net/projects/snpeff/files/askat.zip"> <b><font size=+1>here</font></b> </a><p> or you can get the whole project and a zip file using this command: git clone https://github.com/pcingola/Askat.git (requires git installed).
You can take a look at some <a href="html/AskatWrapper.pdf">slides</a> on the wrapper.

<hr><center><h3>Basics</h3></center>

<b>Requirements</b><p>

<ul>
  <li> R : Make sure you have a resonably up to date version of R installed and availabe in you PATH.
	<li> Java 1.6 : Most modern computers have Java installed.
	<li> Fast-LMM : Pre-compiled versions of fast-Lmm can be found <a href="http://fastlmm.codeplex.com/"> here </a> (Linux and Windows).<br>
			<b>WARNING:</b> You should use relatively new version of Fast-LMM (e.g. version v1.09).
	<li> R Libraries: The ASKAT wrapper will tell you how to to install any missing R libraries.
</ul>

<b>Install</b><p>

Just unzip the file and run the JAR by doing:
<pre>
unzip askat.zip
</pre>

<b>Running the program</b><p>

In order to run the program you need to go into the directory where you installed it an run the JAR file:

<pre>
cd askat
java -Xmx2G -jar Askat.jar [options] genotype
</pre>

E.g.: To show a help message, just run the program without any parameters:

<pre>
$ java -jar Askat.jar
ASKAT algorithm by Karim Oualkacha, optimized from N^3 to N^2 complexity by Stepan Grinek 
Askat wrapper version 1.01 (build 2012-11-21), by Pablo Cingolani

Usage: java -jar Askat.jar [options] genotype
Options:
	-b &lt;num&gt;       : Number of SNPs used for calculating the kinship matrix. Default: 100000
	-d             : Debug mode (implies verbose)
	-d1            : Debug mode. Perform only one sub-block calculation and stop
	-i &lt;bed&gt;       : BED file containing intervals to group SNPs. Default: none
	-maxMaf        : Maximum MAF (minor allelel frequency). Default: 1.0
	-noDep         : Do not perform dependency check.
	-h             : Show this help and exit.
	-kin &lt;type&gt;    : Kinship estimation type. Options {chr, avg, all, block}. Default: CHROMOSOME
	-p &lt;num&gt;       : Number of parallel processes. Default: 8
	-pathBin &lt;dir&gt; : Path to binary programs (e.g. FastLmm). Default: './'.
	-pathR &lt;dir&gt;   : Path to R scripts (ASKAT scripts). Default './r/'.
	-sb &lt;num&gt;      : Number of SNPs used for calculating the ASKAT algorithm. Default: 20
	-pACC &lt;double&gt; : Accuracy parameter for the p-value computation, default is 1e-9.
  -v             : Be verbose.
</pre>

<b>Some comments about the 'intervals' option</b><p>

Whenever command line option '-i file.bed' is specified, Askat wrapper analyses all variants matching those intervals as a single block.
This means that variants are not subdivided into sub-blocks. 
This is effectively like setting the subBlock parameter ('-sb') to infinity.<br>

If a variant hits more than one interval, then it is analyzed on all intervals.
Although overlaping intervals are allowed, each interval is 'treated' as a unique block (this means that no statistical corrections are made in the model).<br>
<p>

<b>Running a full example</b><p>

The genotype 'karim4k' is available for testing the program:

<p>
<b>WARNING!</b> By default Askat will try to use ALL the processors available in the computer. This will produce a significant slow down for all other users in the system.
<p>

The command line fo the example is:
<pre>
$ <b>java -Xmx1g -jar Askat.jar -v karim4k | tee karim4k.out </b>
</pre>
Note that we specified 1G of ram to be used by the program (parameter '-Xmx1G'), but in some other cases, more memory me be required.
<p>

The output looks like this:
<pre>
<font size=-1>
00:00:00.001  ASKAT algorithm by Karim Oualkacha
00:00:00.004  Askat wrapper version 1.0 - epsilon 'almost 1.0' (build 2012-06-29), by Pablo Cingolani

00:00:00.004  Checking dependencies.
00:00:00.005  Checking dependency: Program 'R'
00:00:00.274  OK
00:00:00.275  Checking dependency: Program 'Rscript'
00:00:00.434  OK
00:00:00.434  Checking dependency: Program 'fastlmmc'
00:00:00.535  OK
00:00:00.536  Checking dependency: R library 'GenABEL'
00:00:02.099  Checking dependency: R library 'CompQuadForm'
00:00:02.387  Checking dependency: R library 'nFactors'
00:00:02.982  Checking dependency: R library 'MASS'
00:00:03.286  All dependencies found.

00:00:03.287  Running algorithm.
00:00:03.287  Creating blocks.
00:00:04.423  Creating block 'karim4k.block.1_100.tped'. Number of entries: 4000
00:00:04.511  Running block: karim4k.block.1_100.tped
00:00:04.512  Calculating kinship matrix for block: karim4k.block.1_100
00:00:23.197  Starting block: karim4k.block.1_100
00:00:23.284  Create batches.
      File 'karim4k.block.1_100.tped' has 4000 lines.
      Split up to 500 lines per batch.
00:00:23.287  Batch 1. Line 1. Creating batch : karim4k.block.1_100.1.askat
00:00:23.496  Batch 2. Line 501. Creating batch : karim4k.block.1_100.2.askat
00:00:23.673  Batch 3. Line 1001. Creating batch : karim4k.block.1_100.3.askat
00:00:23.845  Batch 4. Line 1501. Creating batch : karim4k.block.1_100.4.askat
00:00:24.019  Batch 5. Line 2001. Creating batch : karim4k.block.1_100.5.askat
00:00:24.190  Batch 6. Line 2501. Creating batch : karim4k.block.1_100.6.askat
00:00:24.365  Batch 7. Line 3001. Creating batch : karim4k.block.1_100.7.askat
00:00:24.538  Batch 8. Line 3501. Creating batch : karim4k.block.1_100.8.askat
00:01:33.747  ASKAT_RESUTS:  Block:  karim4k.block.1_100.1.askat  Sub-Block:  1 - 20   chr:pos:  1:100 - 1:2000       Id:  snp_1 - snp_20       p-value:  0.0006166988  Q:  38588.06  ...
00:01:35.749  ASKAT_RESUTS:  Block:  karim4k.block.1_100.2.askat  Sub-Block:  1 - 20   chr:pos:  1:50100 - 1:52000    Id:  snp_501 - snp_520    p-value:  0.6740977     Q:  10817.98  ...
00:01:36.750  ASKAT_RESUTS:  Block:  karim4k.block.1_100.3.askat  Sub-Block:  1 - 20   chr:pos:  1:100100 - 1:102000  Id:  snp_1001 - snp_1020  p-value:  4.551971e-07  Q:  57434.71  ...
00:01:36.750  ASKAT_RESUTS:  Block:  karim4k.block.1_100.8.askat  Sub-Block:  1 - 20   chr:pos:  1:350100 - 1:352000  Id:  snp_3501 - snp_3520  p-value:  0.6530349     Q:  11460.14  ...
00:01:37.751  ASKAT_RESUTS:  Block:  karim4k.block.1_100.4.askat  Sub-Block:  1 - 20   chr:pos:  1:150100 - 1:152000  Id:  snp_1501 - snp_1520  p-value:  0.9720287     Q:  5038.441  ...
00:01:37.751  ASKAT_RESUTS:  Block:  karim4k.block.1_100.6.askat  Sub-Block:  1 - 20   chr:pos:  1:250100 - 1:252000  Id:  snp_2501 - snp_2520  p-value:  0.1502854     Q:  18876.09  ...
00:01:38.752  ASKAT_RESUTS:  Block:  karim4k.block.1_100.5.askat  Sub-Block:  1 - 20   chr:pos:  1:200100 - 1:202000  Id:  snp_2001 - snp_2020  p-value:  1.340421e-06  Q:  54710.54  ...
00:01:38.753  ASKAT_RESUTS:  Block:  karim4k.block.1_100.7.askat  Sub-Block:  1 - 20   chr:pos:  1:300100 - 1:302000  Id:  snp_3001 - snp_3020  p-value:  1.385112e-05  Q:  48982.58  ...
00:02:22.786  ASKAT_RESUTS:  Block:  karim4k.block.1_100.1.askat  Sub-Block:  21 - 40  chr:pos:  1:2100 - 1:4000      Id:  snp_21 - snp_40      p-value:  0.2958578     Q:  16597.14  ...
00:02:25.788  ASKAT_RESUTS:  Block:  karim4k.block.1_100.2.askat  Sub-Block:  21 - 40  chr:pos:  1:52100 - 1:54000    Id:  snp_521 - snp_540    p-value:  0.7232833     Q:  10029.99  ...
00:02:26.789  ASKAT_RESUTS:  Block:  karim4k.block.1_100.3.askat  Sub-Block:  21 - 40  chr:pos:  1:102100 - 1:104000  Id:  snp_1021 - snp_1040  p-value:  0.01380279    Q:  28006.57  ...
00:02:27.790  ASKAT_RESUTS:  Block:  karim4k.block.1_100.8.askat  Sub-Block:  21 - 40  chr:pos:  1:352100 - 1:354000  Id:  snp_3521 - snp_3540  p-value:  0.001698173   Q:  37035.88  ...
00:02:28.791  ASKAT_RESUTS:  Block:  karim4k.block.1_100.4.askat  Sub-Block:  21 - 40  chr:pos:  1:152100 - 1:154000  Id:  snp_1521 - snp_1540  p-value:  0.1902052     Q:  18511.61  ...
00:02:28.791  ASKAT_RESUTS:  Block:  karim4k.block.1_100.5.askat  Sub-Block:  21 - 40  chr:pos:  1:202100 - 1:204000  Id:  snp_2021 - snp_2040  p-value:  0.6349479     Q:  10588.31  ...
00:02:28.792  ASKAT_RESUTS:  Block:  karim4k.block.1_100.6.askat  Sub-Block:  21 - 40  chr:pos:  1:252100 - 1:254000  Id:  snp_2521 - snp_2540  p-value:  0.8626755     Q:  8638.986  ...
00:02:28.792  ASKAT_RESUTS:  Block:  karim4k.block.1_100.7.askat  Sub-Block:  21 - 40  chr:pos:  1:302100 - 1:304000  Id:  snp_3021 - snp_3040  p-value:  0.3494253     Q:  15464.29  ...
...
...
...
</font>
</pre>

<b>Input data formats</b><p>
ASKAT wrapper requires the input data to be formatted in two file: TFAM and TPED format (for details, see <a href="http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml">PLINK</a> software package).<p>
<ul>
	<li> TFAM: This is the typical "Transposed FAM". TFAM file have individual and family information, where one row is an individual. <br>
	The columns are:
	<ul>
		<li> Family ID
		<li> Individual ID
		<li> Paternal ID
		<li> Maternal ID
		<li> Sex (1=male; 2=female; other=unknown)
		<li> Phenotype
	</ul> 
	<li> TPED: This is a "Transposed PED". TPED file have SNP information, where one row is a SNP (for all samples). <br>
	The columns are:
	<ul>
		<li> chromosome (1-22, X, Y or 0 if unplaced)
		<li> rs# or snp identifier
		<li> Genetic distance (morgans)
		<li> Base-pair position (bp units)
		<li> Columns 5 and on: Genotype informaton (two bases per sample, assuming )
	</ul> 
</ul>

E.g.: TFAM data:
<pre>
1 1_1 0 0 1 4.91995
1 1_2 0 0 2 7.14442
1 1_3 1_1 1_2 2 5.26482
2 2_1 0 0 1 2.87951
2 2_2 0 0 2 1.52721
2 2_3 2_1 2_2 2 2.45878
2 2_4 2_1 2_2 1 2.59495
3 3_1 0 0 1 2.17923
3 3_2 0 0 2 3.91823
3 3_3 3_1 3_2 1 2.55475
</pre>

<p>
E.g.: TPED data (long lines have been truncated):
<pre>
1 snp_1 0 100 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_2 0 200 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_3 0 300 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_4 0 400 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_5 0 500 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_6 0 600 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A T T A A T A T A T A T A A A T A A A A A T A T A A A ...
1 snp_7 0 700 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_8 0 800 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A ...
1 snp_9 0 900 A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A A T A A A A A T A T A A A A A ...
</pre>