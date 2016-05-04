#POPEase

**Program Purpose**

Are you new to population genetics? Not familiar with command-line programs, or
just want a simple way to run basic pop gen analyses right on your labtop?
Then POPEase is perfect for you. The purpose of the POPEase program is to help
novices obtain reliable pop gen analyses for SNP data in just a few simple steps.

**What The Program Does**

POPEase (and the associated software) completes analyses that are
normally done using a series of GUI or command-line programs that can be complicated
for beginners.  These programs are integrated into the script and run with a few
simple command-line flags detailed in the program documents. The programs integrated
into this script are standard analyses when examining population genetics of your
data set.

This single python script bypasses the need to understand parsing your file into
the correct format for the downstream analyses, it bypasses FSAT and/or Genepop to
provide population genetics statistics, and it bypasses all programs associated with
population structure including Structure, CLUMPP or CLUMPPAK, distruct, and Structure
Harvester.

**Program Caveats**

First and foremost POPEase is written to be compatible with python version
=>3.0. It will not work with python 2.7 and you will need to update your python to the
current version. You can get the downloadable version of python [here](https://www.python.org/downloads/release/python-351/).

POPEase accepts a VCF file of DNA SNP data. The program is currently only
able to handle VCF data for diploid individuals and is defaulted to output
genotypes as missing if the read depth is not >5. EasyPOP also excludes
non-polymorphic SNPs. Phred Quality scores for VCF input are defaulted to only
be accepted if they are >20 (99% accurate base call) or better. However, if desired,
each of these caveats can be altered in the SPID_VCFtoSTRUCT file following the
instructions in the documentation.

POPEase is written to be easy to use on your personal laptop or
desktop. It can take a long time to run depending on the size of your data set
and may require the computer to be on and working for long periods of time.
Please be aware of that when starting the program. You cannot pause it while
Structure (the longest piece) is working.

##Dependencies
All dependencies can be downloaded via anaconda or from their websites listed
in the documentation.
1. Python 3.0
2. Numpy
3. MatPlotLib
4. Pandas
5. Structure software
6. PGDSpider Software
7. Java

##How to Get Started
1. It is highly recommended that you download and use Anaconda as your command line for this program, especially if you are not familiar with working on the command line.  The rest of the directions will be geared towards Anaconda users. Installing Anaconda will automatically install Python 3.0 for you if you do not have it on your system already. Download Anaconda here: https://www.continuum.io/downloads
2. Create a specific folder on your computer that will be used to deposit the required scripts and associated programs to run POPEase.
3. Download all files and scripts from the POPEase repository on GitHub and save them into the folder you created in step 1. THIS IS NECESSARY TO RUN THE PROGRAM CORRECTLY.
    a. Navigate to https://github.com/reige012/POPEase
    b. 	Click the “Download Zip” button on the right of the page
    c.	Save these files into the folder you created in step 2.
4. Download the following programs. Make sure to place the entire folder for each program into the new folder you created in step 2.
    a. PGDSpider- a file converting Software  http://www.cmpg.unibe.ch/software/PGDSpider/
    b.	STRUCTURE- a program to determine population structure from your data        http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/html/structure.html
    c.	JAVA: Mac Users: You may need to update or install Java on your system. http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html
5. Move individual files from each of the step 4 program folders into the main folder you created in step 2. The POPEase process is much smoother if there is one single folder that contains all of the programs and their associated files.
    a.	From PGDSpider: move all files, but not subfolders such as Example folder
    b.	From Console folder (Structure’s program folder) move three files:
        i.	 extraparams
        ii.  mainparams
        iii. structure (unix executable file)
6. Move a copy of your own VCF file to the folder you made in step 2.
    a.	Note: make sure your VCF file has a .vcf extension
7. That’s the basic setup. To run the programs please follow the detailed instructions in Making POPEase Work  in the documentation.
