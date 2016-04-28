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
be accepted if they are >20 (99% accurate base call) or better.

Lastly, POPEase is written to be easy to use on your personal laptop or
desktop. It can take a long time to run depending on the size of your data set
and may require the computer to be on and working for long periods of time.
Please be aware of that when starting the program. You cannot pause it while
Structure (the longest piece) is working.

##Instructions for Use
1. Create a specific folder on your computer to deposit the required scripts and
   associated programs to run POPEase
2. Download all files and scripts from the POPEase repository on GitHub and save
   them into the folder you created in step 1. THIS IS NECESSARY TO RUN THE
   PROGRAM CORRECTLY.
3. Download the following programs and save them into the folder you created:
    - PGDSpider- a file converting Software
        http://www.cmpg.unibe.ch/software/PGDSpider/
    - STRUCTURE- a program to determine population structure from your data
        http://pritchardlab.stanford.edu/structure_software/release_versions/v2.3.4/html/structure.html
    - FSTAT- a program to provide basic stats about your SNP data
4. Move all items from downloaded folders into the main folder you created. The
   process is much smoother if there is one single file that contains all of the
   programs.
5. Move a copy of your VCF file (with a .vcf extension) to the created folder
6. Follow the specific input command line steps (with advanced options) found
   in the documentation for POPEase. 
