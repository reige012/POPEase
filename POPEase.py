#!/usr/bin/env python
# encoding: utf-8

"""


Edited by Alicia Reigel. 15 April 2016.
Copyright Alicia Reigel. Louisiana State University. 15 April 2016. All
rights reserved.

"""

import csv
import re
import glob
import os
import sys
import argparse
import subprocess
import numpy
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages


"""class FullPaths(argparse.Action):
    Expand user- and relative-paths. Obtained from Brant Faircloth.
    https://github.com/biolprogramming/assignment-17/blob/master/answers/brantfaircloth/task1.py
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath((values)))"""


def parser_get_args():
    """Collect the path to the SNP data (forward and reverse), vcf data, and
       base name for output files"""
    parser = argparse.ArgumentParser(
        description="""Input the path to species name file and desired
            output file name"""
        )
    parser.add_argument(
            '--directorypath',
            required=True,
            type=str,
            help='Enter the full path to the folder/directory where the required programs are located.'
        )
    parser.add_argument(
            '--vcfpath',
            required=True,
            type=str,
            help='Enter the full path to the database.'
        )
    parser.add_argument(
            '--outputfolder',
            required=True,
            type=str,
            help='Enter the desired name for the output folder.'
        )
    parser.add_argument(
            '--structurepops',
            required=True,
            type=str,
            help='Enter the maximum populations (K) to test when running structure.'
        )
    parser.add_argument(
            '--loci',
            required=True,
            type=str,
            help='Enter the number of loci in your data set.'
        )
    parser.add_argument(
            '--individuals',
            required=True,
            type=str,
            help='Enter the number of individuals that you will include in your structure run.'
        )
    parser.add_argument(
            '--numruns',
            type=int,
            default=10,
            help='Enter the number of  structure runs for each K value.'
        )
    return parser.parse_args()


def pgd_structure(outputfolder, inputfile, directorypath):
    """Converts the user's VCF file into a STRUCTURE-formatted file"""
    mypath = os.path.join(directorypath + "/STRUCTURE_directory")
    if not os.path.isdir(mypath):
        os.makedirs(mypath)
    path_to_new_directory = os.path.abspath(mypath)
    # gets the absolute path to the new directory created
    PGDStructstdout = os.path.join(path_to_new_directory, "PGDStructstdout.txt")
    PGDStructstderr = os.path.join(path_to_new_directory, "PGDStructstderr.txt" )
    outputfile_pgdstruct = os.path.join(path_to_new_directory, "STRUCTUREformatfile")
    outputfileabspath = os.path.abspath(outputfile_pgdstruct)
    STRUC_PGD_params = ['PGDSpider2-cli.exe', '-inputfile', inputfile,
        '-inputformat', 'VCF', '-outputfile', outputfile_pgdstruct, '-outputformat',
        'STRUCTURE', '-spid', 'SPID_VCFtoSTRUCT.spid']
    with open(PGDStructstdout, 'w') as stdout_file:
        with open(PGDStructstderr, 'w') as stderr_file:
            my_proc1 = subprocess.run(STRUC_PGD_params,
            stdout=stdout_file,
            stderr=stderr_file,
            universal_newlines=True)
    print("\nCreated output file in Structure format: {}".format(outputfile_pgdstruct))
    return path_to_new_directory, outputfileabspath


def structure_run(structK, loci, individuals, numruns, pgd_structformatted_file, outputfile_struct_run, directorypath):
    """Runs STRUCTURE based on user inputs and mainparams file defaults"""
    K = (int(structK) + 1)
    numruns = (int(numruns) + 1)
    for x in range(1, K):
        for i in range(1, numruns):
            x = str(x)
            i = str(i)
            outputfile = os.path.join(outputfile_struct_run + x + '_' + i)
            list_struct_run = ['structure', '-K', x, '-o', outputfile, '-i', pgd_structformatted_file, '-L', loci, '-N', individuals]
            stdoutfile = os.path.join(directorypath, 'STRUCTUREstdout_' + x + '_' + i)
            stderrfile = os.path.join(directorypath, 'STRUCTUREsterr_' + x + '_' + i)
            with open(stdoutfile, 'w') as struct_stdout_file:
                with open(stderrfile, 'w') as struct_stderr_file:
                    my_proc_2 = subprocess.run(list_struct_run,
                    input=pgd_structformatted_file,
                    stdout=struct_stdout_file,
                    stderr=struct_stderr_file,
                    universal_newlines=True)
    print('''\nFiles created for this process are standard out files,
    standard error files and structure run data files.  Each file is identified
    by the number of populations (K) being tested in that run, followed by an
    underscore and the number of the run.\n
    The ln values will now be compiled and the best estimate for K will be
    computed. Please wait.''')


def compute_ln(structK, directorypath):
    os.chdir(directorypath)
    with open('Best_K_Analysis.csv', 'w') as LN_output:
        headers = ['K', 'Mean_LnP(D)', 'StDevLN']
        writer = csv.writer(LN_output)
        writer.writerow(headers)
        K = (int(structK) + 1)
        for x in range(1, K, 1):
            x = str(x)
            filename = os.path.join("*STRUCTURERUN" + x + '*')
            path_name = os.path.join(directorypath, filename)
            # finds the pathnames for any files matching the filename
            file_path_list = glob.glob(path_name)
            # creates a list of the path names associated with the files found
            LN_prob_list = []
            for file in file_path_list:
                with open(file, 'r') as info:
                    for line in info:
                        if re.search("Estimated Ln Prob of Data", line):
                            new_list = line.replace('\n', '').split('   = ')
                            LN_prob_list.append(float(new_list[1]))
            meanLN = (sum(LN_prob_list) / float(len(LN_prob_list)))
            stdevln = numpy.std(LN_prob_list)
            row = [x, meanLN, stdevln]
            writer.writerow(row)
            path = os.path.abspath('Best_K_Analysis.csv')
    return path


def calculate_ln1P(directorypath, path):
    os.chdir(directorypath)
    dataframe = pd.read_csv(path)
    list1 = (-(dataframe.loc[:, 'Mean_LnP(D)'] - dataframe.loc[:, 'Mean_LnP(D)'].shift(-1)))
    list2 = list(list1)
    list2.insert(0, 'NaN')
    list2.pop()
    dataframe["Ln'P(D)"] = pd.DataFrame(list2)
    dataframe['Ln"P(D)'] = (dataframe.loc[1:, "Ln'P(D)"] - dataframe.loc[:, "Ln'P(D)"].shift(-1))
    dataframe['Ln"P(D)'] = dataframe['Ln"P(D)'].abs()
    dataframe['Delta_K'] = (dataframe.loc[1:, 'Ln"P(D)'] / dataframe.loc[:, 'StDevLN'])
    best_ln = dataframe['Mean_LnP(D)'].min()
    best_deltak = dataframe['Delta_K'].max()
    ln_row = dataframe.loc[(dataframe['Mean_LnP(D)'] == best_ln), 'K': 'Mean_LnP(D)']
    deltak_row = dataframe.loc[(dataframe['Delta_K'] == best_deltak), 'K': 'Delta_K']
    with PdfPages('Best_K_Figures.pdf') as pdf:
        dataframe.plot(x='K', y='Mean_LnP(D)')
        """plt.ylim(plt.ylim()[::-1])"""
        # flips y-axis to put lowest value on top. If desired ten simply remove three tick marks on each side of phrase.
        plt.savefig("Mean_LnP(D)_Figure.png", bbox_inches='tight')
        plt.title('Mean_LnP(D)')
        pdf.savefig()  # saves the current figure into a pdf page
        plt.close()

        dataframe.plot(x='K', y='Delta_K')
        plt.savefig("DeltaK_figure.png", bbox_inches='tight')
        plt.title('Delta K')
        pdf.savefig()
        plt.close()
    print("\nThe data for best K value based on Mean_LnP(D) is below:")
    print("{}\n".format(ln_row))
    print("The data from best K value based on Delta K is below:")
    print("{}\n".format(deltak_row))
    print('''These "best" K values are only suggestions based on the lowest
    value for Mean_LnP(D) and the highest value for Delta K. You should always
    check the Structure documentation and consider your data and study system
    carefully before choosing a final K value.\n''')
    with open(path, 'w') as csvfinal:
        dataframe.to_csv(csvfinal)
    print('''The output files include a PDF (Best_K_Figures.pdf) containing
    graphs for Mean_LnP(D) and Delta K as well as a .png figure file for each
    graph. Additionally, a CSV file (Best_K_Analysis.csv) of all related data
    was created.\nThis program is finished running. To obtain a graphical
    display of the Structure results you should see the documentation.''')


"""def pgd_fstat():
    #dkfneaifn """


def main():
    args = parser_get_args()
    inputfile = os.path.abspath(args.vcfpath)
    directorypath = os.path.abspath(args.directorypath)
    numruns = args.numruns
    path_to_new_directory, outputfileabspath = pgd_structure(args.outputfolder, inputfile, directorypath)
    print('''\nPlease check your STRUCTURE-formatted file for missing data
    Missing or bad quality data is indicated by a -9. You may want to remove
    individuals with high amounts of missing data. Please follow instructions
    in the documentation to do so and re-run this program. If your data is
    correct and ready for running through STRUCTURE please enter "Y". If not,
    please enter "N", correct your VCF file and rerun this program.''')
    y = ((str(input())).lower())
    if y != "y" or "n":
        print("Please enter either Y or N.")
        y = (str(input())).lower()
        if y == "y":
            print('''Great, STRUCTURE will run now. This may take a significant
    amount of time depending on the size of your data set. Please do not
    close this window or type. If the curser is blinking the program
    is running.''')
            outputfile_struct_run = os.path.join(path_to_new_directory, args.outputfolder + '_STRUCTURERUN')
            structK = args.structurepops
            loci = args.loci
            individuals = args.individuals
            structure_run(structK, loci, individuals, numruns, outputfileabspath, outputfile_struct_run, path_to_new_directory)
            path = compute_ln(structK, path_to_new_directory)
            calculate_ln1P(path_to_new_directory, path)
        if y == "n":
            print('Exiting program.')
            sys.exit()


if __name__ == '__main__':
    main()
