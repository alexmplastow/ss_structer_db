import sys
sys.path.append('/usr/local/lib/python3.8/site-packages')
import RNA

from multiprocessing import Process
from multiprocessing import Manager

import os

from csv import writer
import csv

import time

import matplotlib.pyplot as plt

import numpy as np

#Input: dictionary of transcripts, transcript ID
#Output: List of strings of nucleotides, juxtaposed with secondary structure and MFE

def fast_fold(transcript_dict, transcript, RNAStructure_boolean):

    #You may want to renegotiate your use of the 'test nucleotide' variable,
        #it may be redundant

    #Defining the variables for the nucleic acid and two lists,
        #one which contains the nucleotides, structure, and free energies
        #and one which contains a lists of those lists
    #TODO: fix the 'test_nucleotide' misnomer
    test_nucleotide = transcript_dict[transcript]
    test_nucleotide_lists = []
    test_nucleotide_list = []

    #this variable roughly defines how many iterations the
        #folding algorithm will need to operate if we assume the nucleotide is
        #divisible by 40
    k = round(len(test_nucleotide) / 40)

    #Saving a few lines of code by writing a function to lengthen bracket
        #outputs
    def structure_extension(ss):
        while len(ss) < 120:
            ss+='.'
        return ss

    #This is configured to handle cases where the nucleotide length is less than 120
    #Even though you have the script for both calibration and folding to typically
    #handle nucleic acids greater than 120 bases
    if len(test_nucleotide) < 120 and not RNAStructure_boolean:
        (ss, mfe) = RNA.fold(test_nucleotide)
        ss = structure_extension(ss)
    elif len(test_nucleotide) < 120 and RNAStructure_boolean:
        nucleotides, ss, mfe = RNAStructure_fold(test_nucleotide)
        ss = structure_extension(ss)

    else:

        #The output of this set of for loops
        for x in range(0, k - 3):
            for ii in range(x * 40, 120 + x * 40):
                test_nucleotide_list.append(test_nucleotide[ii])
            nuke_structures=''.join(test_nucleotide_list)
            test_nucleotide_lists.append(nuke_structures)
            test_nucleotide_list = []
        for ii in range((k - 2) * 40, len(test_nucleotide)):
            test_nucleotide_list.append(test_nucleotide[ii])
        test_nucleotide_lists.append(''.join(test_nucleotide_list))

        #I reset the 'test_nucleotide_list' variable as an empty list because
            #I wanted to reuse the variable name
        test_nucleotide_list = []

        #Performing the folding and placing in the output
        if RNAStructure_boolean:
            for i in test_nucleotide_lists:
                nucleotides, ss, mfe = RNAStructure_fold(i)
                test_nucleotide_list.append(i + '.' + ss + ',' + str(mfe))
        else:
            for i in test_nucleotide_lists:
                (ss, mfe) = RNA.fold(i)
                ss = structure_extension(ss)
                test_nucleotide_list.append(i + ',' + ss + ',' + str(mfe))

        return test_nucleotide_list


#Input: transcript ID, dictionaries of structues, dictionary of transcripts
#Output: None, multiprocessing operates using global variables, so you simply add the information
    #to a dictionary
def slow_fold(transcript, structure_dict, transcript_dict):
    (ss, mfe) = RNA.fold(transcript_dict[transcript])
    structure_dict[transcript] = (ss, round(mfe, 2))

def txt_to_transcript_dict(file_name):
    txt = open(file_name)
    txt = txt.readlines()

    transcript_dict = {}

    for i, j in enumerate(txt):
        if '>' in j:
            key = j
            transcript_dict[key] = ''
        else:
            transcript_dict[key] += j.rstrip('\n')
    return transcript_dict

def calibration_folding(transcript_dict, calibration_arg, sample_size_arg, RNAStructure_boolean):
    transcript_list=list(set(transcript_dict))

    transcript_count = 0

    #Basic rationale of the process:
        #The for loop will start, spit out the transcript nucleotides
        # and get the folding algorithm running
        #If the molecule can be folded within the time period, a statement will be printed
        #to indicate the outcome, otherwise, the loop will finish up after 20 seconds just after
        #the process terminates and spit out a figure
    if calibration_arg:
        slow_fold_nucleotide_count=[]
        slow_fold_nucleotide_time=[]

    try:
        file=open('structure.csv', newline='')
        data=csv.reader(file)
        keys=[i[0] for i in data]
    except:
        keys=[]

    with open('structure.csv', 'a') as txt_object:
        writer_object = writer(txt_object)

        loop_structure_dict={}


        for transcript in transcript_list[0:sample_size_arg]:

            if transcript in keys:
                continue

            start_fold_time = time.perf_counter()

            #Outputs the complete nucleotide
            manager = Manager()
            structure_dict=manager.dict()

            action_process = Process(target=slow_fold, args=(transcript, structure_dict, transcript_dict))
            action_process.start()
            action_process.join(timeout=10)
            action_process.terminate()
            transcript_count += 1

            #Tells you how many transcripts have been spit out
            try:


                loop_structure_dict[transcript]=structure_dict.values()[0]
                #Outputs the dot brackets, not a whole lot else, you may want to add to the
                    #'structure dict' dictionary object
                fold_time=(time.perf_counter()-start_fold_time)
                writer_object.writerow([transcript, 'slow fold',transcript_dict[transcript],
                                        structure_dict.values()[0][0],
                                        structure_dict.values()[0][-1]])
                if calibration_arg:
                    slow_fold_nucleotide_count.append(len(transcript_dict[transcript]))
                    slow_fold_nucleotide_time.append(fold_time)

            except:
                nuke_list=fast_fold(transcript_dict, transcript, RNAStructure_boolean)
                #Contatins the nucleotide strings, their corresponding dot-brackets, and their minimum free energies
                    #all organized into a list
                writer_object.writerow([transcript] + ['fast fold'] + nuke_list)

    if calibration_arg:
        plt.plot(slow_fold_nucleotide_count, slow_fold_nucleotide_time, 'o')
        plt.savefig('test_plot.png')

        model = np.polyfit(slow_fold_nucleotide_count, slow_fold_nucleotide_time, 2)

        with open(calibration_arg, 'a') as txt_object:
            writer_object = writer(txt_object)
            writer_object.writerow(model)

def transcript_folding(sample_size, transcript_dict, nucleic_acid_length, RNAStructure_boolean):

    #Checking the CSV file to verify that a nucleic acid has been folded
    try:
        file = open('structure.csv', newline='')
        data = csv.reader(file)
        keys = [i[0] for i in data]
    except:
        keys = []

    with open('structure.csv', 'a') as txt_object:
        writer_object = writer(txt_object)

        transcripts = list(set(transcript_dict))[0:sample_size]

        for transcript in transcripts:
            if transcript in keys:
                continue

            if len(transcript_dict[transcript]) > nucleic_acid_length:


                #My nucleotide_structure_list is empty, that sucks
                nucleotide_structure_list = fast_fold(transcript_dict, transcript, RNAStructure_boolean)

                writer_object.writerow([transcript, 'fast fold'] + nucleotide_structure_list)

            elif len(transcript_dict[transcript]) <= nucleic_acid_length:

                #choosing between ViennaRNA and RNAStructure depending on user inputs
                if RNAStructure_boolean:
                    nucleotides, ss, mfe = RNAStructure_fold(transcript_dict[transcript])

                else:
                    (ss, mfe) = RNA.fold(transcript_dict[transcript])

                writer_object.writerow([transcript, 'slow fold', transcript_dict[transcript],
                                        ss, mfe])

#input: path to a file containing an RNAStructure dot-bracket prediction
#Ouput: tuple containing the nucleotides, dot-bracket structure, and lowest free energy value

def RNAStructure_from_file(file):
    #Opening the txt file and inserting all /n'ewlines into a list
    txt = open(file)
    txt = txt.readlines()

    #Eliminating the lingering nucleotides
    #I'm making the assumption that all nucleic acids have at least one unbound nucleotide
    bracket_energy_txt = [line for line in txt if '.' in line]
    nucleotides = tuple(set(bracket_energy_txt).symmetric_difference(set(txt)))
    nucleotides = nucleotides[0].rstrip('\n')
    txt = bracket_energy_txt

    #converting the txt list into a tuple of tuples containing information on the minimum free energy
    #and dot-bracket notation
    energy_structures = tuple((txt[index - 1].rstrip('\n'),line.rstrip('\n'))
                         for index,line in enumerate(txt)
                         if 'ENERGY' not in line and '.' in line)

    #Mapping all the tuples to a function so that the MFE quantity can be isolated and sorting by free energy
    MFE_f = lambda input: (input[-1],float(input[0].split(' ')[2]))
    energy_structures = list(map(MFE_f,energy_structures))
    energy_structures.sort(key = lambda x: x[1])

    #collecting the MFE structure, the variables may be misnomers and prepping them for return from the main function
    MFE_structure = energy_structures[-1]
    structure = MFE_structure[0]
    MFE = MFE_structure[-1]

    return nucleotides, structure, MFE

def RNAStructure_fold(nucleotide_sequence):
    #Opening a file for inserting nucleic acid sequence and inserting
    txt = open('RNAStructure_nucleic_acid.txt', "w")
    txt.write(nucleotide_sequence)
    txt = open('RNAStructure_nucleic_acid.txt')

    #Executing bash commands to
    #this doesn't seem to be executing as intended

    os.system('./RNA_structure.sh')

    #Collecting nucleotides, dot bracket structure, and free energy
    nucleotides, structure, MFE = RNAStructure_from_file(
        'RNAStructure_bracket_output.txt')

    return nucleotides, structure, MFE
