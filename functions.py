from multiprocessing import Process
from multiprocessing import Manager

from tqdm import tqdm

import os

from csv import writer
import csv

import time

import matplotlib.pyplot as plt

import numpy as np


# Determinines the algorithm type based on user input
algorithm_type_test = lambda RNAStructure_boolean: \
    'RNAStructure' if RNAStructure_boolean else 'ViennaRNA'

#Input: dictionary of transcripts, transcript ID
#Output: List of strings of nucleotides, juxtaposed with secondary structure and MFE

def fast_fold(transcript_dict, transcript, RNAStructure_boolean):

    #You may want to renegotiate your use of the 'test nucleotide' variable,
        #it may be redundant

    #Defining the variables for the nucleic acid and two lists,
        #one which contains the nucleotides, structure, and free energies
        #and one which contains a lists of those lists
    nucleic_acid = transcript_dict[transcript]
    nucleic_acid_lists = []
    nucleic_acid_list = []

    #this variable roughly defines how many iterations the
        #folding algorithm will need to operate if we assume the nucleotide is
        #divisible by 40
    k = round(len(nucleic_acid) / 40)

    #Saving a few lines of code by writing a function to lengthen bracket
        #outputs
    def structure_extension(ss):
        while len(ss) < 120:
            ss+='.'
        return ss

    #This is configured to handle cases where the nucleotide length is less than 120
    #Even though you have the script for both calibration and folding to typically
    #handle nucleic acids greater than 120 bases
    if len(nucleic_acid) < 120 and not RNAStructure_boolean:
        (ss, mfe) = ViennaRNA_fold(nucleic_acid)
        ss = structure_extension(ss)
    elif len(nucleic_acid) < 120 and RNAStructure_boolean:
        nucleotides, ss, mfe = RNAStructure_fold(nucleic_acid)
        ss = structure_extension(ss)

    else:

        #The output of this set of for loops
        for x in range(0, k - 3):
            for ii in range(x * 40, 120 + x * 40):
                nucleic_acid_list.append(nucleic_acid[ii])
            nuke_structures=''.join(nucleic_acid_list)
            nucleic_acid_lists.append(nuke_structures)
            nucleic_acid_list = []
        for ii in range((k - 2) * 40, len(nucleic_acid)):
            nucleic_acid_list.append(nucleic_acid[ii])
        nucleic_acid_lists.append(''.join(nucleic_acid_list))

        #I reset the 'nucleic_acid_list' variable as an empty list because
            #I wanted to reuse the variable name
        nucleic_acid_list = []

        #Performing the folding and placing in the output
        if RNAStructure_boolean:
            for i in nucleic_acid_lists:
                nucleotides, ss, mfe = RNAStructure_fold(i)
                nucleic_acid_list.append(i + '.' + ss + ',' + str(mfe))
        else:
            for i in nucleic_acid_lists:
                (ss, mfe) = ViennaRNA_fold(i)
                ss = structure_extension(ss)
                nucleic_acid_list.append(i + ',' + ss + ',' + str(mfe))

        return nucleic_acid_list

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
    print('/////////////////////////////////////')
    print('energy_structures')
    print(energy_structures[0])
    print('energy_structures[0]')
    print(energy_structures[0][0])
    print('energy_structures[0].split(' ')')
    print(energy_structures[0][0].split(' '))
    print('/////////////////////////////////////')
    MFE_f = lambda input: (input[-1],float(input[0].split(' ')[2]))
    energy_structures = list(map(MFE_f,energy_structures))
    energy_structures.sort(key = lambda x: x[1])

    #collecting the MFE structure, the variables may be misnomers and prepping them for return from the main function
    MFE_structure = energy_structures[-1]
    structure = MFE_structure[0]
    MFE = MFE_structure[-1]
    return nucleotides, structure, MFE

#Input: path to ViennaRNA output
#Output: secondary structure and minimum free energy
def ViennaRNA_from_file(file):
    txt = open(file)
    txt = txt.readlines()
    #ViennaRNA outputs can have an extra space in the brackets with
        #the minimum free energy values
    ss_mfe = txt[-1].rstrip('\n').split(' ')
    ss, mfe = ss_mfe[0], float(ss_mfe[-1].replace('(','').replace(')',''))
    return ss, mfe

#Input: nucleotide sequence
#Output: nucletides, structure, minimum free energy according to RNAstructure
def RNAStructure_fold(nucleotide_sequence):
    #Opening a file for inserting nucleic acid sequence and inserting
    txt = open('subdirectory/misc_outputs/RNAStructure_nucleic_acid.txt', "w")
    txt.write(nucleotide_sequence)
    txt = open('subdirectory/misc_outputs/RNAStructure_nucleic_acid.txt')

    #Executing bash commands to

    os.system('./RNA_structure.sh')

    #Collecting nucleotides, dot bracket structure, and free energy
    nucleotides, structure, MFE = RNAStructure_from_file(
        'subdirectory/misc_outputs/RNAStructure_bracket_output.txt')

    return nucleotides, structure, MFE

#Input: nucleotide_sequence
#Output: folded secondary structure and minimum free energy
def ViennaRNA_fold(nucleotide_sequence):
    #Opening a file for inserting nucleic acid sequence and inserting
    txt = open('subdirectory/misc_outputs/ViennaRNA_nucleic_acid.txt',"w")
    txt.write(nucleotide_sequence)
    txt = open('subdirectory/misc_outputs/ViennaRNA_nucleic_acid.txt')

    #Executing bash commands
    os.system('RNAfold subdirectory/misc_outputs/ViennaRNA_nucleic_acid.txt > '
              'subdirectory/misc_outputs/ViennaRNA_dot_bracket_output.txt')

    #Acquiring the secondary structure and minimum free energy from the
        #ViennaRNA output file
    ss, mfe = ViennaRNA_from_file('subdirectory/misc_outputs/ViennaRNA_dot_bracket_output.txt')

    return ss, mfe

#Input: transcript ID, dictionaries of structues, dictionary of transcripts
#Output: None, multiprocessing operates using global variables, so you simply add the information
    #to a dictionary
def slow_calibration_fold(transcript, structure_dict, transcript_dict, RNAStructure_boolean, inTandem_boolean):
    #Based on user input, different dictionary outputs will be returned
    if not RNAStructure_boolean and not inTandem_boolean:
        (ss, mfe) = ViennaRNA_fold(transcript_dict[transcript])
        structure_dict[transcript] = (ss, round(mfe, 2))
    elif RNAStructure_boolean:
        nucleotides, ss, mfe = RNAStructure_fold(transcript_dict[transcript])
        structure_dict[transcript] = (ss, round(mfe, 2))
    else:
        (v_ss, v_mfe) = ViennaRNA_fold(transcript_dict[transcript])
        nucleotides, r_ss, r_mfe = RNAStructure_fold(transcript_dict[transcript])
        structure_dict[transcript] = (v_ss, round(v_mfe, 2), r_ss, round(r_mfe))



#Input: path to transcriptome fasta
#Output: a dictionary with transcript names as keys and nucleotides as values
def txt_to_transcript_dict(file_name):
    #Opening the file
    txt = open(file_name)
    txt = txt.readlines()

    #Generating a dictionary with keys given by transcript sequence names
    transcript_dict = {}

    for i, j in enumerate(txt):
        if '>' in j:
            key = j
            #I opted to name the value as '' rather than None, for no particular reason
            transcript_dict[key] = ''
        else:
            #appending the nucleotide row to the value iteratively
            transcript_dict[key] += j.rstrip('\n')
    return transcript_dict

def calibration_folding(transcript_dict, calibration_arg, sample_size_arg,
                        RNAStructure_boolean, inTandem_boolean, calibration_output,
                        fast_fold_time, output_name):
    #Retrieving transcript names
    transcript_list=list(set(transcript_dict))

    #In case one wants to modify the code to indicate the number of transcripts which have
        #been folded since the code started
    transcript_count = 0

    #Basic rationale of the process:
        #The for loop will start, spit out the transcript nucleotides
        # and get the folding algorithm running
        #If the molecule can be folded within the time period, a statement will be printed
        #to indicate the outcome, otherwise, the loop will finish up after 20 seconds just after
        #the process terminates and spit out a figure or some user specified time period

    if calibration_arg:
        slow_fold_nucleotide_count=[]
        slow_fold_nucleotide_time=[]

    #Checking the structure CSV to determine which transcripts are already folded
    try:
        file=open(output_name, newline='')
        data=csv.reader(file)
        keys=[i[0] for i in data]
    #If the file hasn't been generated, it's to move forward without it
        #for every iteration of the transcripts, each transcript is checked to
        #determine if they are in the list of keys from the csv
    #I get an error when I try opening the empty file, so the try, except
        #statement opens up an empty list as a substitute
    except:
        keys=[]

    algorithm_type = algorithm_type_test(RNAStructure_boolean)



    loop_structure_dict={}

    #The sample_size_arg allows the user to select a limited number of
        #transcripts, this can be convenient for calibration
    for transcript in tqdm(transcript_list[0:sample_size_arg]):

        if transcript in keys:
            continue

        #Opening the csv for writing
        with open(output_name, 'a') as txt_object:
            writer_object = writer(txt_object)

            #Timing the amount of time necessary for folding by capturing the
                #present moment quant
            start_fold_time = time.perf_counter()

            #Defining a new process object and naming a dictionary to retain the
                #structure
            manager = Manager()
            structure_dict=manager.dict()

            #Defining the process using the target function
            action_process = Process(target=slow_calibration_fold,
                                     args=(transcript, structure_dict,
                                           transcript_dict, RNAStructure_boolean, inTandem_boolean))
            #starting the process
            action_process.start()
            #Defining the amount of time the process can operate before timing out
            action_process.join(timeout=fast_fold_time)
            #terminating the process if it takes too long
            action_process.terminate()
            #Adding to the total number of transcripts folded since the program
                #began executing
            transcript_count += 1

            #Tells you how many transcripts have been spit out
            try:


                loop_structure_dict[transcript]=structure_dict.values()[0]
                #Outputs the dot brackets, not a whole lot else, you may want to add to the
                    #'structure dict' dictionary object
                #Calculate the fold time
                fold_time=(time.perf_counter()-start_fold_time)

                if not inTandem_boolean:
                    #The calibration still outputs its folded structures to the csv
                    writer_object.writerow([transcript, 'slow fold', algorithm_type, transcript_dict[transcript],
                                        structure_dict.values()[0][0],
                                        structure_dict.values()[0][-1]])

                else:
                    # The calibration still outputs its folded structures to the csv for ViennaRNA
                    writer_object.writerow([transcript, 'slow fold', 'ViennaRNA', transcript_dict[transcript],
                                            structure_dict.values()[0][0],
                                            structure_dict.values()[0][1]])

                    # The calibration still outputs its folded structures to the csv for RNAStructure
                    writer_object.writerow([transcript, 'slow fold', 'RNAStructure', transcript_dict[transcript],
                                            structure_dict.values()[0][2],
                                            structure_dict.values()[0][-1]])

                txt_object.close()
                #Generating two lists, one with the length of the transcripts
                    #the other with the amount of time it takes to fold them
                    #Folding is the most time intensive process and the rest
                    #is nearly instantaneous
                if calibration_arg:
                    slow_fold_nucleotide_count.append(len(transcript_dict[transcript]))
                    slow_fold_nucleotide_time.append(fold_time)

            #If the folding took to long, we switch to the fast_folding mode and
                #return the new list
            except:
                if not inTandem_boolean:
                    nuke_list=fast_fold(transcript_dict, transcript, RNAStructure_boolean)
                    #Contatins the nucleotide strings, their corresponding dot-brackets, and their minimum free energies
                        #all organized into a list
                    writer_object.writerow([transcript, 'fast fold', algorithm_type] + nuke_list)
                else:
                    nuke_list = fast_fold(transcript_dict, transcript, False)
                    # Contatins the nucleotide strings, their corresponding dot-brackets, and their minimum free energies
                    # all organized into a list
                    writer_object.writerow([transcript, 'fast fold', 'ViennaRNA'] + nuke_list)

                    nuke_list = fast_fold(transcript_dict, transcript, True)
                    # Contatins the nucleotide strings, their corresponding dot-brackets, and their minimum free energies
                    # all organized into a list
                    writer_object.writerow([transcript, 'fast fold', 'RNAStructure'] + nuke_list)

                txt_object.close()

    #Generating a png of the plot of the folding time and nucleotide lengths
    if calibration_arg:
        plt.plot(slow_fold_nucleotide_count, slow_fold_nucleotide_time, 'o')
        plt.savefig('subdirectory/misc_outputs/test_plot.png')

        #Fitting a curve to the data to calculate how long it would take to
            #fold the transcriptome
        model = np.polyfit(slow_fold_nucleotide_count, slow_fold_nucleotide_time, 2)

        #If the file containing the polynomial curve values exists, it will be deleted
            #this keeps the processes further down the pipeline simple
        dir_list = os.listdir()
        if calibration_arg in dir_list:
            os.remove(calibration_arg)

        #Opening a file to hold the values for the polynomyial curve
        with open(calibration_output, 'a') as txt_object:
            writer_object = writer(txt_object)
            writer_object.writerow(model)
            txt_object.close()

def transcript_folding(sample_size, transcript_dict, nucleic_acid_length, RNAStructure_boolean, inTandem_boolean, output_name):

    #Checking the CSV file to verify that a nucleic acid has been folded
    try:
        file = open(output_name, newline='')
        data = csv.reader(file)
        keys = [i[0] for i in data]
    except:
        keys = []

 
    #Checking algorithm type
    algorithm_type = algorithm_type_test(RNAStructure_boolean)

    #Setting variables to calculate the estimated completion time
    transcript_count = 0
    total_folding_time = 0

    if sample_size:
        #Defining the sample size for folding
        transcripts = list(set(transcript_dict))[0:sample_size]
    else:
        transcripts = list(set(transcript_dict))
    
    #Skipping the transcripts which have already been folded
    transcripts = [transcript for transcript in transcripts 
                   if transcript not in keys]


    #Looping through the transcript names
    for transcript in tqdm(transcripts):

        start_fold_time=time.perf_counter()

        #counting the number of folded transcripts
        transcript_count+=1

        # Opening the csv to begin adding folded structures
        with open(output_name, 'a') as txt_object:
            #naming a writer object
            writer_object = writer(txt_object)

            #The script is set to change folding modes depending on transcript length
            if len(transcript_dict[transcript]) > nucleic_acid_length:

                if inTandem_boolean:
                    # Fast folding is performed with ViennaRNA
                    nucleotide_structure_list = fast_fold(transcript_dict, transcript, False)

                    # Writing the fast fold to a csv
                    writer_object.writerow([transcript, 'fast fold', 'ViennaRNA'] + nucleotide_structure_list)

                    # Fast folding is performed with RNAStructure
                    nucleotide_structure_list = fast_fold(transcript_dict, transcript, True)

                    # Writing the fast fold to a csv
                    writer_object.writerow([transcript, 'fast fold', 'RNAStructure'] + nucleotide_structure_list)

                    txt_object.close()
                else:
                    #Fast folding is performed
                    nucleotide_structure_list = fast_fold(transcript_dict, transcript, RNAStructure_boolean)

                    #Writing the fast fold to a casv
                    writer_object.writerow([transcript, 'fast fold', algorithm_type] + nucleotide_structure_list)
                    txt_object.close()

            elif len(transcript_dict[transcript]) <= nucleic_acid_length:

                #choosing between ViennaRNA, RNAStructure, or both depending on user inputs
                if RNAStructure_boolean and not inTandem_boolean:
                    nucleotides, ss, mfe = RNAStructure_fold(transcript_dict[transcript])
                    writer_object.writerow \
                        ([transcript, 'slow fold', algorithm_type, transcript_dict[transcript], ss, mfe])
                    txt_object.close()
                elif inTandem_boolean:
                    (ss, mfe) = ViennaRNA_fold(transcript_dict[transcript])
                    writer_object.writerow \
                        ([transcript, 'slow fold', 'ViennaRNA', transcript_dict[transcript], ss, mfe])
                    nucleotides, ss, mfe = RNAStructure_fold(transcript_dict[transcript])
                    writer_object.writerow \
                        ([transcript, 'slow fold', 'RNAStructure', transcript_dict[transcript], ss, mfe])
                    txt_object.close()

                else:
                    
                    (ss, mfe) = ViennaRNA_fold(transcript_dict[transcript])

                    writer_object.writerow \
                        ([transcript, 'slow fold', algorithm_type, transcript_dict[transcript], ss, mfe])
                    txt_object.close()

        fold_time = (time.perf_counter() - start_fold_time)
        total_folding_time += fold_time        
        folding_rate = transcript_count/total_folding_time
        if transcript_count % 10 == 0:
            print('Estimated folding time remaining (hrs): ' 
                  + (str(folding_rate*len(transcripts)/(3600))))



