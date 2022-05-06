#!/bin/python3

import functions

import argparse

#TODO: It doesn't look like the folded molecules are not being ignored, find the problem

#TODO: clean up your usage of the calibration argument, it does not look inherently problematic, but it should only be
    #TODO: included in the parser, but not the body of functions

#TODO: it looks like the subprocess lingers after the code finishes executing, it may not be problematic, so it may
    #TODO: be prudent not to worry about it

#TODO: fix your usage of the term 'nucleotide' in the variable names

#input: transcript ID, structure_dictionary, transcript_dictionary
#Output: technically nothing, you may want to keep the structure dictionary as a global variable, rather
    #than a subprocess

#defining the main function and providing the arguments for the CLI
def main():
    parser=argparse.ArgumentParser(description='A CLI for generating or improving on transcriptome-wide secondary'
                                               'structure datasets')

    parser.add_argument('-i','--input', help='Input transcriptome_file', type=str)

    parser.add_argument('-c', '--calibrate',
                        help='Determine how long folding will take depending on folding time of a set of nucleotides',
                        action='store_true')

    parser.add_argument('-s', '--sample_size', help='nucleic acid count for folding (default 50)', type=int,
                        nargs='?', const=50)

    parser.add_argument('-ff', '--fast_fold_time', help='nucleic acid count for folding (default 50)',
                        type=int,
                        nargs='?', const=50)

    parser.add_argument('-cof', '--calibration_output', help='output file for calibration', type=str, nargs='?',
                        default='calibration_output.csv', const='calibration_output.csv')

    parser.add_argument('-ft', '--fold_transcriptome', help='indicate if whole transcriptome should begin folding',
                        action='store_true')

    parser.add_argument('-l', '--nucleic_acid_length', help='maximum length of nucleic acids, anything above will be '
                                                            'folded in iterations of 120 nucleotides,' 
                                                            'the default is 1800 nucleotides', type=int,
                                                            nargs='?', const=1800, default=1800)

    parser.add_argument('-rs', '--RNAStructure', help = 'switches the folding algorithm to RNAStructure',
                        action='store_true')

    parser.add_argument('-it', '--inTandem', help = 'run RNAStructure and ViennaRNA concurrently',
                        action='store_true')

    args=parser.parse_args()

    #Collecting the transcript dictionary from the user input
    transcript_dict = functions.txt_to_transcript_dict(args.input)

    #Depending on the arguments passed into the program, calibration can begin here
    if args.calibrate:
        functions.calibration_folding(transcript_dict, args.calibrate,
                                      args.sample_size, args.RNAStructure, args.inTandem,
                                      args.calibration_output)

    #If the user has indicate they wish to begin folding the transcriptome
    if not args.calibrate and args.fold_transcriptome:
        functions.transcript_folding(args.sample_size, transcript_dict, args.nucleic_acid_length,
                                     args.RNAStructure, args.inTandem)

#Not sure why I have this here, but all CLIs seem to require these lines
if __name__ == '__main__':
    main()
