#!/bin/python3

import functions

import argparse

#TODO: SET TO CALCULATE TIME BY SUBTRACTING THE FOLDED COLLECTION OF TRANSCRIPTS (TRIVIAL), #SEE IF YOU CAN MAKE

#TODO: Devise a way to operate RNAStructure and ViennaRNA in tandem

#TODO: set the calibration output csv file to avoid adding additional lines

#input: transcript ID, structure_dictionary, transcript_dictionary
#Output: technically nothing, you may want to keep the structure dictionary as a global variable, rather
    #than a subprocess

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
                                                            'folded in iterations of 120 nucleotides', type=int,
                                                            nargs='?', const=1800, default=1800)

    parser.add_argument('-rs', '--RNAStructure', help = 'switches the folding algorithm to RNAStructure',
                        action='store_true')

    args=parser.parse_args()

    transcript_dict=functions.txt_to_transcript_dict(args.input)

    if args.calibrate:
        functions.calibration_folding(transcript_dict, args.calibrate, args.sample_size, args.RNAStructure)

    if not args.calibrate and args.fold_transcriptome:
        functions.transcript_folding(args.sample_size, transcript_dict, args.nucleic_acid_length,
                                     args.RNAStructure)

if __name__ == '__main__':
    main()
