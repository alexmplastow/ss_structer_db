import argparse
import functions
from decimal import Decimal
import matplotlib.pyplot as plt

import statistics

#TODO:You may need to adjust for the processor usage

#Defining the parser arguments
def main():
    parser=argparse.ArgumentParser(description='test CLI')
    parser.add_argument('-it', '--input', help = 'Input transcriptome file', type=str)
    parser.add_argument('-ico', '--input_calibration_output', help = 'Input calibration output', type=str)
    args = parser.parse_args()

    #I'm pretty sure I opted to fold the transcripts here so I could calculate
        #the amount of time it would take to fold the whole
    transcript_dict=functions.txt_to_transcript_dict(args.input)

    transcript_set=set(transcript_dict)

    transcript_length_dict={}
    for transcript in transcript_set:
        transcript_length_dict[transcript]=len(transcript_dict[transcript])

    #Opening the txt file and pulling out the polynomial
    txt = open(args.input_calibration_output)
    txt = txt.readlines()
    txt[-1] = txt[-1].rstrip('\n')

    #splitting the polynomial into a list and altering the data type
    polynomial = txt[-1].split(',')
    polynomial = [float(Decimal(i)) for i in polynomial]

    a=polynomial[0]
    b=polynomial[1]
    c=polynomial[-1]


    #Defining two lists to contain the maximum nucleotide length before
	#switching folding type, the second contains the amount of time needed
	#to fold the transcriptome
    set_of_n_max = []
    set_of_t_N = []
    
    #I do not need to iterate through all 100000 possibilities
    for n_max in range(0,10000,50):
        
        #Defining a vector to hold the times for transcriptome folding
        time_per_folding=[]
        for transcript in tuple(transcript_dict):
            n=transcript_length_dict[transcript]
            
            #Calculating the amount of time needed for the fast folded specimens
            if n > n_max:
                time_per_folding.append(((n/40)*(a*(120**2)+b*(120)+c)))
            
            #Calculating the time needed for the slow folded specimens
            elif n <= n_max:
                time_per_folding.append((a * (n ** 2) + b*(n) + c))
        
        #Summing the folding times and converting to hours
        set_of_t_N.append(sum(time_per_folding)/(60*60))
        
        #appending the folding lengths
        set_of_n_max.append(n_max)
    
    #Plotting
    plt.plot(set_of_n_max, set_of_t_N)
    plt.xlabel('# of bases to be considered before switching into fast fold mode')
    plt.ylabel('transcriptome folding time (hr)')
    plt.savefig('folding_time_plot_test.png')



#Executing the function
if __name__ == '__main__':
    main()
