import argparse
import functions
from decimal import Decimal
import matplotlib.pyplot as plt

import statistics

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


    set_of_t_N=[]
    set_of_n_max=[]
    for n_max in range(0,10000,50):
        time_per_folding=[]
        t_N=0
        for transcript in tuple(transcript_dict):
            n=transcript_length_dict[transcript]
            if n > n_max:
                t_N+=((n/40)*(a*(120**2)+b*(120)+c))
                time_per_folding.append(((n/40)*(a*(120**2)+b*(120)+c)))
            elif n <= n_max:
                t_N+=(a * (n ** 2) + b*(n) + c)
                time_per_folding.append((a * (n ** 2) + b*(n) + c))
        set_of_t_N.append(sum(time_per_folding)/(60*60))
        set_of_n_max.append(n_max)
    plt.plot(set_of_n_max, set_of_t_N)
    plt.xlabel('# of bases to be considered before switching into fast fold mode')
    plt.ylabel('transcriptome folding time (hr)')
    plt.savefig('folding_time_plot_test.png')
    lowest_time = min(set_of_t_N)
    shortest_length = tuple(iterable for iterable 
                      in tuple(zip(set_of_n_max, set_of_t_N)) 
                      if iterable[-1] = lowest_time)[0]
    print('The value which will result in the fastest folding is: ' + shortest_length)



if __name__ == '__main__':
 	   main()
