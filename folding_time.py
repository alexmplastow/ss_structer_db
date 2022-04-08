import argparse
import functions
from decimal import Decimal
import matplotlib.pyplot as plt

def main():
    parser=argparse.ArgumentParser(description='test CLI')
    parser.add_argument('-it', '--input', help='Input transcriptome_file', type=str)
    parser.add_argument('-ico', '--input_calibration_output', help='Input calibration output', type=str)
    args = parser.parse_args()

    transcript_dict=functions.txt_to_transcript_dict(args.input)

    transcript_set=set(transcript_dict)

    transcript_length_dict={}
    for transcript in transcript_set:
        transcript_length_dict[transcript]=len(transcript_dict[transcript])

    txt = open(args.input_calibration_output)
    txt = txt.readlines()
    txt[-1] = txt[-1].rstrip('\n')

    polynomial = txt[-1].split(',')
    polynomial = [float(Decimal(i)) for i in polynomial]

    a=polynomial[0]
    b=polynomial[1]
    c=polynomial[-1]


    set_of_t_N=[]
    set_of_n_max=[]
    for i in range(0,10000,50):
        t_N=0
        n_max=i
        for transcript in set(transcript_dict):
            n=transcript_length_dict[transcript]
            if n > n_max:
                t_N+=((n/40)*(a*(120**2)+b*(120)+c))
            elif n <= n_max:
                t_N+=(a * (n ** 2) + b*(n) + c)
        set_of_t_N.append(t_N)
        set_of_n_max.append(n_max)
    plt.plot(set_of_n_max, set_of_t_N)
    plt.savefig('folding_time_plot_test.png')




if __name__ == '__main__':
    main()
