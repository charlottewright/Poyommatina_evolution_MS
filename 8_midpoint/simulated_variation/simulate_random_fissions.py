#!/usr/bin/env python3
#%%
import random
import matplotlib.pyplot as plt
import numpy as np
import statistics
import argparse


def generate_random_number(minimum_number, max_number):
    # Generate 20 random numbers between 1 and n
    random_number = random.randint(minimum_number, max_number)
    return random_number

def find_chr_for_value(value, ranges_dict):
    for key, (start, end) in ranges_dict.items():
        if start <= value <= end:
            return key
    return None  # If no range matches

def simulate_fissions(chr2accumulative_length, accumulative_length, number_fissions, minimum_chr_length):
    simulated_chr2accumulative_length = chr2accumulative_length.copy()
    count = 0
    while count < number_fissions:
        random_number = generate_random_number(minimum_chr_length, (accumulative_length-minimum_chr_length)) # generate a number of minimum the min chr lengyh and max the total chr lengh minus the min chr length
        count = count + 1
        # Find which chr has been cut!
        hit_chr = find_chr_for_value(random_number, simulated_chr2accumulative_length)
        hit_chr_start = simulated_chr2accumulative_length[hit_chr][0]
        hit_chr_stop = simulated_chr2accumulative_length[hit_chr][1]
        new_chr_1 = (hit_chr_start, random_number)
        new_chr_2 = ((random_number+1), hit_chr_stop)
        chr_1_length = random_number - hit_chr_start
        chr_2_length = hit_chr_stop - random_number
        if chr_1_length >= minimum_chr_length and chr_2_length >= minimum_chr_length:
            count = count + 1
            del simulated_chr2accumulative_length[hit_chr]
            new_chr_1_ID = hit_chr + '_' + str(count) + 'a'
            new_chr_2_ID = hit_chr + '_' + str(count) + 'b'
            simulated_chr2accumulative_length[new_chr_1_ID] = new_chr_1
            simulated_chr2accumulative_length[new_chr_2_ID] = new_chr_2
    return(simulated_chr2accumulative_length)

parser = argparse.ArgumentParser(description="A script that takes four inputs.")

parser.add_argument('-f', '--file', type=str, required=True, help='Path to the input chr lengths bed file')
parser.add_argument('-n', '--number', type=int, required=True, help='Number of breaks to induce in chromosomes')
parser.add_argument('-e', '--exclude', type=str, required=False, help='Chromosome to exclude (optional)')
parser.add_argument('-p', '--prefix', type=str,  required=True, help='Prefix for output files')
parser.add_argument('-m', '--minimum', type=int, default=0, help='Minimum length of a chromosome produced by a split (default:0 ; i.e. no minimum)')

parser.add_argument('-i', '--iterations', type=int, default=100, help='Number of iterations to simulate (default: 100)')
parser.add_argument('-o', '--observed', type=int, required=True, help='Observed scaled proportional stdev')

# Parse arguments
args = parser.parse_args()

chr_lengths_file = args.file
chr_to_exclude = args.exclude
number_fissions = args.number
number_iterations = args.iterations
prefix = args.prefix
observed_scaled_prop_sd = args.observed
minimum_chr_length = args.minimum
#chr_lengths_file = '/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus/Raw_data/chr_lengths/Cyaniris_semiargus.seqlen.bed'
#chr_to_exclude = 'LR994547.1' # Z-chromosome
#number_iterations = 100
#number_fissions = 64  # number of fissions to simulate
#prefix = 'test'
#observed_scaled_prop_sd =  21.28628 
# %%
chr2length, chr2accumulative_length = {}, {}
accumulative_length = 0 
with open(chr_lengths_file, 'r') as file:
    for line in file:
        cols = line.split('\t')
        chr, start, stop = cols[0], int(cols[1]), int(cols[2].replace("\n", ""))
        if chr != chr_to_exclude:
            chr2length[chr] = stop
            accumulative_start = accumulative_length + start
            accumulative_length = accumulative_length + stop
            chr2accumulative_length[chr] = (accumulative_start, accumulative_length)

# accumulative length is the number of fissions to simulate
# pick n random numbers, each of which represents a fission site

results_scaled_stdev_length = []
for _ in range(0, number_iterations):
    simulated_chr2accumulative_length = simulate_fissions(chr2accumulative_length, accumulative_length, number_fissions, minimum_chr_length)
    # now lets convert from tuples of (start, end) back to just chr lengths
    simulated_chr_lengths = {}
    updated_dict = {}
    for chr in simulated_chr2accumulative_length.keys():  
        start, end = simulated_chr2accumulative_length[chr]
        simulated_chr_lengths[chr] = end - start  # calculate length
    simulated_chr_lengths_list = list(simulated_chr_lengths.values())
    simulated_prop_chr_lengths_list = [(x / accumulative_length)*100 for x in simulated_chr_lengths_list] # convert to proportional length
    mean_prop_length = statistics.mean(simulated_prop_chr_lengths_list)
    mean_length = statistics.mean(simulated_chr_lengths_list)
    stdev_prop_length = statistics.stdev(simulated_prop_chr_lengths_list)
    stdev_length = statistics.stdev(simulated_chr_lengths_list)
    scaled_stdev_length = (stdev_prop_length/mean_prop_length)*100
    results_scaled_stdev_length.append(scaled_stdev_length)
    print(f'mean_length',mean_length)
    print(f'stdev_length',stdev_length)
    #plt.hist(simulated_prop_chr_lengths_list)
    #plt.show() 

# %%
plt.hist(results_scaled_stdev_length, bins=10, edgecolor='black') # plot histogram
plt.xlim(left=0) # Set the x-axis to start at 0
plt.axvline(x=observed_scaled_prop_sd, color='red', linestyle='--', linewidth=2) # Add a vertical line at observed freq
output_plot = prefix + '.histogram_exp_vs_obs_variation_in_chr_lengths.pdf'
plt.savefig(output_plot, format='pdf')  # Save the plot to a PDF file

#%%
sim_file = prefix + '.simulated_values_of_scaled_proportional_stdev.tsv'
with open(sim_file, 'w') as file:
    for item in results_scaled_stdev_length:
        file.write(str(item) + '\n')  # Write each item followed by a newline
