
#%%

#%%
import csv
import statistics
import seaborn as sns
import matplotlib.pyplot as plt  # Make sure to import plt from matplotlib
import pandas as pd
#%%
# Read in query chromosome lengths from index
def read_file2(file2, minimum_size):
    size_dict = {}
    with open(file2, mode='r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) == 3:
                col1, _, size = row
                if int(size) > minimum_size:
                    size_dict[col1] = int(size) 
    return size_dict

# Read in query chromosomes and assigned reference chromosomes
def group_sizes_by_ref_chr(file1, size_dict):
    ref_chr_2_lengths = {}
    with open(file1, mode='r') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) == 3:  # Only process valid rows
                query_chr, status, assigned_ref_chr = row
                if status == 'split': # only consider fissions 
                    if query_chr in size_dict:  # Look up size from the dictionary
                        size = size_dict[query_chr]
                        if assigned_ref_chr not in ref_chr_2_lengths:
                            ref_chr_2_lengths[assigned_ref_chr] = []
                        ref_chr_2_lengths[assigned_ref_chr].append(size)
    
    return ref_chr_2_lengths

def max_per_difference(sizes):
    if len(sizes) < 2:
        return float('nan')  # Not enough data to calculate a difference
    differences = [abs(x - y) for i, x in enumerate(sizes) for j, y in enumerate(sizes) if i < j]
    max_per_difference = (max(differences)/sum(sizes))*100
    return max_per_difference

# Step 3: Calculate the mean difference of sizes for each set of query chr assigned to a given ref chr 
def mean_difference(sizes):
    if len(sizes) < 2:
        return str('nan')  # Not enough values to calculate a difference
    differences = [abs(x - y) for i, x in enumerate(sizes) for j, y in enumerate(sizes) if i < j]
    return sum(differences) / len(differences)

# %%
## Inputs:

minimum_size = 1000000 # to filter out unplaced scaffolds and mitos.
#file1 = '../Analysis/LFFF/species_comparisons/Lysandra_bellargus_vs_Polyommatus_icarus.w3_chromosome_assignments.tsv'
#file2 = '../Raw_data/chr_lengths/Lysandra_bellargus.seqlen.bed'

file1 = '../Analysis/LFFF/species_comparisons/Lysandra_hispana_vs_Polyommatus_icarus.w3_chromosome_assignments.tsv'
file2 = '../Raw_data/chr_lengths/Lysandra_hispana.seqlen.bed'


#file1 = '../Analysis/LFFF/species_comparisons/Lysandra_coridon_vs_Polyommatus_icarus.w3_chromosome_assignments.tsv'
#file2 = '../Raw_data/chr_lengths/Lysandra_coridon.seqlen.bed'

#file1 = '../Analysis/LFFF/species_comparisons/Polyommatus_atlantica_vs_Polyommatus_icarus.w3_chromosome_assignments.tsv'
#file2 = '../Raw_data/chr_lengths/Polyommatus_atlantica.seqlen.bed'

#file1 = '../Analysis/LFFF/species_comparisons/Polyommatus_dorylas_vs_Polyommatus_icarus.w3_chromosome_assignments.tsv'
#file2 = '../Raw_data/chr_lengths/Polyommatus_dorylas.seqlen.bed'
Z_chr = 'OW569320.1'
output_prefix = 'Lysandra_hispana.w3.using_P_icarus_ref'

#output_prefix = 'Lysandra_coridon.w3.using_P_icarus_ref'
#%%
size_dict = read_file2(file2, minimum_size)  # Create a dictionary from file2
ref_chr_2_lengths = group_sizes_by_ref_chr(file1, size_dict)  # Group sizes by col3

differences_list = []
number_chr_per_ref_chr_list = []
result = {}
for ref_chr, sizes in ref_chr_2_lengths.items():
    mean_diff = mean_difference(sizes)
    max_diff = max_per_difference(sizes)
    total_length = sum(sizes)
    if (ref_chr != Z_chr) & (len(sizes) > 1): # has undergone a fusion! and if at least two sizes
        result[ref_chr] = {'mean_diff': mean_diff, 'max_diff':max_diff, 'total_lengths': total_length, 'prop_diff':(mean_diff/total_length)*100,'num_query_chr':len(sizes)}
        differences_list.append((mean_diff/total_length)*100)
        number_chr_per_ref_chr_list.append(len(sizes))
#%%
for ref_chr, stats in result.items():
    print(f"ref_chr: {ref_chr}, Mean Difference: {stats['mean_diff']}, 'Max Difference: {stats['max_diff']}, Total size: {stats['total_lengths']},Proportional difference: {stats['prop_diff']}")
#%%
print('Mean difference in prop length:',round(statistics.mean(differences_list),2),'%')
print('Standard deviation in prop length:',round(statistics.stdev(differences_list),2),'%')

#%%
statistics.mean(number_chr_per_ref_chr_list)
#%%
df = pd.DataFrame(result).T

# Reset the index to have a default integer index
df.reset_index(inplace=True)

df.columns = ['ref_chr','mean_diff', 'max_diff','total_lengths', 'prop_diff','num_query_chr']

# Sort the DataFrame by 'ref_chr'
#df_sorted = df.sort_values(by='ref_chr')


#%%
sns.scatterplot(data=df, x='num_query_chr', y='max_diff')

#%%
output_file = output_prefix + '.variation_in_cut_length.tsv'
df.to_csv(output_file, sep='\t', index=False)
