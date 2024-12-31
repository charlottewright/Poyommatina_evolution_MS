

#%%
working_dir = '/lustre/scratch122/tol/teams/blaxter/projects/lepidoptera/cw22/polyommatus/Analysis/busco_paints/'
assignments_file = 'Polyommatus_icarus.window_17_chromosome_assignments.tsv'
assignments_file_path = working_dir + assignments_file

chromosome2merians = {}
with open(assignments_file_path, 'r') as f:
        next(f) # skip header
        for line in f:
            parts = line.strip().split()  # Split line by whitespace
            chromosome, merians = parts[0], parts[2]
            merian_list = merians.split(',')
            chromosome2merians[chromosome] = merian_list

complete_locations_file = 'Polyommatus_icarus_complete_location.tsv'
complete_locations_file_path = working_dir + complete_locations_file

ouput_file = 'Polyommatus_icarus.buscos_on_correct_merian.tsv'
output_filepath = working_dir + ouput_file
correct_n, false_n = 0, 0
with open(output_filepath, 'w') as o:
    with open(complete_locations_file_path, 'r') as f:
            next(f) # skip header
            for line in f:
                parts = line.strip().split()  # Split line by whitespace
                busco, chromosome, merian = parts[0], parts[1], parts[3]  
                dominant_merians = chromosome2merians[chromosome] # get dominant merians for this chr
                if merian in dominant_merians:
                    correct_n = correct_n + 1
                    o.write("%s\t%s" % (busco, chromosome) + "\n")
                else:
                    false_n = false_n + 1
                    print(line)
# %%
false_n
#%%
correct_n
# %%
