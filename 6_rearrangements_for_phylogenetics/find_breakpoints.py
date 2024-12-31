#!/usr/bin/env python3
#%%
import pandas as pd
from collections import Counter
import argparse
#%%
def read_ref_buscos(query_file):
    chr2busco, busco2chr = {}, {} # for each chr, which buscos are on it & for each busco which chr is it on
    chr2busco_pos = {} # for each chr, whats its buscos midpos. Store buscos as a list of tuples.
    with open(query_file, 'r') as query:
        for line in query:
            if 'Complete' in line:
                cols = line.split('\t')
                busco, chr ,start, end = cols[0], cols[2], int(cols[3]), int(cols[4])
                busco2chr[busco] = chr
             #  midpos = (int(start) + int(end)) / 2
                if chr in chr2busco.keys():
                    current_buscos = chr2busco[chr]
                    current_buscos.append(busco)
                    current_busco_pos = chr2busco_pos[chr]
                    current_busco_pos.append((busco, start, end))
                else:
                    chr2busco[chr] = [busco]
                    chr2busco_pos[chr] = [(busco, start, end)]

    return(chr2busco, busco2chr, chr2busco_pos)

def read_query_buscos(query_file, busco2ref_chr):
    chr2busco, busco2chr = {}, {} # for each chr, which buscos are on it & for each busco which chr is it on
    ref2query_chr = {} # for each query chr, which ref chr does each busco on it belong to 
    with open(query_file, 'r') as query:
        for line in query:
            if 'Complete' in line:
                cols = line.split('\t')
                busco, chr = cols[0], cols[2]
                busco2chr[busco] = chr
                try:
                    ref_chr = busco2ref_chr[busco] # look up which ref chr the busco is on - means we're only using query buscos that are also present in ref
                    if ref_chr in ref2query_chr.keys():
                        current_query_chr = ref2query_chr[ref_chr]
                        current_query_chr.append(chr)
                        ref2query_chr[ref_chr] = current_query_chr
                    else:
                        ref2query_chr[ref_chr] = [chr]
                    if chr in chr2busco.keys():
                        current_buscos = chr2busco[chr]
                        current_buscos.append(busco)
                    else:
                        chr2busco[chr] = [busco]
                except KeyError: # if busco isnt in ref then no need to add it into the input dicts
                    continue
    return(chr2busco, busco2chr, ref2query_chr)

def filter_ref_buscos(ref_chr2busco, busco2ref_chr, ref_chr2busco_pos, min_buscos):
    filter_chr = []
    filter_buscos = []
    for ref_chr, buscos in ref_chr2busco.items():
        if len(buscos) < min_buscos:
            filter_chr.append(ref_chr)
            print('Removed ref chr due to low number of markers:', ref_chr)
            filter_buscos = filter_buscos + buscos
    [ref_chr2busco.pop(i) for i in filter_chr] # filter dict
    [ref_chr2busco_pos.pop(i) for i in filter_chr] # filter dict
    [busco2ref_chr.pop(i) for i in filter_buscos] # filter dict
    return(ref_chr2busco, busco2ref_chr, ref_chr2busco_pos)

#filter_query_chr to only keep sequences that have at least X number of markers.
def filter_query_buscos(query_chr2busco, ref2query_chr_dict, min_buscos):
    filter_chr = [] # query_chr to filter out
    for chr, busco_list in query_chr2busco.items():
        if len(busco_list) < min_buscos:
            filter_chr.append(chr)
            print('Removed query chr due to low number of markers:', chr)
    [query_chr2busco.pop(i) for i in filter_chr] # filter dict

    for ref, query_chr_list in ref2query_chr_dict.items():
        for i in filter_chr:
            if i in query_chr_list:
                filtered_query_chr = [j for j in query_chr_list if j != i] 
                ref2query_chr_dict[ref] = filtered_query_chr
    return(query_chr2busco, ref2query_chr_dict)


def filter_ref2_query_chr(ref2query_chr_dict, min_matching_buscos): # using min_matching_buscos as default
    for ref_chr, query_chrs in ref2query_chr_dict.items():
        query_chr_to_remove = []
        if len(set(query_chrs)) > 1:  # if 1 ref chr contains elements of >1 query chr
            for i in set(query_chrs):
                per_ref_chr = (query_chrs.count(i) / len(query_chrs)*100)
            # print(ref_chr, i, per_ref_chr)
                if query_chrs.count(i) < min_matching_buscos: # if under two buscos match to a given a query chr
                    print('Filtering out: ', ref_chr, i, 'as number matching BUSCOs is: ', query_chrs.count(i))
                    query_chr_to_remove.append(i)
           # print('Query chr to remove:', query_chr_to_remove)
            filtered_query_chr = set(query_chrs) - set(query_chr_to_remove)
            ref2query_chr_dict[ref_chr] = list(filtered_query_chr) # update dict
    return(ref2query_chr_dict)

def output_breakpoints(ref_chrs, ref2query_chr_dict, busco2query_chr, ref_chr2busco, prefix):
    output_file = prefix + '.tsv'
    desired_span = 250000 # make an argument
    with open(output_file, 'w') as output:
        output.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % ("ref_chr", "start_query_chr", "current_query_chr", "breakpoint_start_busco", "breakpoint_end_busco", "breakpoint_start", "breakpoint_end", "breakpoint_span", "enlarged_breakpoint_start", "enlarged_breakpoint_end", "enlarged_breakpoint_span"))
        for ref_chr in ref_chrs:
                    query_chr_list = ref2query_chr_dict[ref_chr]
                    query_buscos = []
                    ref_buscos = ref_chr2busco[ref_chr]
                    for i in query_chr_list: # contains all buscos on query_chr that have >=1 busco on the ref_chr
                        query_buscos = query_buscos + query_chr2busco[i]
                    common_buscos = [value for value in ref_buscos if value in query_buscos]
                    ref_busco_pos = ref_chr2busco_pos[ref_chr]
                    ref_busco_pos_filt = [item for item in ref_busco_pos if item[0] in common_buscos] # filter pos to just keep intersecting buscos
                    ref_busco_pos_filt = sorted(ref_busco_pos_filt, key = lambda x: x[1])
                    counter = 1
                    for i in ref_busco_pos_filt:
                        if counter == 1:
                            start_query_chr = busco2query_chr[i[0]]
                            last_busco = i[0]
                            last_pos = i[2]
                        else: 
                            current_query_chr = busco2query_chr[i[0]]
                            if start_query_chr !=  current_query_chr: # if this busco belongs to a different query chr to the last busco, a switch in chr has occured
                                print(ref_chr, '\t', start_query_chr, '\t', last_busco)
                                print(ref_chr, '\t', current_query_chr, '\t', i[0])
                                # breakpoint is between these two chr
                                breakpoint_start_busco = last_busco
                                breakpoint_end_busco = i[0]
                                breakpoint_start = last_pos
                                breakpoint_end = i[1] # get start position of busco
                                breakpoint_span = breakpoint_end - breakpoint_start
                                print(breakpoint_start, breakpoint_end)
                                print('Breakpoint span: ', str(breakpoint_span/1000), ' kb')
                                print(' ')
                                if breakpoint_span < 0: # need to broaden breakpoint, do so by using end position of busco2 instead of start
                                    breakpoint_end = i[2]
                                    breakpoint_span = breakpoint_end - breakpoint_start
                                    print('Adjusted breakpoint span: ', str(breakpoint_span/1000), ' kb')
                                if breakpoint_span < desired_span:
                                    distance_to_add = desired_span - breakpoint_span
                                    enlarged_breakpoint_start = breakpoint_start - int(distance_to_add/2) # need start/end to be integers in downstream uses
                                    enlarged_breakpoint_end = breakpoint_end + int(distance_to_add/2)
                                    enlarged_breakpoint_span = enlarged_breakpoint_end - enlarged_breakpoint_start
                                    print('Enlarged breakpoint span: ', str(enlarged_breakpoint_span/1000), ' kb')
                                else:
                                    enlarged_breakpoint_start = breakpoint_start
                                    enlarged_breakpoint_end = breakpoint_end
                                    enlarged_breakpoint_span = breakpoint_span
                                output.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % (ref_chr, start_query_chr, current_query_chr, breakpoint_start_busco, breakpoint_end_busco, breakpoint_start, breakpoint_end, breakpoint_span, enlarged_breakpoint_start, enlarged_breakpoint_end, enlarged_breakpoint_span))
                                start_query_chr = current_query_chr # update starting query_chr for next iteration
                            last_busco = i[0]
                            last_pos = i[2] # get end position of busco
                        counter = counter + 1
    return()

def get_dominant_query_chr_and_last_pos_of_window(ref_busco_pos_filt, busco2query_chr, window_end_index, window_size, previous_dominant_chr_in_window, query_chr_list):
    for j in range(window_end_index-window_size,window_end_index, 1):
        record = ref_busco_pos_filt[j] # record has chr, busco and end_position of busco
        query_chr = busco2query_chr[record[0]] # change to explicitly using index of dict
        query_chr_list.append(query_chr)
    # Count the occurrences of each term
    query_chr_counts = Counter(query_chr_list) # how many genes per query chr are there in this window?
    # Find the maximum count
    max_count = max(query_chr_counts.values())
    query_chr_with_max_count = [chr for chr, count in query_chr_counts.items() if count == max_count] # get all query_chr with the max count (more than one may have the max..)
    if len(query_chr_with_max_count) == 1: # i.e. only one chr has the max count
        dominant_chr_in_window = query_chr_with_max_count[0]
        count = 1
        for j in range(window_end_index-window_size, window_end_index, 1):
            record = ref_busco_pos_filt[j] # record has chr, busco and end_position of busco
            query_chr = busco2query_chr[record[0]] # change to explicitly using index of dict
            if query_chr == dominant_chr_in_window: 
                if count == 1:
                    first_busco = record[0]
                    first_pos = record[1]
                last_busco = record[0]
                last_pos = record[2]
                count = count + 1
    #    print(j, dominant_chr_in_window, last_busco, last_pos)
    elif previous_dominant_chr_in_window == None:  # if more than one chr has the same max count, and this is the first window, lets assign the last gene belonging to a max query_chr as the dominant chr
        count = 1
        for j in range(window_end_index-window_size, window_end_index, 1):
            record = ref_busco_pos_filt[j] # record has chr, busco and end_position of busco
            query_chr = busco2query_chr[record[0]] # change to explicitly using index of dict
            if query_chr in query_chr_with_max_count: 
                if count == 1: # i.e. if this is the first gene of the dominant chr
                    first_busco = record[0]
                    first_pos = record[1]
                dominant_chr_in_window = query_chr
                last_busco = record[0]
                last_pos = record[2]
                count = count + 1
    else:
        print('more than one chr have same count!')
        dominant_chr_in_window = previous_dominant_chr_in_window # if more than one chr has the same max count, lets consider the current dominant chr to still be dominant
        count = 1
        for j in range(window_end_index-window_size, window_end_index, 1):
            record = ref_busco_pos_filt[j] # record has chr, busco and end_position of busco
            query_chr = busco2query_chr[record[0]] # change to explicitly using index of dict
            if query_chr == dominant_chr_in_window: 
                if count == 1: # .e. if this is the first gene of the dominant chr
                    first_busco = record[0]
                    first_pos = record[1]
                last_busco = record[0]
                last_pos = record[2]
                count = count + 1
    return(dominant_chr_in_window, first_busco, first_pos, last_busco, last_pos, set(query_chr_list))      

def get_last_pos_of_specific_chr_in_window(ref_busco_pos_filt, busco2query_chr, window_end_index, window_size, chr_of_interest):
    for j in range(window_end_index-window_size, window_end_index, 1):
        record = ref_busco_pos_filt[j] # record has chr, busco and end_position of busco
        query_chr = busco2query_chr[record[0]] # change to explicitly using index of dict
        if query_chr == chr_of_interest: 
            last_busco = record[0]
            last_pos = record[2] # this way, automatically the last gene belonging to this query_chr will be assigned as last_busco and last_pos
    return(last_busco, last_pos)      


def get_first_pos_of_specific_chr_in_window(ref_busco_pos_filt, busco2query_chr, window_end_index,window_size, chr_of_interest):
    count = 1
    for j in range(window_end_index-window_size, window_end_index, 1):
        record = ref_busco_pos_filt[j] # record has chr, busco and end_position of busco
        query_chr = busco2query_chr[record[0]] # change to explicitly using index of dict
        if query_chr == chr_of_interest: 
            if count == 1:
                first_busco = record[0]
                first_pos = record[1] # this way, automatically the last gene belonging to this query_chr will be assigned as last_busco and last_pos
            count = count + 1
    return(first_busco, first_pos)  

def output_breakpoints(breakpoint_list, prefix):
    output_file = prefix + '.breakpoints.tsv'
    with open(output_file, 'w') as output:
        output.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % ("ref_chr", "start_query_chr", "current_query_chr", "breakpoint_start_busco", "breakpoint_end_busco", "breakpoint_start", "breakpoint_end", "breakpoint_span"))
        for i in breakpoint_list:
            output.write(("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n") % (i[0], i[1], i[2], i[3], i[4], i[5], i[6], i[7]))


def find_breakpoints(ref_chrs, ref2query_chr_dict, ref_chr2busco, ref_chr2busco_pos, query_chr2busco,window_size=3):
    breakpoint_list = []
    for ref_chr in ref_chrs:
            #ref_chr = 'OW569322.1'
        previous_chr_list = []
        previous_dominant_chr_in_window = None
        query_chr_list = ref2query_chr_dict[ref_chr]
        query_buscos = []
        ref_buscos = ref_chr2busco[ref_chr]
        for i in query_chr_list: # contains all buscos on query_chr that have >=1 busco on the ref_chr
            query_buscos = query_buscos + query_chr2busco[i]
        common_buscos = [value for value in ref_buscos if value in query_buscos]
        ref_busco_pos = ref_chr2busco_pos[ref_chr]
        ref_busco_pos_filt = [item for item in ref_busco_pos if item[0] in common_buscos] # filter pos to just keep intersecting buscos
        ref_busco_pos_filt = sorted(ref_busco_pos_filt, key = lambda x: x[1])
        counter = 1
        for window_end_index in range(window_size, len(ref_busco_pos_filt), window_size): # iterate over each window, using i to specify the end of each window
            query_chr_list = []
            dominant_chr_in_window, first_busco, first_pos, last_busco, last_pos, chr_list = get_dominant_query_chr_and_last_pos_of_window(ref_busco_pos_filt, busco2query_chr, window_end_index, window_size, previous_dominant_chr_in_window, query_chr_list)
            if (dominant_chr_in_window != previous_dominant_chr_in_window) and (previous_dominant_chr_in_window != None): # i.e if the dominant query chr has switched between two windows, we have a breakpoint
             #   print('Switch between:', previous_dominant_chr_in_window, dominant_chr_in_window)
                if (previous_dominant_chr_in_window not in chr_list) and (dominant_chr_in_window not in previous_chr_list): # i.e not a single gene in this window belongs to the previous dominant chr, and not a single gene of the new dominant chr was not in the previous window
                  #  print('yas!')
                    # then we can just take the last position from previous window to be the start of the breakpoint, and the start position from this window to be the end of the breakpoint
                    breakpoint_start_busco = previous_last_busco
                    breakpoint_start_pos = previous_last_pos
                    breakpoint_end_busco = first_busco
                    breakpoint_end_pos = first_pos

                elif previous_dominant_chr_in_window in chr_list: # we need to get last location of this particular chr
                    breakpoint_start_busco, breakpoint_start_pos = get_last_pos_of_specific_chr_in_window(ref_busco_pos_filt, busco2query_chr, window_end_index, window_size, previous_dominant_chr_in_window)
                    breakpoint_end_busco = first_busco
                    breakpoint_end_pos = first_pos

                elif dominant_chr_in_window in previous_chr_list: # the final option is that at least one gene of the new dominant query chr was in the previous window, if so the breakpoint occured in the last window
                    breakpoint_start_busco, breakpoint_start_pos = get_last_pos_of_specific_chr_in_window(ref_busco_pos_filt, busco2query_chr, window_end_index-window_size, window_size, previous_dominant_chr_in_window)
                    breakpoint_end_busco, breakpoint_end_pos = get_first_pos_of_specific_chr_in_window(ref_busco_pos_filt, busco2query_chr, window_end_index-window_size, window_size, dominant_chr_in_window)

                breakpoint_start = (previous_dominant_chr_in_window, breakpoint_start_busco, breakpoint_start_pos)
                breakpoint_end = (dominant_chr_in_window, breakpoint_end_busco, breakpoint_end_pos)
                breakpoint_span = breakpoint_end_pos - breakpoint_start_pos
             #   print('Breakpoint span:', breakpoint_span)
                if breakpoint_span < 0:
                    print('WARNING: breakpoint span under 0 bp, consider altering via start/end locations of genes.')
                    # see previous scripts of way to fix this.
                breakpoint_tuple = (ref_chr, previous_dominant_chr_in_window, dominant_chr_in_window, breakpoint_start_busco, breakpoint_end_busco, breakpoint_start_pos, breakpoint_end_pos, breakpoint_span)
                breakpoint_list.append(breakpoint_tuple)
            previous_dominant_chr_in_window = dominant_chr_in_window # update to dominant query chr
            previous_last_busco = last_busco
            previous_last_pos = last_pos
            previous_chr_list = chr_list
    return(breakpoint_list)
    
 #%%
if __name__ == "__main__":
    SCRIPT = "lep_fusion_fission_finder.py"
    # argument set up
    parser = argparse.ArgumentParser()

    # Add required string arguments
    parser.add_argument("-q", "--query_file", type=str, help="Busco table of query species")
    parser.add_argument("-r", "--reference_file", type=str, help="Busco table of reference species")
    parser.add_argument("-o", "--output_prefix", type=str, help="Output prefix")

    # Add optional integer arguments
    parser.add_argument("-w", "--window_size", type=int, default=3, help="Minimum number of buscos per shared segment (default: 3)")
    parser.add_argument("-m", type=int, default=3, help="Minimum number of buscos per chr (default: 3)")
    parser.add_argument("-s", type=int, default=3, help="Minimum number of buscos shared by pair of chromosomes (default: 3)")

    # Parse the arguments
    args = parser.parse_args()
#%%
#    query_file = '../../../../Scripts/lep_fusion_fission_finder/reference_data/Merian_elements_full_table.tsv'
#    ref_file = 'Polyommatus_surakovi.tsv'
#    output_prefix = 'TEST'
#    min_matching_buscos, min_buscos_per_chr = 3, 3
    ref_chr2busco, busco2ref_chr, ref_chr2busco_pos = read_ref_buscos(args.reference_file)
    ref_chr2busco, busco2ref_chr, ref_chr2busco_pos = filter_ref_buscos(ref_chr2busco, busco2ref_chr, ref_chr2busco_pos, args.m)
    query_chr2busco, busco2query_chr, ref2query_chr_dict = read_query_buscos(args.query_file, busco2ref_chr)
    query_chr2busco, ref2query_chr_dict = filter_query_buscos(query_chr2busco, ref2query_chr_dict, args.m)
    ref2query_chr_dict = filter_ref2_query_chr(ref2query_chr_dict, args.s)
    ref_chrs = set(ref2query_chr_dict.keys())
    query_chrs  = ref2query_chr_dict.values()
    breakpoint_list = find_breakpoints(ref_chrs, ref2query_chr_dict, ref_chr2busco, ref_chr2busco_pos, query_chr2busco,args.window_size)
#%% 
    output_breakpoints(breakpoint_list, args.output_prefix)
