


import yaml
import pandas as pd

with open('deeptools_heatmaps.yaml', 'r') as file:
    obj = yaml.safe_load(file)


    all_info_as_dictionary = {}
for plot, plot_dict in obj.items():
    # get plot name
    plot_name = plot_dict['name']
    all_info_as_dictionary[plot_name] = {}
    beds_list = []
    beds_label_list = []
    bigwig_list = []
    bigwig_label_list = []
    
    
    ################### get bed files  ###################
    # get bed files in pipeline samples
    for sample_id in plot_dict['peaks_to_plot']['samples_in_this_pipeline']:
        bed_file = sample_id + '.bed'
        if plot_dict['peaks_to_plot']['samples_in_this_pipeline'][sample_id]['label'] == '':
            region_label = sample_id
        else:
            region_label = plot_dict['peaks_to_plot']['samples_in_this_pipeline'][sample_id]['label']
        
        if bed_file != '':
            beds_list.append (bed_file)
            beds_label_list.append (region_label)
        
    # get bed files in local path 
    for local_path in plot_dict['peaks_to_plot']['bed_files_in_local_path']:
        bed_file = local_path
        if plot_dict['peaks_to_plot']['bed_files_in_local_path'][local_path]['label'] == '':
            region_label = local_path
        else:
            region_label = plot_dict['peaks_to_plot']['bed_files_in_local_path'][local_path]['label']
        
        if bed_file != '':
            beds_list.append (bed_file)
            beds_label_list.append (region_label)
        
    # get bed files in remote server 
    for remote in plot_dict['peaks_to_plot']['bed_files_in_remote']:
        bed_file = remote
        if plot_dict['peaks_to_plot']['bed_files_in_remote'][remote]['label'] == '':
            region_label = remote
        else:
            region_label = plot_dict['peaks_to_plot']['bed_files_in_remote'][remote]['label']
        
        if bed_file != '':
            beds_list.append (bed_file)
            beds_label_list.append (region_label)
    

    
    ################### get bigwig files  ###################
    # get bigwig files in pipeline samples
    for sample_id in plot_dict['bigwig_files_for_signal']['samples_in_this_pipeline']:
        bigwig_file = sample_id + '.bigwig'
        if plot_dict['bigwig_files_for_signal']['samples_in_this_pipeline'][sample_id]['label'] == '':
            region_label = sample_id
        else:
            region_label = plot_dict['bigwig_files_for_signal']['samples_in_this_pipeline'][sample_id]['label']
        
        if bigwig_file != '':
            bigwig_list.append (bigwig_file)
            bigwig_label_list.append (region_label)
        
    # get bigwig files in local path 
    for local_path in plot_dict['bigwig_files_for_signal']['bigwig_files_in_local_path']:
        bigwig_file = local_path
        if plot_dict['bigwig_files_for_signal']['bigwig_files_in_local_path'][local_path]['label'] == '':
            region_label = local_path
        else:
            region_label = plot_dict['bigwig_files_for_signal']['bigwig_files_in_local_path'][local_path]['label']
        
        if bigwig_file != '':
            bigwig_list.append (bigwig_file)
            bigwig_label_list.append (region_label)
        
    # get bigwig files in remote server 
    for remote in plot_dict['bigwig_files_for_signal']['bigwig_files_in_remote']:
        bigwig_file = remote
        if plot_dict['bigwig_files_for_signal']['bigwig_files_in_remote'][remote]['label'] == '':
            region_label = remote
        else:
            region_label = plot_dict['bigwig_files_for_signal']['bigwig_files_in_remote'][remote]['label']
        
        if bigwig_file != '':
            bigwig_list.append (bigwig_file)
            bigwig_label_list.append (region_label)
    
    # add files and labels as lines with spacing between them (ready for injection in command line of deeptools!)
    all_info_as_dictionary[plot_name]['bedFiles'] = ' '.join(beds_list)
    all_info_as_dictionary[plot_name]['bedFiles_labels'] = ' '.join(beds_label_list)
    all_info_as_dictionary[plot_name]['bigwigFiles'] = ' '.join(beds_list)
    all_info_as_dictionary[plot_name]['bigwigFiles_labels'] = ' '.join(beds_label_list)


data_as_df = pd.DataFrame.from_dict(all_info_as_dictionary, orient='index')
data_as_df.index.name = 'plot_name'

data_as_df.to_csv('data.csv')