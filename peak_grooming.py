import pandas as pd
import numpy as np

peaks = pd.read_table('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Combined_peakCoverage.txt').dropna()

## For normalization with read depth:
# mapped_depth = pd.read_table('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Combined_peakCoverage_mappedDepth.txt')
# for col in mapped_depth.columns:
#     peaks[col] = peaks[col]/mapped_depth[col].values[0]

## Number of elements to remove after '_' delimiting: 
prefix_count = 2

## Define the mapping of old names to new labels:
part_mapping = {'Enrich': 'E', 'Input': 'I', 'Ars':'S', 'Mock':'M', 'Fraction':'F', 'Ctrl':'C', 'Control':'C', 'SG':'G3BP'}

## Modify the column names:
new_columns = []
column_occurrences = {}  # Dictionary to track column name occurrences

for col in peaks.columns:
    parts = col.split('_')
    if len(parts) > prefix_count:
        for i in range(prefix_count, len(parts)):
            if parts[i] in part_mapping:
                parts[i] = part_mapping[parts[i]]
        
        new_col = '_'.join(parts[prefix_count:])
        
        ## Additional operation: Check the number of elements after the prefix_count
        num_elements = len(new_col.split('_'))
        if num_elements == 4:
            new_col_parts = new_col.split('_')
            new_col_base = '_'.join(new_col_parts[:-1])  # Remove the last element
            
            ## Track occurrences of modified column names
            if new_col_base in column_occurrences:
                column_occurrences[new_col_base] += 1
                new_col = f'{new_col_base}_{column_occurrences[new_col_base]}'
            else:
                column_occurrences[new_col_base] = 1
                new_col = f'{new_col_base}_1'
        
        new_columns.append(new_col)
    else:
        new_columns.append(col)

peaks.columns = new_columns

## Remove the specified columns:
columns_to_remove = ['NoRT_C', 'NoTemp_C', 'G3BP_E_C', 'G3BP_I_C']
peaks = peaks.drop(columns=columns_to_remove)

## Define the desired column order:
desired_order = ['chr', 'start', 'end', 'name', 'score', 'strand', 'Nuc_F_M', 'Nuc_F_S', 'Cyto_F_M', 'Cyto_F_S',
                 'NLS_I_M', 'NLS_I_S', 'NES_I_M', 'NES_I_S', 'G3BP_I_M', 'G3BP_I_S', 'NLS_E_M', 'NLS_E_S',
                 'NES_E_M', 'NES_E_S', 'G3BP_E_M', 'G3BP_E_S']

## Create a dictionary to store the replicate counts for each column:
replicate_counts = {'Nuc_F_M': 3, 'Nuc_F_S': 3, 'Cyto_F_M': 3, 'Cyto_F_S': 3, 'NLS_I_M': 2, 'NLS_I_S': 2, 'NES_I_M': 2,
                    'NES_I_S': 2, 'G3BP_I_M': 4, 'G3BP_I_S': 5, 'NLS_E_M': 4, 'NLS_E_S': 4, 'NES_E_M': 4, 'NES_E_S': 4,
                    'G3BP_E_M': 6, 'G3BP_E_S': 7}

new_desired_order = []
for col in desired_order:
    replicate_count = replicate_counts.get(col, 1)
    if replicate_count > 1 and col not in ('chr', 'start', 'end', 'name', 'score', 'strand'):
        new_desired_order.extend([f'{col}_{i+1}' for i in range(replicate_count)])
    else:
        new_desired_order.append(col)

## Reorder columns by new order:
peaks = peaks[new_desired_order]
peaks[peaks.columns[6:]] = peaks[peaks.columns[6:]].astype(int)
peaks['TOTAL_TagCount'] = peaks[peaks.columns[6:]].sum(axis = 1)

## Get unique base names from column names (considering the first three elements after splitting)
unique_base_names = set('_'.join(col.split('_')[:3]) for col in peaks.columns[6:])

column_groups = {}
for base_name in unique_base_names:
    columns = [col for col in peaks.columns if col.startswith(base_name)]
    column_groups[base_name] = columns

## Calculate biological complexity:
BC_df = pd.DataFrame()
for group, columns in column_groups.items():
    new_column_name = f'{group}_BC'
    non_zero_counts = np.sum(peaks[columns] != 0, axis=1)  # Calculate counts using NumPy
    BC_df[new_column_name] = non_zero_counts

BC_df = BC_df[['Nuc_F_M_BC', 'Nuc_F_S_BC', 'Cyto_F_M_BC', 'Cyto_F_S_BC', 'NLS_I_M_BC', 'NLS_I_S_BC', 'NES_I_M_BC', 'NES_I_S_BC', 'G3BP_I_M_BC', 'G3BP_I_S_BC', 'NLS_E_M_BC', 'NLS_E_S_BC', 'NES_E_M_BC', 'NES_E_S_BC', 'G3BP_E_M_BC', 'G3BP_E_S_BC']]
BC_df['TOTAL_BC'] = BC_df.sum(axis = 1)

## concatenate the BC counter to the original df:
peaks = pd.concat([peaks, BC_df], axis = 1)

peaks.to_csv('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Combined_peakCoverage_groomed.txt', sep = '\t')
# peaks.to_csv('/Users/soonyi/Desktop/Genomics/CoCLIP/Analysis/Combined_peakCoverage_groomed_normalized.txt', sep = '\t', index = False)