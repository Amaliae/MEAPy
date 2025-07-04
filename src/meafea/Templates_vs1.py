#!/usr/bin/env python
# coding: utf-8

# In[100]:


import numpy as np
import pandas as pd
import os


# In[101]:


def Rules():
    
    print('In Excel the rule is /n Name should have ending mapMEA.xlsx') 
    print('In benchling formats csv table should have well labels and column labels and name should contain param identifier, it should be in dpath')
   
    print('csv and hdf5 folders should be in the same directory')
    
    print('descritpion excel file should In mainpath of experiment ')
    
    


# In[4]:


def make_template():
    
    wellid=[ 0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  0.,  6.,
        6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  6.,  1.,  1.,
        1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  1.,  7.,  7.,  7.,
        7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  7.,  2.,  2.,  2.,  2.,
        2.,  2.,  2.,  2.,  2.,  2.,  2.,  2.,  8.,  8.,  8.,  8.,  8.,
        8.,  8.,  8.,  8.,  8.,  8.,  8.,  3.,  3.,  3.,  3.,  3.,  3.,
        3.,  3.,  3.,  3.,  3.,  3.,  9.,  9.,  9.,  9.,  9.,  9.,  9.,
        9.,  9.,  9.,  9.,  9.,  4.,  4.,  4.,  4.,  4.,  4.,  4.,  4.,
        4.,  4.,  4.,  4., 10., 10., 10., 10., 10., 10., 10., 10., 10.,
       10., 10., 10.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,  5.,
        5.,  5., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11., 11.,
       11., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17., 17.,
       23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 23., 16.,
       16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 16., 22., 22.,
       22., 22., 22., 22., 22., 22., 22., 22., 22., 22., 15., 15., 15.,
       15., 15., 15., 15., 15., 15., 15., 15., 15., 21., 21., 21., 21.,
       21., 21., 21., 21., 21., 21., 21., 21., 14., 14., 14., 14., 14.,
       14., 14., 14., 14., 14., 14., 14., 20., 20., 20., 20., 20., 20.,
       20., 20., 20., 20., 20., 20., 13., 13., 13., 13., 13., 13., 13.,
       13., 13., 13., 13., 13., 19., 19., 19., 19., 19., 19., 19., 19.,
       19., 19., 19., 19., 12., 12., 12., 12., 12., 12., 12., 12., 12.,
       12., 12., 12., 18., 18., 18., 18., 18., 18., 18., 18., 18., 18.,
       18., 18.]
    
    channelid=[0,   1,   2,   3,   4,   5,   6,   7,   8,   9,  10,  11,  12,
        13,  14,  15,  16,  17,  18,  19,  20,  21,  22,  23,  24,  25,
        26,  27,  28,  29,  30,  31,  32,  33,  34,  35,  36,  37,  38,
        39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,
        52,  53,  54,  55,  56,  57,  58,  59,  60,  61,  62,  63,  64,
        65,  66,  67,  68,  69,  70,  71,  72,  73,  74,  75,  76,  77,
        78,  79,  80,  81,  82,  83,  84,  85,  86,  87,  88,  89,  90,
        91,  92,  93,  94,  95,  96,  97,  98,  99, 100, 101, 102, 103,
       104, 105, 106, 107, 108, 109, 110, 111, 112, 113, 114, 115, 116,
       117, 118, 119, 120, 121, 122, 123, 124, 125, 126, 127, 128, 129,
       130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142,
       143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 153, 154, 155,
       156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168,
       169, 170, 171, 172, 173, 174, 175, 176, 177, 178, 179, 180, 181,
       182, 183, 184, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194,
       195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207,
       208, 209, 210, 211, 212, 213, 214, 215, 216, 217, 218, 219, 220,
       221, 222, 223, 224, 225, 226, 227, 228, 229, 230, 231, 232, 233,
       234, 235, 236, 237, 238, 239, 240, 241, 242, 243, 244, 245, 246,
       247, 248, 249, 250, 251, 252, 253, 254, 255, 256, 257, 258, 259,
       260, 261, 262, 263, 264, 265, 266, 267, 268, 269, 270, 271, 272,
       273, 274, 275, 276, 277, 278, 279, 280, 281, 282, 283, 284, 285,
       286, 287]
    
    welllabel=['A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1', 'A1',
       'A1', 'B1', 'B1', 'B1', 'B1', 'B1', 'B1', 'B1', 'B1', 'B1', 'B1',
       'B1', 'B1', 'A2', 'A2', 'A2', 'A2', 'A2', 'A2', 'A2', 'A2', 'A2',
       'A2', 'A2', 'A2', 'B2', 'B2', 'B2', 'B2', 'B2', 'B2', 'B2', 'B2',
       'B2', 'B2', 'B2', 'B2', 'A3', 'A3', 'A3', 'A3', 'A3', 'A3', 'A3',
       'A3', 'A3', 'A3', 'A3', 'A3', 'B3', 'B3', 'B3', 'B3', 'B3', 'B3',
       'B3', 'B3', 'B3', 'B3', 'B3', 'B3', 'A4', 'A4', 'A4', 'A4', 'A4',
       'A4', 'A4', 'A4', 'A4', 'A4', 'A4', 'A4', 'B4', 'B4', 'B4', 'B4',
       'B4', 'B4', 'B4', 'B4', 'B4', 'B4', 'B4', 'B4', 'A5', 'A5', 'A5',
       'A5', 'A5', 'A5', 'A5', 'A5', 'A5', 'A5', 'A5', 'A5', 'B5', 'B5',
       'B5', 'B5', 'B5', 'B5', 'B5', 'B5', 'B5', 'B5', 'B5', 'B5', 'A6',
       'A6', 'A6', 'A6', 'A6', 'A6', 'A6', 'A6', 'A6', 'A6', 'A6', 'A6',
       'B6', 'B6', 'B6', 'B6', 'B6', 'B6', 'B6', 'B6', 'B6', 'B6', 'B6',
       'B6', 'C6', 'C6', 'C6', 'C6', 'C6', 'C6', 'C6', 'C6', 'C6', 'C6',
       'C6', 'C6', 'D6', 'D6', 'D6', 'D6', 'D6', 'D6', 'D6', 'D6', 'D6',
       'D6', 'D6', 'D6', 'C5', 'C5', 'C5', 'C5', 'C5', 'C5', 'C5', 'C5',
       'C5', 'C5', 'C5', 'C5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5', 'D5',
       'D5', 'D5', 'D5', 'D5', 'D5', 'C4', 'C4', 'C4', 'C4', 'C4', 'C4',
       'C4', 'C4', 'C4', 'C4', 'C4', 'C4', 'D4', 'D4', 'D4', 'D4', 'D4',
       'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'D4', 'C3', 'C3', 'C3', 'C3',
       'C3', 'C3', 'C3', 'C3', 'C3', 'C3', 'C3', 'C3', 'D3', 'D3', 'D3',
       'D3', 'D3', 'D3', 'D3', 'D3', 'D3', 'D3', 'D3', 'D3', 'C2', 'C2',
       'C2', 'C2', 'C2', 'C2', 'C2', 'C2', 'C2', 'C2', 'C2', 'C2', 'D2',
       'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2', 'D2',
       'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1', 'C1',
       'C1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1', 'D1',
       'D1', 'D1']
    
    channellabel=[23, 24, 12, 13, 21, 22, 31, 32, 33, 42, 34, 43, 23, 24, 12, 13, 21,
       22, 31, 32, 33, 42, 34, 43, 23, 24, 12, 13, 21, 22, 31, 32, 33, 42,
       34, 43, 23, 24, 12, 13, 21, 22, 31, 32, 33, 42, 34, 43, 23, 24, 12,
       13, 21, 22, 31, 32, 33, 42, 34, 43, 23, 24, 12, 13, 21, 22, 31, 32,
       33, 42, 34, 43, 23, 24, 12, 13, 21, 22, 31, 32, 33, 42, 34, 43, 23,
       24, 12, 13, 21, 22, 31, 32, 33, 42, 34, 43, 23, 24, 12, 13, 21, 22,
       31, 32, 33, 42, 34, 43, 23, 24, 12, 13, 21, 22, 31, 32, 33, 42, 34,
       43, 23, 24, 12, 13, 21, 22, 31, 32, 33, 42, 34, 43, 23, 24, 12, 13,
       21, 22, 31, 32, 33, 42, 34, 43, 31, 32, 42, 43, 33, 34, 23, 24, 13,
       22, 12, 21, 31, 32, 42, 43, 33, 34, 23, 24, 13, 22, 12, 21, 31, 32,
       42, 43, 33, 34, 23, 24, 13, 22, 12, 21, 31, 32, 42, 43, 33, 34, 23,
       24, 13, 22, 12, 21, 31, 32, 42, 43, 33, 34, 23, 24, 13, 22, 12, 21,
       31, 32, 42, 43, 33, 34, 23, 24, 13, 22, 12, 21, 31, 32, 42, 43, 33,
       34, 23, 24, 13, 22, 12, 21, 31, 32, 42, 43, 33, 34, 23, 24, 13, 22,
       12, 21, 31, 32, 42, 43, 33, 34, 23, 24, 13, 22, 12, 21, 31, 32, 42,
       43, 33, 34, 23, 24, 13, 22, 12, 21, 31, 32, 42, 43, 33, 34, 23, 24,
       13, 22, 12, 21, 31, 32, 42, 43, 33, 34, 23, 24, 13, 22, 12, 21]
    
    
    
    template_df=pd.DataFrame({'Well Label':welllabel, 'Well ID': wellid, 'Channel ID': channelid, 'Channel Label':channellabel})
    
    return template_df


# In[5]:


def make_excel_experiment_old(path, Name):
    
    cols=['Well Label', 'Compound', 'Concentration_µM', 'Vehicle', 'Treatment',
       'Experiment_Excel', 'hdf5_data_folder', 'raw_data_folder', 'save_folder',
       'PlateID', 'Date_neuronal_plating', 'Date_compound_treatment',
       'Cell_origin', 'Genotype', 'CellLine',
       'Cell_culture_type', 'Neuron_density_k_per_well',
       'Astrocyte_density_k_per_well', 'Media_type', 'Comment', 'Well']
    
    dfexcel=pd.DataFrame(columns=cols)
    
   
    dfexcel['Well Label']=[let+str(num) for let in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] for num in np.arange(1, 13, 1).astype('int')]
    
    
    
    dfexcel.to_excel(os.path.join(path, Name+'mapMEA.xlsx'))
    
    
    return dfexcel


# In[6]:


def make_excel_experiment(path, Name):
   # Define the column names
   cols = ['Well Label', 'Compound', 'Concentration_µM', 'Vehicle', 'Treatment',
           'Experiment_Excel', 'hdf5_data_folder', 'raw_data_folder', 'save_folder',
           'PlateID', 'Date_neuronal_plating', 'Date_compound_treatment',
           'Cell_origin', 'Genotype', 'CellLine',
           'Cell_culture_type', 'Neuron_density_k_per_well',
           'Astrocyte_density_k_per_well', 'Media_type', 'Comment', 'Well']
   # Create a DataFrame with five empty rows
   empty_rows = pd.DataFrame({col: ['']*5 for col in cols})
   # Create a DataFrame with column names as the sixth row
   col_names_row = pd.DataFrame([cols], columns=cols)
   # Create another DataFrame with the required data starting from the 7th row
   data = {'Well Label':[let+str(num) for let in ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'] for num in np.arange(1, 13, 1).astype('int')]}
   df_data = pd.DataFrame(data)
   # Concatenate the empty rows, the column names row, and the data
   dfexcel = pd.concat([empty_rows, col_names_row, df_data], ignore_index=True)
   # Save to Excel file
    
    
   dfexcel.to_excel(os.path.join(path, Name+'mapMEA.xlsx'), index=False, header=False)


   return dfexcel.iloc[6:, :]


# In[7]:


def find_excel(dpath, identifier):
    
    '''mapMEA for excel and any other param for identifier'''

    excelpath, excel= '', pd.DataFrame()
    excel_files=[item for item in list(os.listdir(dpath)) if identifier in item]
    
    if len(excel_files)>0:
            
        excelpath=os.path.join(dpath, excel_files[0])
        
        try:

            excel=pd.read_excel(excelpath, skiprows=5)
            
        except:
            
            try:
                
                excel=pd.read_csv(excelpath)
            except:
                print('No exel or benchling')
                
                
            
            
        
    return excelpath, excel
        
        
            
            
    


# In[15]:


def template_excel(path, Name):
    '''
    path: The path to the excel file folder without the file-specific name
    Name: The name a newly created excel file should have if no file was found under the specifed directory
    '''
    
    template_df=make_template()
    
    excels=[item for item in list(os.listdir(path)) if 'mapMEA' in item]

    
    if len(excels)>0:
        
        excelpath, excel=find_excel(path, 'mapMEA')
        
        print(excel.columns)
        
        print('found excel')
                                        
    else:
        
        make_excel_experiment(path, Name)
        
        print('made an excel')
        
        excelpath, excel=find_excel(path, 'mapMEA')
        
        
        
        
        
       
        
    
        
    constant_cols=['Experiment_Excel', 'hdf5_data_folder', 'raw_data_folder', 'Save_folder', 'PlateID', 'Date_neuronal_plating', 'Date_compound_treatment',
                   'Cell_origin', 'Genotype', 'CellLine',
                   'Cell_culture_type', 'Neuron_density_k_per_well',
                   'Astrocyte_density_k_per_well', 'Media_type', 'Comment']

    var_cols=['Well Label', 'Compound', 'Concentration_µM', 'Vehicle', 'Treatment', 'Experiment_Excel', 'PlateID', 'Cell_origin', 'Genotype', 'CellLine',
           'Cell_culture_type', 'Neuron_density_k_per_well',
           'Astrocyte_density_k_per_well']
    
    # Specification of  column read-in where row value of template is taken only if it is present, otherwise the first column value will be used as default. 
    # This gives the flexibility to only enter one default value at the top but to still modify individual values below if desired.
    col_to_fill_list = [col for col in excel.columns if col not in ["Well Label", "Well"]] # Define the list of columns to check >> We only omit the grey columns that are not to be modified.

    
    # Loop over each row in the DataFrame
    for index, row in excel.iterrows():
        # Loop over each column in the green_col_list
        for column in col_to_fill_list:
            # Define the default value
            default_value = excel.loc[0, column]
            # Check if the cell is empty
            if pd.isnull(row[column]):
                # Fill the cell with the default value
                excel.loc[index, column] = default_value  

    
    
    merged_df=pd.merge(template_df, excel, on='Well Label', how='left')
    #print(merged_df, 'merged_df')       
    return merged_df


# In[9]:


def Benchling_layout(dpath, excel, Param):
    
    '''if thete is enchling csv file in dir, you can use it for feeding info to excel file'''
    
    benchlingpath, fg=find_excel(dpath, Param)
    
    if len(benchlingpath)>0:
        
        fg=pd.read_csv(benchlingpath, index_col='Unnamed: 0')

        for value in np.unique(fg.values):

            rowindex, colindex=np.where(fg==value)

            welllab=fg.index[rowindex]+fg.columns[colindex].astype('str')



            excel.loc[excel['Well Label'].isin(welllab)==True, Param]=value

    
    return excel
    
    


# In[ ]:





# TODO
# 
# Fill in info on culture, manually or automoatically 

# In[ ]:




