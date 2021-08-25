import pandas as pd
import numpy as np
import os 
import argparse

def get_files(directory):
    current_folder = directory
    dir_list = os.listdir(directory)
    for folder in dir_list:
        if 'RNA' in folder:
            mrna = pd.read_csv(directory + '/' + folder + '/' + 'suntag.loc3', header = None, sep = '  ', engine = 'python')
            print('mrna:' + folder)
        if 'single' in folder:
            singles = pd.read_csv(directory + '/' + folder + '/' + 'ng.loc3', header = None, sep = '  ', engine = 'python')
            print('singles:' + folder)           
        if 'T' in folder:
            bright = pd.read_csv(directory + '/' + folder + '/' + 'ng.loc3', header = None, sep = '  ', engine = 'python')
            print('bright:' + folder)
           
    print('xy_voxels:')
    xy_voxels = float(input())
    print('z_voxels:')
    z_voxels = float(input())
    
    mrna.columns = ['x', 'y', 'z', 'intensity', 'distance_to_rna']
    mrna = mrna[mrna.intensity != -8]
    mrna['x'] = mrna['x']*xy_voxels
    mrna['y'] = mrna['y']*xy_voxels
    mrna['z'] = mrna['z']*z_voxels
    print('mrna:')
    print(mrna.head(10))
    print(mrna.shape)
    
    singles.columns = ['x', 'y', 'z', 'intensity', 'distance_to_rna']
    singles = singles[singles.intensity != -8]
    singles['x'] = singles['x']*xy_voxels
    singles['y'] = singles['y']*xy_voxels  
    singles['z'] = singles['z']*z_voxels
    print('singles:')
    print(singles.head(10))
    print(singles.shape)
    
    
    bright.columns = ['x', 'y', 'z', 'intensity', 'distance_to_rna']
    bright['x'] = bright['x']*xy_voxels
    bright['y'] = bright['y']*xy_voxels
    bright['z'] = bright['z']*z_voxels
    print('brights:')
    print(bright.head(10))
    print(bright.shape)
    
    return(mrna, singles, bright, current_folder)
    
    
def drop_brights(singles, bright):    
    singles = singles[~singles.isin(bright)].dropna()
    print('just singles:')
    print(singles.head(10))
    print(singles.shape)
    return(singles)
        

def get_distances(mrna, singles, bright):
    for i in singles.index:
        all_dis = {}
        for j in mrna.index:
            dis_x = np.power((singles.loc[i, 'x'] - mrna.loc[j, 'x']), 2)
            dis_y= np.power((singles.loc[i, 'y'] - mrna.loc[j, 'y']), 2)
            dis_z= np.power((singles.loc[i, 'z'] - mrna.loc[j, 'z']), 2)
            distance = np.sqrt(dis_x + dis_y + dis_z)
            all_dis[distance] = j
        min_distance = min(all_dis.keys())
        singles.loc[i, 'distance_to_rna'] = min_distance
        min_distance_rna = all_dis[min_distance]
        singles.loc[i, 'closest_rna'] = min_distance_rna
    
    for i in bright.index:
        all_dis = {}
        for j in mrna.index:
            dis_x = np.power((bright.loc[i, 'x'] - mrna.loc[j, 'x']), 2)
            dis_y= np.power((bright.loc[i, 'y'] - mrna.loc[j, 'y']), 2)
            dis_z= np.power((bright.loc[i, 'z'] - mrna.loc[j, 'z']), 2)
            distance = np.sqrt(dis_x + dis_y + dis_z)
            all_dis[distance] = j
        min_distance = min(all_dis.keys())
        bright.loc[i, 'distance_to_rna'] = min_distance
        min_distance_rna = all_dis[min_distance]
        bright.loc[i, 'closest_rna'] = min_distance_rna
    
    return(singles, bright)
    

def remove_edges(singles, bright):
    singles = singles[(singles.x > 0.5) & (singles.y > 0.5) & (singles.z > 0.5)]
    singles = singles[(singles.x < (max(singles.x) - 0.5)) & (singles.y < (max(singles.y) - 0.5)) & (singles.z < (max(singles.z) - 0.5))]
    print('singles without edges:')
    print(singles.head(10))
    print(singles.shape)
    
    bright = bright[(bright.x > 0.5) & (bright.y > 0.5) & (bright.z > 0.5)]
    bright = bright[(bright.x < (max(bright.x) - 0.5)) & (bright.y < (max(bright.y) - 0.5)) & (bright.z < (max(bright.z) - 0.5))]
    print('brights without edges:')
    print(bright.head(10))
    print(bright.shape)
    
    return(singles, bright)
    
    
def get_colocs(singles, bright):
    coloc_singles = singles[singles.distance_to_rna <= 0.3]
    print('colocalised singles:')
    print(coloc_singles.head(10))
    print(coloc_singles.shape)
    non_coloc_singles = singles[singles.distance_to_rna > 0.3]
    print('non colocalised singles:')
    print(non_coloc_singles.head(10))
    print(non_coloc_singles.shape)
    coloc_bright = bright[bright.distance_to_rna <= 0.3]
    print('colocalised brights:')
    print(coloc_bright.head(10))
    print(coloc_bright.shape)
    non_coloc_bright = bright[bright.distance_to_rna > 0.3]
    print('non colocalised brights:')
    print(non_coloc_bright.head(10))
    print(non_coloc_bright.shape)
    
    return(coloc_singles, non_coloc_singles, coloc_bright, non_coloc_bright)
    
    
def get_median(non_coloc_singles):
    single_median = non_coloc_singles['intensity'].median()
    print('median intensity of non-colocalised singles:')
    print(single_median)
    return(single_median)
    
def get_ribosome_numbers(coloc_bright, coloc_singles, single_median):
    conversion_factor = 1
    ribosomes_bright = (coloc_bright['intensity']/single_median) / conversion_factor
    ribosomes_singles = (coloc_singles['intensity']/single_median)/conversion_factor
    number_of_ribosomes = np.concatenate((np.array(ribosomes_bright), np.array(ribosomes_singles)))
    print('number of ribosomes:')
    print(number_of_ribosomes)
    print(number_of_ribosomes.shape)
    return(number_of_ribosomes)
    
def make_file(number_of_ribosomes, current_folder):
    df = pd.DataFrame(number_of_ribosomes, columns = ['num_ribo'])
    df['embryo'] = current_folder
    return df
    

parser = argparse.ArgumentParser()
parser.add_argument('directory', help = 'folder with Airlocalise files in')
args = parser.parse_args()

mrna, singles, bright, current_folder = get_files(args.directory)
singles = drop_brights(singles, bright)
singles, bright = get_distances(mrna, singles, bright)
singles, bright = remove_edges(singles, bright)
coloc_singles, non_coloc_singles, coloc_bright, non_coloc_bright = get_colocs(singles, bright)
single_median = get_median(non_coloc_singles)
number_of_ribosomes = get_ribosome_numbers(coloc_bright, coloc_singles, single_median)
df = make_file(number_of_ribosomes, current_folder)
df.to_csv(current_folder + '/num_ribo.csv')



