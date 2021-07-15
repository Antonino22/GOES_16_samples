import os
import matplotlib as mpl
import matplotlib.pyplot as plt
os.environ['PROJ_LIB'] = '/Users/anthonycrespo/anaconda3/pkgs/proj4-5.2.0-h6de7cb9_1006/share/proj'
from mpl_toolkits.basemap import Basemap
import numpy as np
from remap_copy import remap
from netCDF4 import Dataset
from osgeo import gdal
import glob
import datetime
import shutil
import warnings
import scandir

#Ignore deprecated warnings
warnings.filterwarnings("ignore")


def ext(extent1):
    extent = extent1
    if extent == 'Nation':
        zoom = [-76, -58, -51, -20]
        return zoom
    elif extent == 'Feb8-NC':
        zoom = [-64.65, -31.76, -64.275, -31.24] #Hail event for Dec 13
        return zoom
    elif extent == 'Dec14':
        zoom = [-67, -32.7, -64.1, -32.25] #Hail event for Dec 13
        return zoom
    elif extent == 'Oct24':
        zoom = [-64.6, -31.62, -63.8, -31.35] #Hail event for Dec 13
        return zoom        
    else:
        pass    

def abi_data(fname, extent, resolution, variable, Band): 
    
    H, a, b, f, longitude, x1, y1, x2, y2 = sat_variables(fname)

    grid = remap(fname, variable, extent, resolution, H, a, b, f, longitude, x1, y1, x2, y2)

    if Band >= 7:
        data = grid.ReadAsArray()
        Unit = "Brightness Temperature [K]"
    else:
        data = grid.ReadAsArray()
        a = data

        # set any negative values to zero\n"
        data[data < 0.] = 0.

        # normalize data to 1.16
        data = data/1.16

        # set any values greater than 1.0 equal to 1.0
        data[a>1.]=1.

        Unit = "Reflectance"
     
    return data     


def sat_variables(path):
    nc = Dataset(path, 'r')

    # Calculate the image extent 
    H = nc.variables['goes_imager_projection'].perspective_point_height 
    x1 = nc.variables['x_image_bounds'][0] * H
    x2 = nc.variables['x_image_bounds'][1] * H
    y1 = nc.variables['y_image_bounds'][1] * H
    y2 = nc.variables['y_image_bounds'][0] * H
    
    # Read the central longitude
    longitude = nc.variables['goes_imager_projection'].longitude_of_projection_origin

    # Read the semi major axis
    a = nc.variables['goes_imager_projection'].semi_major_axis

    # Read the semi minor axis
    b = nc.variables['goes_imager_projection'].semi_minor_axis
    
    # Flattening factor
    f = 1/nc.variables['goes_imager_projection'].inverse_flattening
    
    return H, a, b, f, longitude, x1, y1, x2, y2


def min_tb(extent1, resolution, variable, Ch13, Tb_min):
    Band = int(file_list[Ch13 + number][-57:-55])
    nc1 = Dataset(file_list[Ch13 + number],'r')
    extent = ext(extent1)

    import datetime
    time_var = nc1.time_coverage_start

    import calendar
    itimehr = time_var[11:13]
    itimemn = time_var[14:16]
    time = itimehr + itimemn

    #Get data
    data = abi_data(file_list[Ch13 + number], extent, resolution, variable, Band)

    Min = np.min(data)

    if number == 0:
        Der = 0
    else:
        Der = Min - Tb_min[number - 1]

    source = file_list[Ch13 + number] + '.aux.xml'
    destination =  file_list[Ch13 + number][0:45] + 'z_aux_files' + file_list[Ch13 + number][59:200] + '.aux.xml'

    # Move a file from the directory d1 to d2
    shutil.move(source, destination)

    # Delete aux.xml
    os.remove(destination) 

    return Min, time , Der

def Tri_diff(extent1, resolution, variable, number, Ch11, Ch14, Ch15):
    Band11 = int(file_list[Ch11 + number][-57:-55])
    Band14 = int(file_list[Ch14 + number][-57:-55])
    Band15 = int(file_list[Ch15 + number][-57:-55])
    extent = ext(extent1)

    data11 = abi_data(file_list[Ch11 + number], extent, resolution, variable, Band11)
    data14 = abi_data(file_list[Ch14 + number], extent, resolution, variable, Band14)      
    data15 = abi_data(file_list[Ch15 + number], extent, resolution, variable, Band15)   

    d11_14 = data11 - data14
    d14_15 = data14 - data15
    Tri = np.max(d11_14 - d14_15)
    max_diff = np.max(d11_14)

    source1 = file_list[Ch11 + number] + '.aux.xml'
    source2 = file_list[Ch14 + number] + '.aux.xml'
    source3 = file_list[Ch15 + number] + '.aux.xml'
    destination1 =  file_list[Ch11 + number][0:45] + 'z_aux_files' + file_list[Ch11 + number][59:200]+ '.aux.xml'
    destination2 =  file_list[Ch14 + number][0:45] + 'z_aux_files' + file_list[Ch14 + number][59:200]+ '.aux.xml'
    destination3 =  file_list[Ch15 + number][0:45] + 'z_aux_files' + file_list[Ch15 + number][59:200]+ '.aux.xml'

    # Move a file from the directory d1 to d2
    shutil.move(source1, destination1)
    shutil.move(source2, destination2)
    shutil.move(source3, destination3)

    # Delete aux.xml
    os.remove(destination1) 
    os.remove(destination2) 
    os.remove(destination3) 

    return Tri, max_diff
    

def plot_temp(Tb, Time, file_res, Der, Tri, Diff, Ch11):
    print('Now plotting!')

    # Save figure as a .png file with "save_name" to specified directory
    date = datetime.datetime.strptime(file_list[Ch11][-49:-42], '%Y%j').date().strftime('%Y%m%d')
    save_name = 'G' + file_list[Ch11][-53:-51] + '_ABI_' + date + '_' + 'rib'

    fig, ax = plt.subplots(2, 2, sharex=True)
    ax[0,0].plot(Time, Tb, color='black')
    ax[0,1].plot(Time, Der, color='black')
    ax[1,0].plot(Time, Tri, color='black')
    ax[1,1].plot(Time, Diff, color='black')
    ax[0,0].set_title("Ch13 Min Tb")
    ax[0,1].set_title("Updraft Strength")
    ax[1,0].set_title("Maximum Tri-spectral Difference")
    ax[1,1].set_title("Particle Size (Ice Crystals)")
    ax[1,0].tick_params(axis='x', which='major', labelsize=10, labelrotation=45)
    ax[1,1].tick_params(axis='x', which='major', labelsize=10, labelrotation=45)
    ax[0,0].set_ylabel('Temperature (K)')
    if date == '20201024' :
        ax[0,1].set_ylabel('Temperature (K 10 min$^{-1}$)')
    else:
        ax[0,1].set_ylabel('Temperature (K 15 min$^{-1}$)')
    ax[1,0].set_ylabel('Temperature (K)')
    ax[1,1].set_ylabel('Temperature (K)')
    ax[0,0].set_xlabel('Time (UTC)')
    ax[0,1].set_xlabel('Time (UTC)')
    ax[1,0].set_xlabel('Time (UTC)')
    ax[1,1].set_xlabel('Time (UTC)')
    plt.tight_layout()

    if date == '20180208':
        #ax[0,0].axvspan('1830', '1845', color='yellow', alpha=0.5)
        #ax[0,1].axvspan('1830', '1845', color='yellow', alpha=0.5)
        #ax[1,0].axvspan('1830', '1845', color='yellow', alpha=0.5)
        #ax[1,1].axvspan('1830', '1845', color='yellow', alpha=0.5)
        #ax[0,0].axvspan('1900', '1915', color='red', alpha=0.5)
        #ax[0,1].axvspan('1900', '1915', color='red', alpha=0.5)
        #ax[1,0].axvspan('1900', '1915', color='red', alpha=0.5)
        #ax[1,1].axvspan('1900', '1915', color='red', alpha=0.5)
        #ax[0,0].axvspan('1915', '1930', color='blue', alpha=0.5)
        #ax[0,1].axvspan('1915', '1930', color='blue', alpha=0.5)
        #ax[1,0].axvspan('1915', '1930', color='blue', alpha=0.5)
        #ax[1,1].axvspan('1915', '1930', color='blue', alpha=0.5)
        #ax[0,0].axvspan('1930', '1945', color='green', alpha=0.5)
        #ax[0,1].axvspan('1930', '1945', color='green', alpha=0.5)
        #ax[1,0].axvspan('1930', '1945', color='green', alpha=0.5)
        #ax[1,1].axvspan('1930', '1945', color='green', alpha=0.5)
        pass
    elif date == '20181214':
        ax[0,0].axvspan('0200', '0215', color='yellow', alpha=0.5)
        ax[0,1].axvspan('0200', '0215', color='yellow', alpha=0.5)
        ax[1,0].axvspan('0200', '0215', color='yellow', alpha=0.5)
        ax[1,1].axvspan('0200', '0215', color='yellow', alpha=0.5)
        #ax[0,0].axvspan('0245', '0300', color='red', alpha=0.5)
        #ax[0,1].axvspan('0245', '0300', color='red', alpha=0.5)
        #ax[1,0].axvspan('0245', '0300', color='red', alpha=0.5)
        #ax[1,1].axvspan('0245', '0300', color='red', alpha=0.5)
    elif date == '20201024':
        ax[0,0].axvspan('1810', '1820', color='yellow', alpha=0.5)
        ax[0,1].axvspan('1810', '1820', color='yellow', alpha=0.5)
        ax[1,0].axvspan('1810', '1820', color='yellow', alpha=0.5)
        ax[1,1].axvspan('1810', '1820', color='yellow', alpha=0.5)
        ax[0,0].axvspan('1820', '1830', color='red', alpha=0.5)
        ax[0,1].axvspan('1820', '1830', color='red', alpha=0.5)
        ax[1,0].axvspan('1820', '1830', color='red', alpha=0.5)
        ax[1,1].axvspan('1820', '1830', color='red', alpha=0.5)
        ax[0,0].axvspan('1910', '1920', color='green', alpha=0.5)
        ax[0,1].axvspan('1910', '1920', color='green', alpha=0.5)
        ax[1,0].axvspan('1910', '1920', color='green', alpha=0.5)
        ax[1,1].axvspan('1910', '1920', color='green', alpha=0.5)
    else:
        pass    

    plt.show()
    #plt.savefig(save_path + save_name, bbox_inches='tight', dpi=file_res)

##########################################################################################################################

# File settings
file_path = os.getcwd() + '/' + '02_08_2018' + '/'
save_path = os.getcwd() + '/'  # Directory where figures will be saved
figures = 'multiple'  # Number of figures to create: 'single' (one data file) or 'multiple' (multiple data files)

# Plot settings
file_res = 150  # DPI setting for image resolution

# Mapping settings
variable = 'CMI'
extent1 = 'Feb8-NC'
resolution = 1.
Ch11 = 0
Ch13 = 23 # 23,28,31
Ch14 = 46 # 46,56,62
Ch15 = 69 # 69,84,93
num_of_files = 23 #23,28,31 

##########################################################################################################################

if __name__ == "__main__":

    if figures == 'multiple':
        Tb_min = []
        Der = []
        Tri = []
        Diff = []
        Time = []
        file_list = []
        number = 0

        for paths, dirs, files in scandir.walk(file_path):
        #for (paths, dirs, files) in os.walk(folder):
            for file in files:
                if file.endswith(".nc"):
                    file_list.append(os.path.join(paths, file))

        # Loop through data files, making/saving a figure for each data file
        while number < num_of_files:
            print('Run # :',number)
            TB, TIME, Deriv =  min_tb(extent1, resolution, variable, Ch13, Tb_min)
            Tb_min.append(TB)
            Time.append(TIME) 
            Der.append(Deriv)
            tri, max_d =  Tri_diff(extent1, resolution, variable, number, Ch11, Ch14, Ch15)
            Tri.append(tri)
            Diff.append(max_d)
            number = number+1
                

        plot_temp(Tb_min, Time, file_res, Der, Tri, Diff, Ch11)

        print('Time')
        print(Time)
        print('Tb')
        print(Tb_min)
        print('Der')
        print(Der)
        print('Tri')
        print(Tri)
        print('Diff')
        print(Diff)
        
        print('Done!')

    else:
        pass

