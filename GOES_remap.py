import os
import matplotlib as mpl
import matplotlib.pyplot as plt
os.environ['PROJ_LIB'] = '/Users/anthonycrespo/anaconda3/pkgs/proj4-5.2.0-h6de7cb9_1006/share/proj'
from mpl_toolkits.basemap import Basemap
import numpy as np
from remap_copy import remap
from cpt_convert import loadCPT
from matplotlib.colors import LinearSegmentedColormap
from netCDF4 import Dataset
from matplotlib.patches import Rectangle
from osgeo import gdal
import glob
import pandas as pd
import datetime
import shutil
import warnings

#Ignore deprecated warnings
warnings.filterwarnings("ignore")

def sat_variables(path):
    nc = Dataset(path, 'r')

    # Calculate the image extent 
    H = nc.variables['goes_imager_projection'].perspective_point_height + nc.variables['goes_imager_projection'].semi_major_axis
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

def cpt(Band):
    if 8<= Band <= 10:
        cpt = {'red': ((0.0, 0.0, 0.0),(0.290, 0.263, .263),(0.385, 1.0, 1.0),(0.475, 0.443, .443),
                      (0.515, 0.0, 0.0),(0.575, 1.0, 1.0),(0.664, 1.0, 1.0),(1.0, 0.0, 0.0)),
              'green': ((0.0, 0.0, 0.0),(0.290, .513, .513),(0.385, 1.0, 1.0),(0.475, .443, .443),
                      (0.515, 0., 0.0),(0.575, 1.0, 1.0),(0.664, 0.0, 0.0),(1.0, 0.0, 0.0)),
               'blue': ((0.0, 0.0, 0.0),(0.290, .137, .137),(0.385, 1.0, 1.0),(0.475,0.694, 0.694),
                      (0.515, .451, .451),(0.552, 0.0, 0.0),(0.664, 0.0, 0.0),(1.0, 0.0, 0.0))}

    elif Band >= 11:
        cpt = {'red':((0.0, 0.0, 0.0),(.001, 1.00, 1.00),(.107, 1.00, 1.00),(.113, 0.498, 0.498),
                     (.173, 1.00, 1.00),(.179, 0.902, 0.902),(.227, 0.102, 0.102),(.233, 0.00, 0.00),
                     (.287, 0.902, 0.902),(.293, 1.00, 1.00),(.346, 1.00, 1.00),(.352, 1.00, 1.00),
                     (.406, 0.101, 0.101),(.412, 0.00, 0.00),(.481, 0.00, 0.00),(.484, 0.00, 0.00),
                     (.543, 0.00, 0.00),(.546, 0.773, 0.773),(.994, 0.012, 0.012),(.997, 0.004, 0.004),
                     (1.0, 0.0, 0.0)),
             'green':((0.0, 0.0, 0.0),(.001, 1.00, 1.00),(.107, 1.00, 1.00),(.113, 0.00, 0.00),
                     (.173, 0.498, 0.498),(.179, 0.902, 0.902),(.227, 0.102, 0.102),(.233, 0.00, 0.00),
                     (.287, 0.00, 0.00),(.293, 0.00, 0.00),(.346, 0.902, 0.902),(.352, 1.00, 1.00),
                     (.406, 1.00, 1.00),(.412, 1.00, 1.00),(.481, 0.00, 0.00),(.484, 0.00, 0.00),
                     (.543, 1.00, 1.00),(.546, 0.773, 0.773),(.994, 0.012, 0.012),(.997, 0.004, 0.004),
                     (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.00, 0.00),(.001, 1.00, 1.00),(.107, 0.00, 0.00),(.113, 0.498, 0.498),
                     (.173, 0.786, 0.786),(.179, 0.902, 0.902),(.227, 0.102, 0.102),(.233, 0.00, 0.00),
                     (.287, 0.00, 0.00),(.293, 0.00, 0.00),(.346, 0.00, 0.00),(.352, 0.00, 0.00),
                     (.406, 0.00, 0.00),(.412, 0.00, 0.00),(.481, 0.451, 0.451),(.484, 0.451, 0.451),
                     (.543, 1.00, 1.00),(.546, 0.773, 0.773),(.994, 0.012, 0.012),(.997, 0.004, 0.004),
                     (1.0, 0.0, 0.0))}
    else:
        cpt = {'red':((0.0, 0.0, 0.0),(.001, 1.00, 1.00),(.107, 1.00, 1.00),(.113, 0.498, 0.498),
                     (.173, 1.00, 1.00),(.179, 0.902, 0.902),(.227, 0.102, 0.102),(.233, 0.00, 0.00),
                     (.287, 0.902, 0.902),(.293, 1.00, 1.00),(.346, 1.00, 1.00),(.352, 1.00, 1.00),
                     (.406, 0.101, 0.101),(.412, 0.00, 0.00),(.481, 0.00, 0.00),(.484, 0.00, 0.00),
                     (.543, 0.00, 0.00),(.546, 0.773, 0.773),(.994, 0.012, 0.012),(.997, 0.004, 0.004),
                     (1.0, 0.0, 0.0)),
             'green':((0.0, 0.0, 0.0),(.001, 1.00, 1.00),(.107, 1.00, 1.00),(.113, 0.00, 0.00),
                     (.173, 0.498, 0.498),(.179, 0.902, 0.902),(.227, 0.102, 0.102),(.233, 0.00, 0.00),
                     (.287, 0.00, 0.00),(.293, 0.00, 0.00),(.346, 0.902, 0.902),(.352, 1.00, 1.00),
                     (.406, 1.00, 1.00),(.412, 1.00, 1.00),(.481, 0.00, 0.00),(.484, 0.00, 0.00),
                     (.543, 1.00, 1.00),(.546, 0.773, 0.773),(.994, 0.012, 0.012),(.997, 0.004, 0.004),
                     (1.0, 0.0, 0.0)),
             'blue': ((0.0, 0.00, 0.00),(.001, 1.00, 1.00),(.107, 0.00, 0.00),(.113, 0.498, 0.498),
                     (.173, 0.786, 0.786),(.179, 0.902, 0.902),(.227, 0.102, 0.102),(.233, 0.00, 0.00),
                     (.287, 0.00, 0.00),(.293, 0.00, 0.00),(.346, 0.00, 0.00),(.352, 0.00, 0.00),
                     (.406, 0.00, 0.00),(.412, 0.00, 0.00),(.481, 0.451, 0.451),(.484, 0.451, 0.451),
                     (.543, 1.00, 1.00),(.546, 0.773, 0.773),(.994, 0.012, 0.012),(.997, 0.004, 0.004),
                     (1.0, 0.0, 0.0))}
        
    return cpt

def abi_plot_title(fname,Band):
    # Pull Julian date from file name, convert to Gregorian date, and reformat
    julian = datetime.datetime.strptime(fname[-49:-42], '%Y%j').date()
    date = julian.strftime('%d %b %Y')
    
    Band1 = str(Band)
    
    # Create plot title
    abi_title = 'GOES-'+ fname[-53:-51] + '/ABI\n' + 'Band ' + Band1 + '\n' + fname[-42:-40] + ':' + fname[-40:-38] + ' UTC, ' + date
        
    return abi_title

def map_settings(data, extent, Band, fig):

    bmap = Basemap(resolution='h', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)

    #map_settings(ax, lon_ticks, lat_ticks, domain)

    bmap.readshapefile('/Users/anthonycrespo/Desktop/New_code_for_plotting/arg_adm1/ARG_adm1','ARG_adm1',linewidth=.5,color='black')
    bmap.readshapefile('/Users/anthonycrespo/Desktop/New_code_for_plotting/Countries_Shape/ne_10m_admin_0_countries','ne_10m_admin_0_countries',linewidth=.7,color='black')

    bmap.drawparallels(np.arange(-90.0, 90.0, 4.), linewidth=0, dashes=[4, 4], color='white', labels=[True,
    True, False, False], fmt='%g', labelstyle="+/-",xoffset=0.05, yoffset=-.50, size=15)

    bmap.drawmeridians(np.arange(0.0, 360.0, 4.), linewidth=0,dashes=[4, 4], color='white', labels=[False,
    True, False, True], fmt='%g', labelstyle="+/-",xoffset=-0.50, yoffset=0.05, size=15)

    cpt_1 = cpt(Band)

    # Makes a linear interpolation
    cpt_convert = LinearSegmentedColormap('cpt', cpt_1, 2048)

    # Plot the GOES-16 channel with the converted CPT colors
    if Band <= 6:
        bmap.imshow(data, origin='upper', cmap='Greys_r', vmin=0., vmax=1.)

    else:
        # Plot the GOES-16 channel with the converted CPT colors (you may alter the min and max to match your preference)
        bmap.imshow(data, origin='upper', cmap=cpt_convert, vmin=162, vmax=330) 
        
    # Insert the colorbar at the bottom    
    cb = bmap.colorbar(location='bottom', size = '2%', pad = '2.5%') #size adjust thickness and pad how far away from x-axis
    cb.set_label(label='Brightness temperatures (K)', size=10, weight='bold')

    
def abi_data(fname, Band, extent, resolution, variable): 
    
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

def plot_abi(save_path, file_res, extent, resolution, variable):
    
    Band = int(fname[-57:-55])
    
    #Set up figure
    DPI = file_res
    fig = plt.figure(figsize=(2000/float(DPI), 2000/float(DPI)), frameon=True, dpi=DPI, edgecolor='w')
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax = plt.axis('off')
    
    #Get data
    data = abi_data(fname, Band, extent, resolution, variable)
    
    #Format map
    map_settings(data, extent, Band, fig)
    
    #Add title
    plt.title(abi_plot_title(fname,Band), pad=10, ma='center', size=12, weight='bold')
    
    
    #Show figure
    plt.show()
    
    # Save figure as a .png file with "save_name" to specified directory
    date = datetime.datetime.strptime(fname[-49:-42], '%Y%j').date().strftime('%Y%m%d')
    save_name = 'G' + fname[-53:-51] + '_ABI_' + date + '_' + fname[-42:-38]
    fig.savefig(save_path + save_name, bbox_inches='tight', dpi=file_res)

    # Close file
    #name.close()

    # Erase plot so we can build the next one
    plt.close()


# File settings
file_path = os.getcwd() + '/'  # Directory where ABI data files are located
save_path = os.getcwd() + '/'  # Directory where figures will be saved
figures = 'single'  # Number of figures to create: 'single' (one data file) or 'multiple' (multiple data files)
file_name = 'OR_ABI-L2-CMIPF-M3C13_G16_s20180391930384_e20180391941162_c20180391941239.nc'  # For plotting a single file

# Plot settings
file_res = 150  # DPI setting for image resolution

# Mapping settings
variable = 'CMI'
extent = [-76, -58, -51, -20]
resolution = 1.

##########################################################################################################################

if __name__ == "__main__":

    if figures == 'single':
        # Open single data file and set file name
        fname = file_path + file_name
        plot_abi(save_path, file_res, extent, resolution, variable)
        print('Done!')

    if figures == 'multiple':
        # Collect all of the .nc files in specified directory
        file_list = sorted(glob.glob(file_path + '*.nc')) 
        # Loop through data files, making/saving a figure for each data file
        for fname in file_list:
            plot_abi(save_path, file_res)
        print('Done!')

