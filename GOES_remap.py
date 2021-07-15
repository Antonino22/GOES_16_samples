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
from osgeo import gdal
import glob
import pandas as pd
import datetime
import shutil
import warnings

#Ignore deprecated warnings
warnings.filterwarnings("ignore")


all_reports = {'lat':[-31.470844, -31.422844, -31.431458, -31.421580, -31.435851, -32.201393, -31.4201, 
                      -31.434703, -31.492307, -31.411644, -31.335028],
               'lon':[-64.541380, -64.497452, -64.499612, -64.500763, -64.456595, -64.452606, -64.1888, 
                      -64.504108, -64.544306, -64.506675, -64.279508],
               'cities': ["Icho Cruz","San Lorenzo","Tupungato", "Calle Tokio", "San Nicolas", 
                          "Villa del Dique", "Cordoba", "Villa Carlos Paz", "San Antonion de Redondo", 
                          "Villa Carlos Paz N", "Villa Allende"],
               'date': ['Feb 8 18', 'Feb 8 18', 'Feb 8 18', 'Feb 8 18', 'Feb 8 18', 'Dec 13 18', "Dec 13 18", 
                        "Oct 24 20", "Oct 24 20", "Oct 24 20", "Oct 24 20"],
               'time': ["18:50", "19:20", "19:20", "19:30", "19:45","02:20", "03:00","18:20", 
                        "18:32", "18:35", "19:27"]}
report = pd.DataFrame (all_reports, columns = ['lat','lon','cities','date','time'])

def ext(extent1):
    extent = extent1
    if extent == 'Nation':
        zoom = [-76, -58, -51, -20]
        return zoom
    elif extent == 'Dec13':
        zoom = [-68, -35.5, -60, -31] #Hail event for Dec 13
        return zoom
    elif extent == 'Feb8-NC':
        zoom = [-65.25, -32.25, -63., -31.15]
        return zoom  
    elif extent == 'Feb8':
        zoom = [-66., -34.5, -62.25, -31.]
        return zoom      
    elif extent == 'Oct24':
        zoom = [-64.75, -32., -62., -31.]
        return zoom      
    else:
        pass    


def sat_variables(path, dis):
    nc = Dataset(path, 'r')

    # Calculate the image extent 
    H = nc.variables['goes_imager_projection'].perspective_point_height 
    x1 = nc.variables['x_image_bounds'][0] * H
    x2 = nc.variables['x_image_bounds'][1] * H
    y1 = nc.variables['y_image_bounds'][1] * H
    y2 = nc.variables['y_image_bounds'][0] * H
    
    # Read the central longitude
    longitude = nc.variables['goes_imager_projection'].longitude_of_projection_origin
         
    if dis == '4':
        # Read the semi major axis
        a = nc.variables['goes_imager_projection'].semi_major_axis + 4572

        # Read the semi minor axis
        b = nc.variables['goes_imager_projection'].semi_minor_axis + 9144
    elif dis == '9':
        # Read the semi major axis
        a = nc.variables['goes_imager_projection'].semi_major_axis + 9144

        # Read the semi minor axis
        b = nc.variables['goes_imager_projection'].semi_minor_axis + 9144
    elif dis == '15':
        # Read the semi major axis
        a = nc.variables['goes_imager_projection'].semi_major_axis + 15240

        # Read the semi minor axis
        b = nc.variables['goes_imager_projection'].semi_minor_axis + 15240
    else:
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
    #abi_title =  date + ',  ' + fname[-42:-40] + ':' + fname[-40:-38] + ' UTC' + '             ' + 'GOES-' + fname[-53:-51] + '/ABI' + ' ' + 'Band ' + Band1 
    Name = 'GOES-' + fname[-53:-51] + '/ABI' + ' ' + 'Band ' + Band1
    Date = date + ',  ' + fname[-42:-40] + ':' + fname[-40:-38] + ' UTC' 

    plt.title(' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + Name + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' '  + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' +
    ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' + ' ' 
    + Date + ' ' + ' ' + ' ' + ' ', y=1, size=20, weight='bold')

    #return abi_title

def map_settings(data, extent, extent1, Band, fig):

    bmap = Basemap(resolution='h', llcrnrlon=extent[0], llcrnrlat=extent[1], urcrnrlon=extent[2], urcrnrlat=extent[3], epsg=4326)

    #map_settings(ax, lon_ticks, lat_ticks, domain)

    bmap.readshapefile('/Users/anthonycrespo/Desktop/New_code_for_plotting/arg_adm1/ARG_adm1','ARG_adm1',linewidth=.5,color='black')
    bmap.readshapefile('/Users/anthonycrespo/Desktop/New_code_for_plotting/Countries_Shape/ne_10m_admin_0_countries','ne_10m_admin_0_countries',linewidth=.7,color='black')

    if extent1=='Nation':
        bmap.drawparallels(np.arange(-90.0, 90.0, 2.), linewidth=.5, dashes=[4, 4], color='white', labels=[True,
        True, False, False], fmt='%g', labelstyle="+/-",xoffset=0.05, yoffset=-.50, size=15)
        bmap.drawmeridians(np.arange(0.0, 360.0, 2.), linewidth=1,dashes=[4, 4], color='white', labels=[False,
        True, False, True], fmt='%g', labelstyle="+/-",xoffset=-0.50, yoffset=0.05, size=15)
    elif extent1=='Dec13':
        bmap.drawparallels(np.arange(-90.0, 90.0, .5), linewidth=0.3,dashes=[4, 4], color='white', labels=[True,
        True, False, False], fmt='%g', labelstyle="+/-", xoffset=0.05, yoffset=-.50, size=15)
        bmap.drawmeridians(np.arange(0.0, 360.0, .5), linewidth=1,dashes=[4, 4], color='white', labels=[False,
        True, False, True], fmt='%g', labelstyle="+/-", xoffset=-0.50, yoffset=0.05, size=15) 
    elif extent1 == 'Feb8-NC':
        bmap.drawparallels(np.arange(-90.0, 90.0, .15), linewidth=0.3,dashes=[4, 4], color='white', labels=[True,
        True, False, False], fmt='%g', labelstyle="+/-", xoffset=0.05, yoffset=-.50, size=15)
        bmap.drawmeridians(np.arange(0.0, 360.0, .15), linewidth=1,dashes=[4, 4], color='white', labels=[False,
        True, False, True], fmt='%g', labelstyle="+/-", xoffset=-0.50, yoffset=0.05, size=15)
    elif extent1 == 'Feb8':
        bmap.drawparallels(np.arange(-90.0, 90.0, .5), linewidth=0.3, dashes=[4, 4], color='white', 
        labels=[True,True, False, False], fmt='%g', labelstyle="+/-", xoffset=0.05, yoffset=-.50, size=15)
        bmap.drawmeridians(np.arange(0.0, 360.0, .5), linewidth=1,dashes=[4, 4], color='white', 
        labels=[False,False, False, True], fmt='%g', labelstyle="+/-",xoffset=-0.50, yoffset=0.05, size=15)   
    elif extent1 == 'Oct24':
        bmap.drawparallels(np.arange(-90.0, 90.0, .15), linewidth=0.3,dashes=[4, 4], color='white', labels=[True,
        True, False, False], fmt='%g', labelstyle="+/-", xoffset=0.05, yoffset=-.50, size=15)
        bmap.drawmeridians(np.arange(0.0, 360.0, .2), linewidth=1,dashes=[4, 4], color='white', labels=[False,
        True, False, True], fmt='%g', labelstyle="+/-", xoffset=-0.50, yoffset=0.05, size=15)        
    else:
        pass       

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
    if extent1 == 'Nation':    
        cb = bmap.colorbar(location='bottom', size = '2%', pad = '2.5%') #size adjust thickness and pad how far away from x-axis
        cb.set_label(label='Brightness temperatures (K)', size=10, weight='bold')
    elif extent1 == 'Dec13':  
        cb = bmap.colorbar(location='bottom', size = '2%', pad = '5%') #size adjust thickness and pad how far away from x-axis
        cb.set_label(label='Brightness temperatures (K)', size=10, weight='bold')
    elif extent1 == 'Feb8-NC':
        cb = bmap.colorbar(location='bottom', size = '2%', pad = '10%') #size adjust thickness and pad how far away from x-axis
        cb.set_label(label='Brightness temperatures (K)', size=10, weight='bold')
    elif extent1 == 'Feb8':
        cb = bmap.colorbar(location='bottom', size = '3%', pad = '4.5%') #size adjust thickness and pad how far away from x-axis
        cb.set_label(label='Brightness temperatures (K)', size=10, weight='bold')     
    elif extent1 == 'Oct24':
        cb = bmap.colorbar(location='bottom', size = '2%', pad = '11%') #size adjust thickness and pad how far away from x-axis
        cb.set_label(label='Brightness temperatures (K)', size=10, weight='bold')        
    else:
        pass    

    date = datetime.datetime.strptime(fname[-49:-42], '%Y%j').date().strftime('%Y%m%d') 
                     

    
def abi_data(fname, Band, extent, resolution, variable, dis): 
    
    H, a, b, f, longitude, x1, y1, x2, y2 = sat_variables(fname, dis)

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
     
    return data, grid      

def plot_abi(save_path, file_res, extent1, resolution, variable, report, image, geotif, dis):
    
    Band = int(fname[-57:-55])

    extent = ext(extent1)

    #Set up figure
    DPI = file_res
    fig = plt.figure(figsize=(2000/float(DPI), 2000/float(DPI)), frameon=True, dpi=DPI, edgecolor='w')
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)
    ax = plt.axis('off')

    #Get data
    data, grid = abi_data(fname, Band, extent, resolution, variable, dis)
    
    #Format map
    map_settings(data, extent, extent1, Band, fig)
    
    #Add title
    #plt.title(abi_plot_title(fname,Band), pad=10, ma='center', size=12, weight='bold')
    abi_plot_title(fname,Band)
    
    # Save figure as a .png file with "save_name" to specified directory
    date = datetime.datetime.strptime(fname[-49:-42], '%Y%j').date().strftime('%Y%m%d')
    save_name = 'G' + fname[-53:-51] + '_ABI_' + date + '_' + fname[-42:-38]

    #Show figure
    #plt.show()

    if image == 'yes':
        plt.savefig(save_path + save_name, bbox_inches='tight', dpi=file_res)
    else:
        pass

    # Export the result to GeoTIFF
    if geotif == 'yes':
        driver = gdal.GetDriverByName('GTiff')
        driver.CreateCopy(save_path + save_name + '.tif', grid, 0)
    else:
        pass

    source = fname + '.aux.xml'
    destination =  fname[0:45] + 'z_aux_files' + fname[59:200]+ '.aux.xml'

    # Move a file from the directory d1 to d2
    shutil.move(source, destination)

    # Delete aux.xml
    os.remove(destination) 

    # Close the plot window
    plt.close()


# File settings
#file_path = os.getcwd() + '/'  # Directory where ABI data files are located
file_path = os.getcwd() + '/' + '02_08_2018' + '/' + 'C11' + '/'
save_path = os.getcwd() + '/'  # Directory where figures will be saved
figures = 'single'  # Number of figures to create: 'single' (one data file) or 'multiple' (multiple data files)
file_name = 'OR_ABI-L2-CMIPF-M3C11_G16_s20180391900384_e20180391911151_c20180391911240.nc'  # For plotting a single file

# Plot settings
file_res = 150  # DPI setting for image resolution

# Mapping settings
image = 'yes'
geotif = 'yes'
variable = 'CMI'
extent1 = 'Feb8'
dis = '15'
resolution = 1.

##########################################################################################################################

if __name__ == "__main__":

    if figures == 'single':
        # Open single data file and set file name
        fname = file_path + file_name
        plot_abi(save_path, file_res, extent1, resolution, variable, report, image, geotif, dis)
        print('Done!')

    if figures == 'multiple':
        # Collect all of the .nc files in specified directory
        file_list = sorted(glob.glob(file_path + '*.nc')) 
        # Loop through data files, making/saving a figure for each data file
        for fname in file_list:
            plot_abi(save_path, file_res, extent1, resolution, variable, report, image, geotif. dis)
        print('Done!')

