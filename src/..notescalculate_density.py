from PIL import Image
import re
from matplotlib import gridspec
import csv
import os

# Define constants
SEM_FOLDER = '../SEM/'
# DENSITY_PLOTS_FOLDER = '../Density Plots/20150216_K3_895C_L/'
# DENSITY_RESULTS = '../Density Results/Density results (K3,K8).csv'
# SAVE_ALL = False # if True, overrides keyword save in disp_results
CHIP_SIZE = (17, 17) # used for displaying the relative coordinates of an image

# Read an SEM file
# pre: a .tif file from the Keck SEM
# post: an array containing image data, the pixel size (in nm), and the stage coordinates (in mm)
def read_file(file, display = False):
    # Store the image data in an array
    im = Image.open(file)
    imarray = array(im, dtype = float)
    imarray = imarray[:,:,0]
    
    # If True, display the image
    if display:
        figure(figsize=(10,10))
        fig = imshow(imarray, cm.Greys_r)
    
    # Find  image header
    f = open(file, 'r')
    for line in f:
        if line[0:10] == 'Stage at X': y = float(re.findall("\d+.\d+", line)[0])
        if line[0:10] == 'Stage at Y': x = float(re.findall("\d+.\d+", line)[0])
        if line[0:10] == 'Pixel Size': 
            pixel_size = float(re.findall("\d+.\d+", line)[0])
            units = re.findall(".m", line)[0]
            if units == 'nm': pixel_size /= 1000
    f.close()

    return imarray, pixel_size, (x,y)


# Return a line cut of an image
# pre: the image array and the endpoints of the line cut (in pixels)
# post: a 1xN array containing values across the line cut
def line_cut_backup(imarray, p1, p2):
    line = []
    
    # Find the starting and ending coordinates relative to the two points
    imin = min(p1[0],p2[0])
    imax = max(p1[0],p2[0])
    jmin = min(p1[1],p2[1])
    jmax = max(p1[1],p2[1])
    
    # If the line is mostly horizontal, increment the x-value 
    if (jmax - jmin) > (imax - imin):
        for j in range(jmin,jmax):
            i = (imax-imin)*(j-jmin)/(jmax - jmin) + imin
            line += [imarray[i,j]]
            
    # If it's mostly vertical, increment the y-value 
    else:
        for i in range(imin,imax):
            j = (jmax-jmin)*(i-imin)/(imax - imin) + jmin
            line += [imarray[i,j]]

    return array(line)


# Return a line cut of an image
# pre: the image array and the endpoints of the line cut (in pixels)
# post: a 1xN array containing values across the line cut
def line_cut(imarray, p1_list, p2_list):
    line = []
    
    if type(p1_list) == tuple: 
        p1_list = [p1_list]
        p2_list = [p2_list]
    
    for i in range(0, len(p1_list)):
        p1 = p1_list[i]
        p2 = p2_list[i]
        
        # Find the starting and ending coordinates relative to the two points
        imin = min(p1[0],p2[0])
        imax = max(p1[0],p2[0])
        jmin = min(p1[1],p2[1])
        jmax = max(p1[1],p2[1])

        # If the line is mostly horizontal, increment the x-value 
        if (jmax - jmin) > (imax - imin):
            for j in range(jmin,jmax):
                i = (imax-imin)*(j-jmin)/(jmax - jmin) + imin
                line += [imarray[i,j]]

        # If it's mostly vertical, increment the y-value 
        else:
            for i in range(imin,imax):
                j = (jmax-jmin)*(i-imin)/(imax - imin) + jmin
                line += [imarray[i,j]]

    return array(line)


# Return the average over several parallel line cuts
# pre: the image array, the endpoints of the center line cut, 
#      and the number of parallel line cuts above and below the center line cut 
# post: a (nlines) x N array containing all of the parallel line cuts;
#       a 1xN array containing their average values
def rect_cut(imarray, p1, p2, nlines):    
    rect = [] # a 2D array for storing all line cuts
    avgs = zeros(len(line_cut(imarray,p1,p2))) # a 1D array for the average of each line cut
    
    if type(p1) == tuple: 
        p1 = [p1]
        p2 = [p2]
    
    for i in range(-nlines/2+1,nlines/2+1):
        p1_shift = list(p1)
        p2_shift = list(p2)
        for j in range(0,len(p1)):
            p1_shift[j] = (p1[j][0]+i,p1[j][1])
            p2_shift[j] = (p2[j][0]+i,p2[j][1])
        
        cut = line_cut(imarray,p1_shift,p2_shift)
        
        # Stack line cuts vertically
        if len(rect) == 0:
            rect = cut
        else:
            rect = vstack((rect, cut))
            
        # Add line cuts together to find average
        avgs += cut
#     # Loop through parallel line cuts
#     for i in range(-dx,dx+1):
#         if abs(p2[1] - p1[1]) > abs(p2[0] - p1[0]):
#             cut = line_cut(imarray,(p1[0]+i,p1[1]),(p2[0]+i,p2[1]))
#         else:
#             cut = line_cut(imarray,(p1[0],p1[1]+i),(p2[0],p2[1]+i))
        
#         # Stack line cuts vertically
#         if len(rect) == 0:
#             rect = cut
#         else:
#             rect = vstack((rect, cut))
        
#         # Add line cuts together to find average
#         avgs += cut
    
    if nlines > 0: avgs = avgs / nlines
        
    return array(rect), avgs


# Locate the tubes in a line cut
# pre: the 2nd derivative of the line cut and a threshold to control the sensitivity (usually 20-30)
# post: a list containing the tube locations along the line cut
def find_tubes(diff2, threshold):
    tubes = []
    peak = False
    
    for i in range(0, len(diff2)):
        # Peaks must be larger than the threshold
        if diff2[i] >= threshold:
            if not peak:
                peak = True
                tubes += [i]
            else:
                # Record the location of the peak
                if diff2[i] > diff2[tubes[-1]]: tubes[-1] = i
        else:
            peak = False
            
    return tubes


# Convert the position of an image from absolute to relative coordinates
# pre: the stage position of an image and up to 2 corners
# post: if no corners given, returns (nan, nan)
#       if 1 corner is given, assumes that corners are aligned along the axes
#       if 2 corners are given, determines their rotation and calculates the image coordinates accordingly
def abs_to_rel_coords(image_pos, corner_pos, chip_size):
    if type(corner_pos) == tuple and len(corner_pos) == 2:
        return (abs(image_pos[0] - corner_pos[0]), abs(image_pos[1] - corner_pos[1])) 
    elif type(corner_pos) == list and len(corner_pos) == 2:
        rot = arctan((corner_pos[1][1]-corner_pos[0][1]) / (corner_pos[1][0]-corner_pos[0][0]))
        
        dy = corner_pos[1][1]-image_pos[1]
        dx = corner_pos[1][0]-image_pos[0]
        theta = arctan(dy / dx)
        dist = sqrt(dx**2 + dy**2)
        
        pos = (abs(dist * sin(theta - rot)), abs(dist * cos(theta - rot)))
        
        if chip_size[0] < 0:
            pos = (abs(chip_size[0]) - pos[0], pos[1])
        
        if chip_size[1] < 0:
            pos = (pos[0], abs(chip_size[1]) - pos[1])
        
        return pos
    else:
        return (nan, nan)


# Display the starting image, a zoom-in of the line cut, the location of the tubes, 
#  and a histogram of their relative spacing
# pre: an SEM file, the endpoints of the line cut, the number of lines to average, and a threshold
# post: display the results in a single figure and saves them to a file if keyword save is set
def disp_results(file, p1, p2, pixel_size = 0, threshold = 20, nlines = 20, 
                 corner_pos = None, chip_size = CHIP_SIZE, save = False):
    
    if type(p1) == tuple: 
        p1 = [p1]
        p2 = [p2]
    
    # Read the image file
    imarray, pixel_size, pos = read_file(SEM_FOLDER + file, display = False)
    
    # Convert its position to relative coordinates
    pos = abs_to_rel_coords(pos, corner_pos, chip_size)
    
    # Take a line cut and look for tubes
    rect, avgs = rect_cut(imarray, p1, p2, nlines)
    diff2 = -diff(diff(avgs))
    tubes = find_tubes(diff2, threshold)
    
    fig = figure(figsize = (17,12))
    
    if pos == (nan, nan):
        gs1 = gridspec.GridSpec(1, 1)
        gs1.update(right = 0.55)
    else:
        gs1 = gridspec.GridSpec(2, 1, height_ratios = [1.5,1])
        gs1.update(right = 0.55, hspace = 0.1)

    gs2 = gridspec.GridSpec(2, 1, height_ratios = [1, 5])
    gs2.update(left = 0.6, bottom = 0.50, hspace = 0.05)
    
    gs3 = gridspec.GridSpec(1, 1)
    gs3.update(left = 0.6, top = 0.43)
    
    # gs = gridspec.GridSpec(3, 2, 
    #         width_ratios=[1.5, 1], height_ratios=[1, 5, 5]) 
    # gs.update(wspace = 0.1, hspace = 0.2)
    
    # Show the starting image
    ax0 = subplot(gs1[0,0])
    ax0.imshow(imarray, cm.Greys_r)

    # Outline the region used for the line cut
    dx = nlines / 2
    for i in range(0, len(p1)):
        y = [p1[i][0]-dx, p1[i][0]+dx, p2[i][0]+dx, p2[i][0]-dx, p1[i][0]-dx]
        x = [p1[i][1], p1[i][1], p2[i][1], p2[i][1], p1[i][1]]
        
        ax0.plot(x, y, 'b-', lw = 2)
#     if abs(p2[1] - p1[1]) > abs(p2[0] - p1[0]):
#         y = [p1[0]-dx, p1[0]+dx, p2[0]+dx, p2[0]-dx, p1[0]-dx]
#         x = [p1[1], p1[1], p2[1], p2[1], p1[1]]
#     else:
#         x = [p1[1]-dx, p1[1]+dx, p2[1]+dx, p2[1]-dx, p1[1]-dx]
#         y = [p1[0], p1[0], p2[0], p2[0], p1[0]]
#     ax0.plot(x, y, 'b-', lw = 2)
    ax0.margins(0,0)
    
    if len(tubes) >= 10:
        mean_spacing = mean(diff(tubes))*pixel_size
        std_spacing = std(diff(tubes))*pixel_size
    elif len(tubes) > 0:
        mean_spacing = len(avgs)*pixel_size / len(tubes)
        std_spacing = nan
    else:
        mean_spacing = nan;
        std_spacing = nan;
    
    # t = ax0.set_title(os.path.basename(file), fontsize = 16)
    # t.set_y(1.02)
    # if pixel_size > 0 and len(tubes) > 0: 
    #     title = 'Avg. Spacing: ' + "%.2f" % spacing + ' um'
    # else:
    #     title = 'No. Tubes: ' + str(len(tubes))
    # t = ax0.set_title(title, fontsize = 14)
    # t.set_y(1.01)
    
    # Show a zoomed in view of the line cut
    ax1 = subplot(gs2[0,0])
    ax1.imshow(rect, cm.Greys_r)
    ax1.get_xaxis().set_visible(False)
    ax1.get_yaxis().set_visible(False)
    
    # Mark the tubes with red x's
    ax1.plot([x+1 for x in tubes],[nlines/2]*len(tubes),'rx')
    ax1.axis('off')
    ax1.set_aspect('auto')
    ax1.margins(0,0)
    
    # Plot the 2nd derivative of the line cut, and also mark the tubes with x's
    ax2 = subplot(gs2[1,0])
    ax2.plot(linspace(0,len(diff2)*pixel_size,len(diff2)),diff2)
    ax2.plot([x*pixel_size for x in tubes],diff2[tubes],'rx')
    ax2.plot(linspace(0,len(diff2)*pixel_size,len(diff2)),[threshold]*len(diff2),'k--')
    ax2.set_xlim([0, len(diff2)*pixel_size])
    ax2.set_xlabel('Position (um)', fontsize = 14)
      
    if pos != (nan, nan):
        ax3 = subplot(gs1[1,0])
        ax3.set_aspect(chip_size[0]/chip_size[1]/min(chip_size[0]/chip_size[1],2.5))
        ax3.plot(pos[0],pos[1],'sk', ms = 10)
        
        ax3.annotate('(%0.1f, %0.1f)' % pos, xy = pos, xytext=(0, 15), 
                     ha = 'center', textcoords='offset points', fontsize = 14)
        
        ax3.set_xlim([0,abs(chip_size[0])])
        ax3.set_ylim([0,abs(chip_size[1])])
        # ax3.get_xaxis().set_visible(False)
        # ax3.get_yaxis().set_visible(False)
        ax3.set_xlabel('Image Location (mm)', fontsize = 14)
        
    
    ax4 = subplot(gs3[0,0])
    if len(tubes) > 2:
        ax4.hist([x*pixel_size for x in diff(tubes)],bins = 20, color = 'g')
        ax4.set_xlabel('Tube Spacing (um)', fontsize = 14)
    ax4.annotate('avg = %0.1f' % mean_spacing, xy = (1,1), xytext=(-10, -20),
                    xycoords = "axes fraction", ha = 'right', 
                    textcoords='offset points', fontsize = 14)
    ax4.annotate('std = %0.1f' % std_spacing, xy = (1,1), xytext=(-10, -42),
                    xycoords = "axes fraction", ha = 'right', 
                    textcoords='offset points', fontsize = 14)
    
#     if save or SAVE_ALL:
#         filename = os.path.splitext(os.path.split(file)[1])[0]
#         savefig(DENSITY_PLOTS_FOLDER + filename + '_(%0.f_%0.f)_' % pos + '.png',
#             bbox_inches = 'tight', pad_inches= 0.2)
        
#         chip = ''; temp = ''
#         for s in filename.upper().split('_'):
#             for x in ['I','J','K','M','N']:
#                 if x in s: chip = s
#             if 'C' in s: temp = s.replace('C','')
#             if temp == '' and s.isdigit() and int(s) > 800 and int(s) < 1100: temp = s
        
#         csvfile = open(DENSITY_RESULTS,'a')
#         writer = csv.writer(csvfile, lineterminator='\n', delimiter = ';')
#         writer.writerow([filename, chip, temp, '%.2f' % mean_spacing, '%.2f' % std_spacing, 
#                          '(%0.1f, %0.1f)' % pos, str(p1), str(p2), str(nlines), ''])
    
    # gs.tight_layout(fig)
    
    show()

