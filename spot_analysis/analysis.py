import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
import cv2
import json
from sklearn.mixture import GaussianMixture

import circles
import html_utils

per_base_normalization = False
per_spot_normalization = False
per_cycle_normalization = False
verbose = 1

argcc = 1
argc = len(sys.argv)
while argcc < argc:
    if sys.argv[argcc] == '--per_base_normalization':
        per_base_normalization = True
    if sys.argv[argcc] == '--per_spot_normalization':
        per_spot_normalization = True
    if sys.argv[argcc] == '--per_cycle_normalization':
        per_cycle_normalization = True
    if sys.argv[argcc] == '-v':
        verbose += 1
    argcc += 1


# create file names
filenames = glob.glob('*.tif')
if len(filenames) < 1:
    filenames = glob.glob('*.png')
filenames = sorted(filenames)
num = len(filenames)
cycles = []
spot_file = -1
for i in range(num):
    # first UV image is our spot file
    if spot_file == -1 and '_365_' in filenames[i]:
        spot_file = i
    # 4 image files preceeding the UV image are the intensity captures at the end of a cycle
    if '_365_' in filenames[i]:
        cycles.append(i-4)
cycles.append(num-4)

num_cycles = len(cycles)
if verbose > 0:
    print('found %d cycles' % num_cycles)

# find spots
#   load spot image
#   convert to grayscale
#   find spots
num_spots = 0
spot_image = cv2.imread(filenames[spot_file])
spot_image_gray = cv2.cvtColor(spot_image, cv2.COLOR_BGR2GRAY)

spots = circles.findCircles2(spot_image_gray, minr=30)
num_spots = len(spots)

def get_intensity(image, x, y, r):
    #print('finding intensity at (%d,%d) rad %d' % (x, y, r))

    hist = np.zeros(256) # the intensities
    circle = np.zeros((r*2-1,r*2-1)).astype(np.uint8)
    cv2.circle(circle, (r-1,r-1), r, 1, -1)
    count = 0
    for cy in range(r*2-1):
        for cx in range(r*2-1):
            if circle[cy,cx] == 1:
                val = image[y-r+cy+1, x-r+cx+1]
                hist[val] += 1
                count += 1
    #print('points considered: %d' % count)

    med_count = 0
    target = count / 2
    median_intensity = 0
    for i in range(256):
        med_count += hist[i]
        if med_count >= target:
            median_intensity = i
            break

    return median_intensity

# generate json intensitiy information from each cycle
intensity = np.zeros((num_cycles, num_spots, 4))
for cycle in range(num_cycles):
    if verbose > 0:
        print('processing cycle: %d' % cycle)
    for base in range(4):
        image = cv2.imread(filenames[cycles[cycle]+base])
        image_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
        spot_num = 0
        for spot in spots:
            x = int(spot[0])
            y = int(spot[1])
            r = int(spot[2])
            if (x - r) < 0 or (y - r) < 0:
                continue
            intensity[cycle,spot_num,base] = get_intensity(image_gray, x, y, r)
            spot_num += 1

# generate a 1D array of all min & max values (will histogram this later)
minMaxList = np.concatenate([np.amin(intensity[0,:,:], axis=1), np.amax(intensity[0,:,:], axis=1)])

# find the median dark and median bright intensity across all spots, per base
# then use the mean of each distribution to normalize
if per_cycle_normalization:
    for base in range(4):
        baseIntensities = intensity[0,:,base].reshape(-1, 1)
        if verbose > 1:
            print('base intensities:\n%s\nshape: %s' % (baseIntensities, baseIntensities.shape))
        # this should be a bi-modal distribution, but we don't know how many spots are the bright and how many are the low
        gm = GaussianMixture(n_components=2, random_state=0).fit(baseIntensities)
        if verbose > 1:
            print(f'mu1={gm.means_[0]}, mu2={gm.means_[1]}')
        minIntensity = gm.means_[0]
        maxIntensity = gm.means_[1]
        if maxIntensity < minIntensity:
            temp = minIntensity
            minIntensity = maxIntensity
            maxIntensity = temp
        intensityDelta = (maxIntensity - minIntensity)
        intensity[:, :, base] -= minIntensity
        if intensityDelta > 0:
            intensity[:, :, base] /= intensityDelta

if per_spot_normalization:
    for spot in range(num_spots):
        minIntensity = np.min(intensity[:, spot, :])
        maxIntensity = np.max(intensity[:, spot, :])
        intensity[:, spot, :] -= minIntensity
        intensityDelta = (maxIntensity - minIntensity)
        if intensityDelta > 0:
            intensity[:, spot, :] /= intensityDelta

if per_base_normalization:
    for spot in range(num_spots):
        for base in range(4):
            minIntensity = np.min(intensity[:, spot, base])
            maxIntensity = np.max(intensity[:, spot, base])
            intensityDelta = (maxIntensity - minIntensity)
            intensity[:, spot, base] -= minIntensity
            if intensityDelta > 0:
                intensity[:, spot, base] /= intensityDelta


with open('run.data', 'w') as f:
    for spot in range(num_spots):
        for cycle in range(num_cycles):
            f.write('%d,%d,%.3f,%.3f,%.3f,%.3f\n' % (spot, cycle,
                intensity[cycle, spot, 0],
                intensity[cycle, spot, 1],
                intensity[cycle, spot, 2],
                intensity[cycle, spot, 3]))

#
# Generate plots
#

'''
fig, axs = plt.subplots(4, 1)
for base in range(4):
    baseIntensities = intensity[0,:,base]
    axs[base].hist(baseIntensities, 21)
    axs[base].set_title('base %d intensites' % base)
'''


# this is how to vary the plot size
#fig.set_size_inches(12.8,10.24) # generates a 1280x1024 image

#
# Generate HTML Report
#

html = html_utils.HTMLUtils('report.html')
html.add_header('Analysis Report')

# show the grayscale image and the one with the spot overlay
fig, ax = plt.subplots()
ax.imshow(spot_image_gray, cmap='gray')
plt.savefig('cycle0UV.jpg')
plt.close()

fig, ax = plt.subplots()
ax.imshow(spot_image_gray, cmap='gray')
for spotnum, spot in enumerate(spots):
    circle = plt.Circle((spot[0], spot[1]), spot[2], color='r')
    ax.add_patch(circle)
    plt.text(spot[0]-20, spot[1]+10, str(spotnum), color='white')
plt.savefig('spots.jpg')
plt.close()

html.add_images([
    {'img':'cycle0UV.jpg', 'title': 'Cycle 0 UV Image', 'subtitle': '(gray scale of cycle 0 UV image)'},
    {'img':'spots.jpg', 'title': 'Cycle 0 spots found'}])

# show a histogram of cycle zero intensities
fig = plt.figure('min max hist')
plt.hist(minMaxList, 21)
plt.savefig('min_max_histogram.jpg')
plt.close()
html.add_images([{'img': 'min_max_histogram.jpg', 'title': 'cycle 0 intensity distribution'}])

# show per-spot intensity bargraphs

html.start_div()
html.add_text('<h2>Spot Information</h2>')

base_colors = ['green', 'yellow', 'blue', 'red']
#fig, axs = plt.subplots(int(num_spots/4+0.75), 4) # row,col
for spot in range(num_spots):
    fig, ax = plt.subplots()
    if per_base_normalization or per_spot_normalization or per_cycle_normalization:
        ax.set_ylim([0, 1.0])
    else:
        ax.set_ylim([0, 200])
    for cycle in range(num_cycles):
        for base in range(4):
            ax.bar(cycle + base*0.1, intensity[cycle, spot, base], color = base_colors[base], width = 0.1)
    fig.set_size_inches(3.2,2.0)
    plt.savefig('spot_%d.jpg' % spot)
    plt.close()

for spot in range(0,num_spots,4):
    spot_info = []
    for k in range(4):
        spot_k = spot+k
        if spot_k < num_spots:
            spot_info.append({'img': 'spot_%d.jpg' % spot_k, 'title': 'Spot %d' % spot_k})
    html.add_images(spot_info)
    html.add_text('<br>\n') # add a little space between rows
html.end_div()

html.close()

