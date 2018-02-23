from __future__ import division
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
import numpy as np
from numpy import cos, sin, linspace, pi, shape
import sys

low_energy = fits.open('A907_0.5-2.0.fits')
low_data = low_energy[0].data
data = low_data   # Kathy just said to use the low energy counts
# Too much background in higher energy?

wcs = WCS(low_energy[0].header)

x_fk5 = 149.5916
y_fk5 = -11.0642
r500_fk5 = 0.1065690271986475


angle = linspace(0, 2*pi, 360)
xcir_fk5 = []
ycir_fk5 = []
for i in angle:
    xcir_fk5.append(x_fk5 + (r500_fk5 * cos(i)))
    ycir_fk5.append(y_fk5 + (r500_fk5 * sin(i)))

xcir_pix = []
ycir_pix = []
for i in range(0, len(xcir_fk5)):
    plc_hld = wcs.wcs_world2pix(np.array([[xcir_fk5[i], ycir_fk5[i]]], np.float), 1)
    xcir_pix.append(plc_hld[0][0])
    ycir_pix.append(plc_hld[0][1])


centre = wcs.wcs_world2pix(np.array([[x_fk5, y_fk5]], np.float), 1)
x_pix = centre[0][0]
y_pix = centre[0][1]


def annular_mask(index, radius, array):
    a, b = index
    nx, ny = array.shape
    y, x = np.ogrid[-a:nx-a, -b:ny-b]
    mask = x*x + y*y >= radius*radius
    return mask


r = round(abs(x_pix - wcs.wcs_world2pix(np.array([[x_fk5 + r500_fk5, y_fk5]], np.float), 1)[0][0]))
ylim, xlim = shape(data)
if (x_pix + r) > xlim:
    print('CRITICAL ERROR: R500 OUTSIDE OF IMAGE ON THE RIGHT SIDE')
    sys.exit()
elif(x_pix - r) < 0:
    print('CRITICAL ERROR: R500 OUTSIDE OF IMAGE ON THE LEFT SIDE')
    sys.exit()
elif (y_pix + r) > ylim:
    print('CRITICAL ERROR: R500 OUTSIDE OF IMAGE ON THE TOP')
    sys.exit()
elif(y_pix - r) < 0:
    print('CRITICAL ERROR: R500 OUTSIDE OF IMAGE ON THE BOTTOM')
    sys.exit()

data[annular_mask([y_pix, x_pix], r, data)] = 0
sums = []
rf = []
how_many = 200
r_frac = r/how_many
for i in range(1, 201):
    data_copy = np.copy(data)
    rf.append(i * r_frac)
    data_copy[annular_mask([y_pix, x_pix], rf[i-1], data_copy)] = 0
    sums.append(np.sum(data_copy))

if sums[-1] < 2000:
    print('ERROR: THIS CLUSTER HAS LESS THAN 2000 counts')
    sys.exit()

rings = []
for i in range(0, len(sums)):
    if i == 0:
        rings.append(sums[i])
    else:
        rings.append(sums[i] - sums[i - 1])


ann_tot = 0  # Finds outer radius of N-1th annulus
N = 8
target = 2000
count = 1
annuli_id = N - 1
while ann_tot < target:
    ann_tot += rings[-count]
    count += 1
rad_minus_1 = how_many - (count - 1)  # Small ring no

F = (how_many / rad_minus_1) ** (1 / (N - annuli_id))  # how_many is outer_rad in number of small rings
radii = [int(float(how_many)/F**(N-i-1)) - 1 for i in range(N)]

final_cnts = []
for i in range(0, len(radii)):
    if i == 0:
        final_cnts.append(sums[radii[i]])
    else:
        final_cnts.append(sums[radii[i]] - sums[radii[i - 1]])
pix_rad = [i*r_frac for i in radii]

thing = [round(i) for i in pix_rad]
width = [thing[0]]
centre = [width[0]/2]
for i in range(1, len(thing)):
    width.append(thing[i] - thing[i - 1])
    centre.append(width[i]/2 + thing[i - 1])


plt.bar(centre, final_cnts, width=width)

#plt.plot(pix_rad, final_cnts, color='black')
for i in range(0, N):
    plt.axvline(round(pix_rad[i]), color='red', ls='dashed')
plt.xlabel('Annular Radius [PIX]')
plt.ylabel('Annular Count')
plt.show()


# This now just needs writing into FCtrlA
# Maybe in the FCtrlA rewrite I am planning
# No doubt this code can be made better anyhow

