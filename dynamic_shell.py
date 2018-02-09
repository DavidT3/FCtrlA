from __future__ import division
from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import pyplot as plt
import numpy as np
from numpy import cos, sin, linspace, pi

low_energy = fits.open('0.50-2.00keV_img.fits')
high_energy = fits.open('2.00-10.0keV_img.fits')
low_data = low_energy[0].data
high_data = high_energy[0].data
data = low_data + high_data

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


"""plt.imshow(low_data, origin='lower')
plt.plot(x_pix, y_pix, 'x', color='black')
plt.plot(xcir_pix, ycir_pix, color='red')
plt.show()"""


def annular_mask(index, radius, array):
    a, b = index
    nx, ny = array.shape
    y, x = np.ogrid[-a:nx-a, -b:ny-b]
    mask = x*x + y*y >= radius*radius
    return mask


r = round(abs(x_pix - wcs.wcs_world2pix(np.array([[x_fk5 + r500_fk5, y_fk5]], np.float), 1)[0][0]))
data[annular_mask([y_pix, x_pix], r, data)] = 0

sums = []
rf = []
r_frac = r/50
for i in range(1, 51):
    data_copy = np.copy(data)
    rf.append(i * r_frac)
    data_copy[annular_mask([y_pix, x_pix], rf[i-1], data_copy)] = 0
    sums.append(np.sum(data_copy))

rings = []
for i in range(0, len(sums)):
    if i == 0:
        rings.append(sums[i])
    else:
        rings.append(sums[i] - sums[i - 1])


"""plt.plot(radii_pix, final_cnts, 'x', color='black')
plt.axhline(y=target, color='red')
plt.show()"""

# This should work, just now need to loop it half a million times to find the best combo of N and target
# Will have to use chi square or something similar to assess which is the best.
# Just give every count value the same error, or might as well just go with residuals

perc_list = np.arange(0.05, 0.3, 0.01)
n_list = []
f_list = []
resid_list = []
p_list = []
for N in range(3, 9):
    for perc in perc_list:
        ann_tot = 0  # Finds size of first annulus
        target = sums[-1] * perc
        count = 0
        while ann_tot < target:
            ann_tot += rings[count]
            count += 1
        rad_1 = count * r_frac

        F = (r / rad_1) ** (1 / (N - 1))
        radii = [int(rad_1 / r_frac)]  # In ring units (not actual radius)
        for i in range(2, N + 1):
            radii.append(int(round(radii[0] * F ** (i - 1))))

        final_cnts = []
        for i in range(0, len(radii)):
            if i == 0:
                final_cnts.append(sums[radii[i] - 1])
            else:
                final_cnts.append(sums[radii[i] - 1] - sums[radii[i - 1] - 1])

        resid = sum([(cnt - target)**2 for cnt in final_cnts])
        n_list.append(N)
        f_list.append(F)
        p_list.append(perc)
        resid_list.append(resid)

print('The best N is ' + str(n_list[resid_list.index(min(resid_list))]))
print('The best F is ' + str(f_list[resid_list.index(min(resid_list))]))
print('The percentage target was ' + str(p_list[resid_list.index(min(resid_list))]))


# radii_pix = [x * r_frac for x in radii]



