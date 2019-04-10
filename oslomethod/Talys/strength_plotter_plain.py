from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
import seaborn as sns
import matplotlib.ticker as ticker

Bn=6.534000 # MeV

sns.set()
# sns.set_context("paper")
sns.set_context("paper",
                font_scale=1.4,
                rc={'legend.loc' : 'best', 'figure.figsize':(5.5,5.5)})
#sns.set(rc={'figure.figsize':(5.5,6.5)})
sns.set_style("ticks", { 'axes.grid': False})
# plt.rcParams['legend.loc'] = 'best'

# Cross section factor for formula strength = xsec*factor/Egamma
xsec_factor = 8.68e-8

def readGSF(folder):

	# import fgteo.rsg to get the calibration
	folded1, cal, E_array, tmp = read_mama_2D(folder+"/fgteo.rsg")

	# Constants for energy binning
	if(cal["a0x"] != cal["a0y"]
		or cal["a1x"] != cal["a1y"]
		or cal["a2x"] != cal["a2y"]):
		raise ValueError("Calibration coefficients for the axes don't match")

	a0 =  cal["a0x"]/1e3 # in MeV
	a1 = cal["a1x"]/1e3 # in MeV

	strengthfile = open(folder+'/strength.nrm', 'r')
	strengthlines = strengthfile.readlines()
	strengthfile.close()
	N = len(strengthlines)
	strength = np.zeros((N,3))
	for i in range(N):
		words = strengthlines[i].split()
		if i < int(N/2):
			strength[i,0] = a0 + i*a1 # Energy coordinate
			strength[i,1] = float(words[0]) # Strength coordinate
		else:
			strength[i-int(N/2),2] = float(words[0]) # Strength uncertainty (this way due to horrible file format!)

	transextfile = open(folder+'/transext.nrm')
	transextlines = transextfile.readlines()
	transextfile.close()
	Next = len(transextlines)
	strengthext = np.zeros((Next, 2))
	for i in range(Next):
		transext_current = float(transextlines[i].split()[0])
		energy_current = a0 + i*a1 + 1e-8
		strengthext[i,:] = ( energy_current, transext_current/(2*np.pi*energy_current**3) )

	return strength, strengthext

def read_mama_2D(filename):
    # Reads a MAMA matrix file and returns the matrix as a numpy array,
    # as well as a list containing the calibration coefficients
    # and 1-D arrays of calibrated x and y values for plotting and similar.
    matrix = np.genfromtxt(filename, skip_header=10, skip_footer=1)
    cal = {}
    with open(filename, 'r') as datafile:
        calibration_line = datafile.readlines()[6].split(",")
        # a = [float(calibration_line[2][:-1]), float(calibration_line[3][:-1]), float(calibration_line[5][:-1]), float(calibration_line[6][:-1])]
        # JEM update 20180723: Changing to dict, including second-order term for generality:
        # print("calibration_line =", calibration_line, flush=True)
        cal = {"a0x":float(calibration_line[1]), "a1x":float(calibration_line[2]), "a2x":float(calibration_line[3]),
             "a0y":float(calibration_line[4]), "a1y":float(calibration_line[5]), "a2y":float(calibration_line[6])}
    # TODO: INSERT CORRECTION FROM CENTER-BIN TO LOWER EDGE CALIBRATION HERE.
    # MAKE SURE TO CHECK rebin_and_shift() WHICH MIGHT NOT LIKE NEGATIVE SHIFT COEFF.
    # (alternatively consider using center-bin throughout, but then need to correct when plotting.)
    Ny, Nx = matrix.shape
    y_array = np.linspace(0, Ny-1, Ny)
    y_array = cal["a0y"] + cal["a1y"]*y_array + cal["a2y"]*y_array**2
    x_array = np.linspace(0, Nx-1, Nx)
    x_array = cal["a0x"] + cal["a1x"]*x_array + cal["a2x"]*x_array**2
    # x_array = np.linspace(cal["a0x"], cal["a0x"]+cal["a1x"]*Nx, Nx) # BIG TODO: This is probably center-bin calibration,
    # x_array = np.linspace(a[2], a[2]+a[3]*(Ny), Ny) # and should be shifted down by half a bin?
                                                    # Update 20171024: Started changing everything to lower bin edge,
                                                    # but started to hesitate. For now I'm inclined to keep it as
                                                    # center-bin everywhere.
    return matrix, cal, y_array, x_array # Returning y (Ex) first as this is axis 0 in matrix language

# Plotting:
# fig, axes = plt.subplots(figsize=(10,9))
fig, axes = plt.subplots()

axes.tick_params("x", top="off")
axes.tick_params("y", right="off")
axes.yaxis.set_minor_locator(ticker.NullLocator())

strength, strengthext = readGSF(folder=".")
iBn = (np.abs(strengthext[:,0] - Bn)).argmin()
plt.errorbar(strength[:,0], strength[:,1], yerr=strength[:,2], label="Present work", fmt='d', color='red', zorder=0)
plt.plot(strengthext[0:iBn,0], strengthext[:iBn,1], ":", color='red'
         , label='and extrapolation'
         )

strength_export = strength[strength[:,1]>0]
# np.savetxt("strength_export.txt",strength_export, header="E y yerr")

# plt.xlim([0,15])
# plt.ylim([1e-8,5e-6])
plt.yscale('log')
plt.xlabel(r'$\gamma$-ray energy $E_\gamma$ [MeV]')
plt.ylabel(r'$\gamma$-ray strength [MeV$^{-3}$]')

handles, labels = axes.get_legend_handles_labels()
axes.legend(handles, labels, numpoints=1, fancybox=True, framealpha=0.5)

fig.tight_layout()

# plt.savefig('strength_pyplot.png')

plt.show()



