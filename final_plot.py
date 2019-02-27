import numpy as np
import matplotlib.pyplot as plt
import yt
import scipy.optimize as opt
import pdb
import scipy.interpolate
import density_calculator as dc
from matplotlib.patches import Rectangle

def get_data():
    
    """
    ds = yt.load('DD1150/DD1150')

    one_fifty_cut_sf, one_fifty_cut_g = dc.sigmaCalc(ds, 500, 500, 500, 150, (-2000, -2000, 0), (2001, 2001, 1))

    
    two_hundred_cut_sf, two_hundred_cut_g = dc.sigmaCalc(ds, 500, 500, 500, 200, (-2000, -2000, 0), (2001, 2001, 1))
    two_fifty_cut_sf, two_fifty_cut_g = dc.sigmaCalc(ds, 500, 500, 500, 250, (-2000, -2000, 0), (2001, 2001, 1))
    """
    GSD = np.loadtxt('GSD_step500_cutodd_150')
    
    SFR_150 = np.loadtxt('SFR_step500_cutoff_150')

    SFR_200 = np.loadtxt('SFR_step500_cutoff_200')

    SFR_250 = np.loadtxt('SFR_step500_cutoff_250')
    
    logGSD = no.log10(GSD)

    logSFR_150 = np.log10(SFR_150)

    logSFR_200 = np.log10(SFR_200)

    logSFR_250 = np.log10(SFR_250)


def final_plot():

    GSD = np.loadtxt('GSD_step500_cutoff_150')

    SFR_150 = np.loadtxt('SFR_step500_cutoff150')

    SFR_200 = np.loadtxt('SFR_step500_cutoff200')

    SFR_250 = np.loadtxt('SFR_step500_cutoff250')

    logGSD = np.log10(GSD)

    logSFR_150 = np.log10(SFR_150)

    logSFR_200 = np.log10(SFR_200)

    logSFR_250 = np.log10(SFR_250)

    SigmaSFR = np.loadtxt('for_jaimee/bolatto/fig6/logSsfrlogSgas_200pc.dat')
    logSigma_SFR = np.loadtxt('for_jaimee/bolatto/fig6/logSsfrlogSgas_200pc_x.dat')
    logSigmaSFR_SFR = np.loadtxt('for_jaimee/bolatto/fig6/logSsfrlogSgas_200pc_y.dat')
    extent_SFR = (logSigma_SFR[0]-0.5*(logSigma_SFR[1]-logSigma_SFR[0]),
                  logSigma_SFR[-1]+0.5*(logSigma_SFR[1]-logSigma_SFR[0]),
                  logSigmaSFR_SFR[0]-0.5*(logSigmaSFR_SFR[1]-logSigmaSFR_SFR[0]),
                  logSigmaSFR_SFR[-1]+0.5*(logSigmaSFR_SFR[1]-logSigmaSFR_SFR[0]))
    plt.imshow(SigmaSFR, extent=extent_SFR, origin='lower', interpolation='nearest', cmap='gray_r', aspect='auto')
    pgray = Rectangle((0,0), 1, 1, fc='#888888')


    plt.plot(logGSD, logSFR_150, 'o', label = 'cutoff = 150 MYR', alpha = 0.7, markersize = 4); plt.plot(logGSD, logSFR_200, 's', label = 'cutoff = 200 MYR', alpha = 0.7, markersize = 4); plt.plot(logGSD, logSFR_250, 'v', label = 'cutoff = 250 MYR', alpha = 0.7, markersize = 4)


    plt.xlabel(r'$\log\,\Sigma_{\mathrm{gas}}$ [$M_\odot$ pc$^{-2}$]')
    plt.ylabel(r'$\log\,\Sigma_{\mathrm{SFR}}$ [$M_\odot$ pc$^{-2}$ Myr$^{-1}$]')
    plt.gca().get_xaxis().set_ticks([.2,.4,.6,.8,1.0])
    plt.gca().get_xaxis().set_ticklabels(['0.2','0.4','0.6','0.8','1.0'])
    plt.xlim([0,1.2])
    plt.ylim([-7,0])
    plt.legend(loc = 'best')
    plt.show()
