import numpy as np
import matplotlib.pyplot as plt
import yt
import scipy.optimize as opt
import pdb
import scipy.interpolate
#import density_calculator as dc
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
    
    
    logGSD = np.log10(GSD)

    logSFR_150 = np.log10(SFR_150)

    logSFR_200 = np.log10(SFR_200)

    logSFR_250 = np.log10(SFR_250)


def final_plot():
    #Load in simulation data
    GSD = np.loadtxt('GSD_step500_cutoff_150')

    SFR_150 = np.loadtxt('SFR_step500_cutoff150')

    SFR_200 = np.loadtxt('SFR_step500_cutoff200')

    SFR_250 = np.loadtxt('SFR_step500_cutoff250')
    
    #Gas Based
    eSF = np.loadtxt('eSFR_step500_range2000.txt')

    #Load in analytic data
    Kdata = np.loadtxt('Krumholz_step500_range2000')
    Odata = np.loadtxt('Ostriker_step500_range2000')
    FGdata = np.loadtxt('FG_step500_range2000.txt')

    #Convert to logspace
    logGSD = np.log10(GSD)

    logSFR_150 = np.log10(SFR_150)

    logSFR_200 = np.log10(SFR_200)

    logSFR_250 = np.log10(SFR_250)

    logeSF = np.log10(eSF)

    logK = np.log10(Kdata)
    logO = np.log10(Odata)
    logFG = np.log10(FGdata)

    plt.figure(figsize = (8,8))
    
    #Plot Bolatto data
    SigmaSFR = np.loadtxt('for_jaimee/bolatto/fig6/logSsfrlogSgas_200pc.dat')
    logSigma_SFR = np.loadtxt('for_jaimee/bolatto/fig6/logSsfrlogSgas_200pc_x.dat')
    logSigmaSFR_SFR = np.loadtxt('for_jaimee/bolatto/fig6/logSsfrlogSgas_200pc_y.dat')
    extent_SFR = (logSigma_SFR[0]-0.5*(logSigma_SFR[1]-logSigma_SFR[0]),
                  logSigma_SFR[-1]+0.5*(logSigma_SFR[1]-logSigma_SFR[0]),
                  logSigmaSFR_SFR[0]-0.5*(logSigmaSFR_SFR[1]-logSigmaSFR_SFR[0]),
                  logSigmaSFR_SFR[-1]+0.5*(logSigmaSFR_SFR[1]-logSigmaSFR_SFR[0]))
    plt.imshow(SigmaSFR, extent=extent_SFR, origin='lower', interpolation='nearest', cmap='gray_r', aspect='auto')
    pgray = Rectangle((0,0), 1, 1, fc='#888888')

    #Overplot depletion times of 1, 10, 100 Gyr
    logSigma_g = np.linspace(0,1.2,100)
    plt.plot(logSigma_g, logSigma_g-3, 'k--', lw=2, alpha=0.5)
    plt.plot(logSigma_g, logSigma_g-4, 'k--', lw=2, alpha=0.5)
    plt.plot(logSigma_g, logSigma_g-5, 'k--', lw=2, alpha=0.5)
    idx = 70
    plt.text(logSigma_g[idx], logSigma_g[idx]-2.6,
             r'$t_{\mathrm{dep}} = 1$ Gyr',
             rotation=10)
    plt.text(logSigma_g[idx], logSigma_g[idx]-3.6,
             r'$t_{\mathrm{dep}} = 10$ Gyr',
             rotation=10)
    plt.text(logSigma_g[idx], logSigma_g[idx]-4.6,
             r'$t_{\mathrm{dep}} = 100$ Gyr',
             rotation=10)
    
    #actual plotting
    plt.plot(logGSD, logSFR_150, 'o', label = 'cutoff = 250 MYR', alpha = 0.7); plt.plot(logGSD, logSFR_200, 'o', label = 'cutoff = 200 MYR', alpha = 0.7); plt.plot(logGSD, logSFR_250, 'o', label = 'cutoff = 150 MYR', alpha = 0.7)
    
    plt.plot(logGSD, logeSF, 'o', label = 'Gas Based', alpha = 0.7)

    plt.plot(logGSD, logK, 's', label = 'Krumholz', alpha = 0.7)
    plt.plot(logGSD, logFG, 's', label = 'Faucher-Giguere', alpha = 0.7)
    plt.plot(logGSD, logO, 's', label = 'Ostriker', alpha = 0.7)
    
    #Labeling/Tweaking
    plt.style.use('classic')
    plt.xlabel(r'$\log\,\Sigma_{\mathrm{gas}}$ [$M_\odot$ pc$^{-2}$]')
    plt.ylabel(r'$\log\,\Sigma_{\mathrm{SFR}}$ [$M_\odot$ pc$^{-2}$ Myr$^{-1}$]')
    plt.gca().get_xaxis().set_ticks([.2,.4,.6,.8,1.0])
    plt.gca().get_xaxis().set_ticklabels(['0.2','0.4','0.6','0.8','1.0'])
    plt.xlim([0,0.9])
    plt.ylim([-7,0])
    plt.legend(loc = 'upper left')
    
    
    #ax = plt.gca()
    #ax.set_aspect('=')

    #Save/Show
    plt.show()
    #plt.savefig('master_plot_v2.png')
