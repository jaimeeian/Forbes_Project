import numpy as np
import matplotlib.pyplot as plt
import emcee
import math
from glob import glob
import pickle
import yt

def complete_example():
    ''' In this example we'll generate a bunch of fake datasets, then fit them and plot the results'''
    to_test = 6 # how many datasets should we generate?
    trueMs = np.random.uniform(size=to_test)*2 + 0.5 # true 
    truebs = np.random.uniform(size=to_test)*4 + -7 # true 
    truesigmasqs = np.power( np.random.uniform(size=to_test)*0.5+0.01, 2.0)
    ul = 2.0e-8 # the upper limit

    colors = ['k', 'r', 'b', 'g', 'orange', 'purple', 'gray', 'yellow', 'lightblue', 'maroon', 'lightgreen', 'pink'] * 20 # put together a bunch of random colors
    assert len(colors) >= to_test # make sure we have enough colors

    fig,ax = plt.subplots(nrows=2, ncols=2)
    for i in range(to_test):
        # for each set of parameters, generate a dataset...
        x = np.loadtxt('GSD_local/GSD_step500_cutoff250')
        y = np.loadtxt('SFR_local/SFR_step500_cutoff250')

        # ... fit the dataset ..
        # chain is a set of samples from the posterior distribution
        chain = fit_sim(x,y,upper_limit_value=ul)
        # from the chain, we can compute "68% credible intervals" i.e. the range in which we expect the answer to fall 68% of the time, i.e. a point and its asymmetric errorbars.
        mInterval, bInterval, sigmasqInterval = get_intervals(chain)
        # we can also plot a full set of diagnostics for each fit, which can be found in various files examples_*.pdf once the code is run
        # it's a good idea to plot these things, but it's not an absolute necessity.
        diagnostics(x,y,chain,upper_limit_value=ul, mTrue=trueMs[i], bTrue=truebs[i], sigmasqTrue=truesigmasqs[i], fn='step55_cutoff250'+str(i))

        # and plot the resulting best fit values and their uncertainties
        # here we make use of a function defined in this code, plotErrorbar, which plots a point and its asymmetric errorbars in both directions given the 68% credible intervals mentioned above.
        plotErrorbar( mInterval, bInterval, ax[0,0], c=colors[i] ) 
        plotErrorbar( mInterval, sigmasqInterval, ax[1,0], c=colors[i] )
        plotErrorbar( bInterval, sigmasqInterval, ax[1,1], c=colors[i] )

        # since we know the real values of m, b, and sigma in this case, plot those too
        ax[0,0].scatter(trueMs, truebs, c=colors) 
        ax[1,0].scatter(trueMs, truesigmasqs, c=colors) 
        ax[1,1].scatter(truebs, truesigmasqs, c=colors) 

        # finish up labelling the plot.
        ax[0,0].set_xlabel('m')
        ax[1,0].set_xlabel('m')
        ax[1,1].set_xlabel('b')

        ax[0,0].set_ylabel('b')
        ax[1,0].set_ylabel(r'$\sigma^2$')
        ax[1,1].set_ylabel(r'$\sigma^2$')
        fig.delaxes(ax[0,1])
        plt.tight_layout()
        plt.savefig('step500_cutoff250.pdf')
        plt.close(fig)


def test():
    ul = 4.0e-7 # upper limit
    m,b,sigma = 1.9, -7.0, 0.3 # true values of the model parameters
    x,y = generate_test_data(m=m, b=b, sigma=sigma, ul=ul) # generate some fake data

    chain = fit_sim(x,y,upper_limit_value=ul) # fit to the model

    # make some plots to see how well we recover the true answer.
    # The plots are saved as testfit_*pdf
    diagnostics(x,y,chain, upper_limit_value=ul, mTrue=m, bTrue=b, sigmasqTrue=sigma*sigma, fn='testfit')



def get_intervals(chain):
    # first truncate the chain *** Careful here, because the chain needs to have run long enough to be converged, and here we're assuming that has happened without checking.
    nwalkers, niter, ndim = np.shape(chain)
    chain = chain[:,int(niter/2)::20, -3:] # throw out the first half of the chain, and only keep every 20th sample from the rest. Also only keep the last three parameters, m,b, and sigma, and ignore the u's.
    samples = chain.reshape(-1,3) # This should now represent a good number of samples of the joint distribution of m, b, and sigma^2.

    # compute simple summary statistics and return them. 
    # Other things to add here might include the correlation matrix between quantities
    # These are the 68% credible intervals and medians of each of these quantities. Suitable for plotting errorbars
    mInterval = np.percentile( samples[:,0], [16, 50, 84] )
    bInterval = np.percentile( samples[:,1], [16, 50, 84] )
    sigmasqInterval = np.percentile( samples[:,2], [16, 50, 84] )

    return mInterval, bInterval, sigmasqInterval

def plotErrorbar(interval1, interval2, ax, **kwargs):
    ''' A convenience function to plot errorbars given intervals '''
    assert len(interval1)==3
    assert len(interval2)==3
    ax.errorbar( interval1[1], interval2[1], xerr=[[interval1[1]-interval1[0], interval1[2]-interval1[1]]], yerr=[[interval2[1]-interval2[0], interval2[1] - interval2[0]]], **kwargs )




def generate_test_data(N=64, m=1.9, b=-7.0, sigma=0.3, ul=4.0e-7, xmin=0.2, xmax=0.8):
    ''' generate some fake data to test how well the fit recovers the true parameters. 
        The defaults values will follow something very similar to the actual model we're using to fit the data.
        Returns an array of x and y values.
        The parameters are:
            N, the number of data points
            m, the slope of the linear fit in log-x vs log-y space
            b, the offset of the linear fit in log-x vs log-y space
            sigma, the 1 \sigma scatter in the relation in log-x vs log-y space, measured in dex,
            ul, the upper limit to impose on the data (in the same units as y, i.e. NOT logarithmic)
            xmin and xmax - the range of log-x-values to use to generate the data. These ARE logarithmic'''
    sigmasq = sigma*sigma
    x = np.power(10.0, xmin + (xmax-xmin)*np.random.random(size=N) )
    us = np.random.randn(N) * np.sqrt(sigmasq)
    lambdas = np.power(10.0, m*np.log10(x) + b  + us) / ul
    y = np.random.poisson( lambdas ) * ul

    return x,y

def diagnostics(x,y,chain,upper_limit_value=None, mTrue=None, bTrue=None, sigmasqTrue=None, fn='testfit'):
    # at this point x and y should be extremely analogous to our datasets of \Sigma and \dot{\Sigma}_{SFR} respectively.
    # Let's try out our fitting procedure!
    #chain = fit_sim( x, y, upper_limit_value = upper_limit_value ) ### chain is an array of size ~ # of walkers X # of iterations X # of dimensions. In our case this is 300 x 500 x 67
    nwalkers, niter, ndim = np.shape(chain)
    mMLE,bMLE,sigmasqMLE, residualsMLE = fit_mle(x,y,upper_limit_value = upper_limit_value)


    # Some diagnostic plots:

    # First up is a trace plot. This shows how well the MCMC has converged. It shows the ensemble's path for each variable
    # as a function of time. For the MCMC to be converged, the distribution of this ensemble must stop changing.
    fig,ax = plt.subplots(nrows=3)
    N = len(x) # the number of data points
    for k in range(nwalkers):
        ax[0].plot( chain[k, :, N], c='k', alpha=0.02 )
        ax[1].plot( chain[k, :, N+1], c='k', alpha=0.02 )
        ax[2].plot( chain[k, :, N+2], c='k', alpha=0.02 )

    if not mTrue is None:
        ax[0].plot( range(niter), [mTrue]*niter, c='r' ) # the true value
    ax[0].plot( range(niter), [mMLE]*niter, c='b' ) # the MLE value

    if not bTrue is None:
        ax[1].plot( range(niter), [bTrue]*niter, c='r' ) # the true value
    ax[1].plot( range(niter), [bMLE]*niter, c='b' ) # the MLE value

    if not sigmasqTrue is None:
        ax[2].plot( range(niter), [sigmasqTrue]*niter, c='r' ) # the true value
    ax[2].plot( range(niter), [sigmasqMLE]*niter, c='b' ) # the MLE value

    ax[2].set_xlabel('Iterations')
    ax[0].set_ylabel('m')
    ax[1].set_ylabel('b')
    ax[2].set_ylabel(r'$\sigma^2$')

    plt.savefig(fn+'_trace.pdf')
    plt.close(fig)


    # Next up, let's take a look at the data and the different fits.
    fig,ax = plt.subplots()
    xThis = np.logspace( np.log10(np.min(x)), np.log10(np.max(x)), 10 )
    for k in range(nwalkers):
        if k%2==0:
            mThis, bThis = chain[k,-1,N:N+2]
            yThis = np.power( 10.0, mThis*np.log10(xThis) + bThis )
            ax.plot( xThis, yThis, c='k', alpha=0.1 ) # draws from the posterior distribution
    if not mTrue is None and not bTrue is None:
        yThis = np.power( 10.0, mTrue*np.log10(xThis) + bTrue) # true line used to generate the data
    ax.plot( xThis, yThis, c='r', alpha=1 )
    yThis = np.power( 10.0, mMLE*np.log10(xThis) + bMLE) # maximum likelihood estimate
    ax.plot( xThis, yThis, c='b', alpha=1 )

    ax.scatter( x,y, marker='o', c='k' ) # plot the artificial data points
    zeros = y==0
    if not upper_limit_value is None:
        ax.scatter( x[zeros],[upper_limit_value]*np.sum(zeros), marker='v', c='k' ) # plot the artificial upper limits.
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$\Sigma$')
    ax.set_ylabel(r'$\dot{\Sigma}_\mathrm{SFR}$')
    plt.savefig(fn+'_data.pdf')
    plt.close(fig)


    burnIn = niter/2 # an initial guess - may require some adjustment based on inspecting the traceplot.
    import corner
    fig = corner.corner( chain[:, int(burnIn)::5, N:].reshape((-1,3)), labels=[r'm',r'b',r'$\sigma^2$'], truths = [mTrue,bTrue,sigmasqTrue])
    fig.savefig(fn+'_corner.pdf')
    plt.close(fig)




def fit_sim(x, y, upper_limit_value = None, niter=20000, nwalkers=300):
    ''' Fit x vs y data. The function expects x and y to be the same length and positive 
        i.e. don't take the logarithm of the data before passing it to this function.
        The function operates by assuming that the data points to be fit are of the form
            y_i = N_i * const
        where N_i is an integer of counts, e.g. particles, cells, etc. We are going to assume that
        the upper limit value is the same as the value of the constant, i.e. we are going to say
            N_i = y_i / upper_limit_value
        If no upper_limit_value is provided by the user, we'll guess that it's the lowest non-zero
        value in the dataset.'''
    assert len(x)==len(y)
    assert np.all(x>0)
    assert np.all(y>=0)
    if upper_limit_value is None:
        ul = np.min(y[y>0]) # guess that the upper limit is the lowest non-zero value in the dataset.
    else:
        ul = upper_limit_value

    Ns = np.around( y/ul )

    # ok, so now the idea is that our likelihood function, that is the probability distr. of 
    # y given the model parameters, or p(y|m,b,\vec{u}, sigma), will be 
    # prod_i^N 1/y_i! * (10^{mx+b+u_i})^y_i * exp( - 10^{mx+b+u_i} ) * exp(-(1/2) u_i^2/\sigma^2)
    # i.e. each value of N is taken to be drawn from a Poisson distribution whose expected value is 
    # 10^(mx+b+u_i), i.e. a straight line in log space with an extra scatter given by the u_i
    # u_i is taken to be normally distributed with an unknown variance \sigma^2. \sigma^2's prior is 
    # a standard non-informative inverse gamma distribution!

    m,b,sigmasq, residuals = fit_mle(x,y,upper_limit_value=upper_limit_value) # initial guess

    ndim, nwalkers = (len(x)+3, nwalkers)
    mle_est = np.zeros(ndim)
    mle_est[:len(x)] = residuals
    mle_est[len(x)] = m
    mle_est[len(x)+1] = b - np.log10(ul) # the MLE is done on the real data, but our likelihood uses Ns \equiv y/ul
    mle_est[len(x)+2] = sigmasq
    pos = [ mle_est * ( 1.0 + 1.0e-4*np.random.randn(ndim)) for i in range(nwalkers) ]

    sampler = emcee.EnsembleSampler( nwalkers, ndim, lnprob, args=(x,Ns) )
    sampler.run_mcmc( pos, niter)

    chain = sampler.chain
    ## recall that the shape of the chain should be such that: nwalkers, niter, ndim = np.shape(chain)
    chain[:,:,len(x)+1] += np.log10(ul) # factor the 'upper limit' value back in, i.e. convert back from N's, to y's.

    return chain


def lnprob(theta, x, y):
    ''' the natural log of the posterior probability density. 
            theta is the vector of parameters, and x and y are the data suitably transformed.
            In particular, x is still just \Sigma, the gas surface density,
            but y should be things that are very close to integers because in the calling 
            method -- fit_sim -- we have divided \dot{\Sigma}_SFR by the limiting value, i.e. we are
            roughly speaking counting particles or cells that have recently or will soon form stars.'''
    assert len(x)==len(y)
    assert len(theta)==len(x)+3 # random factors + m,b, and sigma.
    assert np.all(x>0)
    assert np.all(y>=0) # in this case y is treated as counts, so should all be ~integers with a minimum of zero.


    uvec = theta[:len(x)]
    m = theta[len(x)] # linear model in logx-logy:  logy = m * logx + b
    b = theta[len(x)+1]
    sigmasq = theta[len(x)+2]
    # this is the only hard boundary in the prior. Everything else is totalled up below.
    if sigmasq<=0:
        return - np.inf

    alpha,beta = (0.01, 0.01) # parameters for the prior on sigma^2
    lambdavec = np.power( 10.0, np.log10(x) * m + b + uvec )

    total = -0.5 * np.sum(uvec*uvec)/(sigmasq) 
    total += - (len(x)/2.0)*np.log(sigmasq)  
    total += np.sum(y*np.log(lambdavec)) 
    total += - np.sum( lambdavec )
    total += (alpha+1)*np.log(1.0/sigmasq) 
    total += - beta/sigmasq

    return total



def fit_mle(x, y, upper_limit_value = None):
    if upper_limit_value is None:
        ul = np.min(y[y>0])
    else:
        ul = upper_limit_value

    # in this case we're doing the simple thing, and just treating the upper limits as the value to be fit.
    logy = np.log10(y)
    logy[y<=ul] = np.log10(ul)

    logx = np.log10(x)

    coeff_array = np.ones( (len(x), 2 ) )
    coeff_array[:,0] = logx[:]
    params, _, _, _ = np.linalg.lstsq( coeff_array, logy[:] )
    m,b = params

    residuals = logy - m*logx - b
    sigmasq = np.std(residuals) ** 2

    return m,b,sigmasq, residuals

"""
if __name__=='__main__':
    #test() # run the test problem.
    complete_example()
"""


def chain_compiler(jump=1, start_point=0, upper_limit=None, niter=20000, gFile='GSD_massadj_noage.txt', sfDir='SFR_massadj_noage'):

    GSD= np.loadtxt(gFile)
    fnames = glob(sfDir+'/*')
    fnames.sort()
    SFR_list = np.asarray([np.loadtxt(f) for f in fnames])


    print('SFR_list: ', SFR_list)
    for i in np.arange(start_point, len(SFR_list)+jump, jump):
        SFR = SFR_list[i]
        print('SFR: ', SFR)
        chain = fit_sim(GSD, SFR, upper_limit_value = upper_limit, niter=niter)
        print(chain)
        ind = fnames[i][-7:-4]
        print('ind: ', ind)
        cshape = chain.shape
        mInterval, bInterval, sigmaInterval = get_intervals(chain)


        diagnostics(GSD, SFR, chain, upper_limit_value=upper_limit, fn='diagonostics_massadj_{ind}'.format(ind=ind))
        pickle.dump((mInterval, bInterval, sigmaInterval), open('intervals_massadj_{ind}.p'.format(ind=ind), 'wb'))


        del chain

        del mInterval, bInterval, sigmaInterval


def slope_v_time(directory='intervals', fn='m_v_time.pdf'):
    fnames = glob(directory+'/*')
    fnames.sort()
    fig,ax = plt.subplots()
    time = np.array([np.arange(10, 350, 10)]).T
    timeInterval = np.zeros((len(time), 3))
    for i in np.arange(0, len(time)):
        timeInterval[i] = time[i]


    for k,f in enumerate(fnames):
        ind = f[-5:-2]
        intervals = pickle.load(open(f, 'rb'))
        mInterval, bInterval, sigmaInterval = intervals


        print('mInterval: ', mInterval)

        t = np.zeros(3)+int(ind)

        #plotErrorbar(timeInterval[k], mInterval, ax)
        plotErrorbar(t, mInterval, ax)


        del mInterval, bInterval, sigmaInterval

        
    plt.xlabel('time')
    plt.ylabel('slope')
    #plt.ylim((0, 6))
    plt.savefig(fn)


def slope_v_radius(directory='intervals_r', fn='m_v_r.pdf'):

    fnames = glob(directory+'/*')
    fnames.sort()
    fig,ax = plt.subplots()



    for k,f in enumerate(fnames):
        ind = f[-6:-2]
        intervals = pickle.load(open(f, 'rb'))
        mInterval, bInterval, sigmaInterval = intervals


        print('mInterval: ', mInterval)

        t = np.zeros(3)+int(ind)

        #plotErrorbar(timeInterval[k], mInterval, ax)                                                                                        
        plotErrorbar(t, mInterval, ax)


        del mInterval, bInterval, sigmaInterval


    plt.xlabel('r')
    plt.ylabel('slope')
    #plt.ylim((0, 6))                                                                                                                        
    plt.savefig(fn)

    
    
    
def radial_values(jump=1, start_point=0, upper_limit=None, niter=20000, gFile='GSD_massadj_noage.txt', sfFile='SFR_massadj_noage/SFR_massadj_noage.txt'):

    GSD= np.loadtxt(gFile)
    SFR = np.loadtxt(sfFile)
    coords = np.loadtxt('coordinates.txt')
    radii = []
    for i in coords:
        x, y, z = i
        r = np.sqrt((x**2)+(y**2)+(z**2))
        radii.append(r)

    rad = []
    for r in radii:
        if r not in rad:
            rad.append(r)
            
    for r in rad:
        inds = np.where(radii<=r)[0]
        print(r, inds)
        ind = int(r)
        gsd_r = GSD[inds]
        sfr_r = SFR[inds]
        print('gsd: ', gsd_r)
        print('sfr: ', sfr_r)
        print(gsd_r.shape, sfr_r.shape, gsd_r.shape==sfr_r.shape)
        np.savetxt('GSD_radii_inclusive/GSD_r_inclusive_{0:04d}.txt'.format(ind), gsd_r)
        np.savetxt('SFR_radii_inclusive/SFR_r_inclusive_{0:04d}.txt'.format(ind), sfr_r)
        


def chain_compiler_r(jump=1, start_point=0, upper_limit=None, niter=20000, gDir='GSD_radii', sfDir='SFR_radii'):
    gNames = glob(gDir+'/*')
    gNames.sort()
    sfNames = glob(sfDir+'/*')
    sfNames.sort()
    print(gNames)
    print(sfNames)

    for i in np.arange(start_point, len(gNames)+jump, jump):
        GSD = np.loadtxt(gNames[i])
        SFR = np.loadtxt(sfNames[i])
        if GSD.size==1:
            GSD = np.asarray([GSD.tolist()])
            SFR = np.asarray([SFR.tolist()])
        else:
            pass


        print('GSD: ', GSD)
        print('SFR: ', SFR)
        print(GSD.size)
        print(type(GSD))

        chain = fit_sim(GSD, SFR, upper_limit_value = upper_limit, niter=niter)
        print(chain)
        ind = gNames[i][-8:-4]
        print('ind: ', ind)
        cshape = chain.shape
        mInterval, bInterval, sigmaInterval = get_intervals(chain)


        diagnostics(GSD, SFR, chain, upper_limit_value=upper_limit, fn='diagonostics_radial_inclusive_{ind}'.format(ind=ind))
        pickle.dump((mInterval, bInterval, sigmaInterval), open('intervals_radial_inclusive_{ind}.p'.format(ind=ind), 'wb'))


        del chain

        del mInterval, bInterval, sigmaInterval


        
