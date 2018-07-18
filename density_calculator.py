import numpy as np
import matplotlib.pyplot as plt
import yt
import scipy.optimize as opt

def cell_patch(data, cent, size, height):

    x, y, z = cent


    #convert inputs from pc to code length
    x = float(((data.quan(x, 'pc') + data.quan(0.5, 'code_length')).in_units('code_length')).v)
    y = float(((data.quan(y, 'pc') + data.quan(0.5, 'code_length')).in_units('code_length')).v)
    z = float(((data.quan(z, 'pc') + data.quan(0.6, 'code_length')).in_units('code_length')).v)


    size = float((data.quan(size, 'pc').in_units('code_length')).v)
    height = float((data.quan(height, 'pc').in_units('code_length')).v)


    #create data object with input limits
    cell_box = data.box(left_edge = (x-size/2, y-size/2, z-height), right_edge = (x+size/2, y+size/2, z+height))

    """
    #Center coordinates
    x_cells = (ad['x'] - data.arr(0.5, 'code_length')).in_units('pc')

    y_cells = (ad['y'] - data.arr(0.5, 'code_length')).in_units('pc')

    z_cells = (ad['z'] - data.arr(0.6, 'code_length')).in_units('pc')
    #Find the cells in the patch

    xbounds_cells = np.logical_and(x_cells > x - size/2, x_cells < x + size/2)
    ybounds_cells = np.logical_and(y_cells > y - size/2, y_cells < y + size/2)
    zbounds_cells = np.logical_and(z_cells > 0 - height, z_cells < 0 + height)
    width_box = np.logical_and(xbounds_cells, ybounds_cells)
    cell_box = np.logical_and(width_box, zbounds_cells)

    print('total cells in x mask: ', np.sum(xbounds_cells))
    print('total cells in y mask: ', np.sum(ybounds_cells))
    print('total cells: ', np.sum(cell_box))
    """

    return cell_box

"""
def particle_patch(data, center, size, height):
    
    ad = data.all_data()
    x, y, z = center

    x = yt.YTQuantity(x, 'pc')
    y = yt.YTQuantity(y, 'pc')
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')


    x_particles = (ad['particle_position_x'] - data.arr(0.5, 'code_length')).in_units('pc')

    y_particles = (ad['particle_position_y'] - data.arr(0.5, 'code_length')).in_units('pc')

    z_particles = (ad['particle_position_z'] - data.arr(0.6, 'code_length')).in_units('pc')

    #Find all the particles in the patch    
    xbounds_particles = np.logical_and(x_particles > x - size/2, x_particles < x + size/2)
    ybounds_particles = np.logical_and(y_particles > y - size/2, y_particles < y + size/2)
    zbounds_particles = np.logical_and(z_particles > 0 - height, z_particles < 0 + height)
    particle_width = np.logical_and(xbounds_particles, ybounds_particles)
    particle_box = np.logical_and(particle_width, zbounds_particles)
    print('total particles in x mask: ', np.sum(xbounds_particles))
    print('total particles in y mask: ', np.sum(ybounds_particles))
    print('total particles: ', np.sum(particle_box))
    return particle_box

"""

def GasDensity(ds, center, size, height, temp):

    x,y,z = center

    #Create patch of cells
    box = cell_patch(ds, center, size, height)

    #Make inputs into YT Quantites
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')

    #Find cells in patch
    mass = box['cell_mass'].in_units('Msun')
#    patch_cells = cell_mass[box]

    #Create Temperature Cut
    cell_temp = box['temperature'].in_units('K')
#    cell_temp_patch = cell_temp[box]
    if temp == 'cold':
        cold_gas = cell_temp < 10**4

        #Add up mass of cells in patch under temp cut
        cold_gas_mass = mass[cold_gas]
        total_mass = np.sum(cold_gas_mass)

    elif temp == 'warm':
        warm_gas = cell_temp > 10**4
        warm_gas_mass = mass[warm_gas]
        total_mass = np.sum(warm_gas_mass)
    elif temp == 'all':
        total_mass = np.sum(mass)
    
    print(total_mass, size**2)
    #Calculate density (Divide by area)
    density = total_mass/(size**2)
    return density

def SFDensity(ds, center, size, height, max_age):

    #define center, establish YT Quantities
    x, y, z = center
    size = yt.YTQuantity(size,'pc')
    height = yt.YTQuantity(height, 'pc')
    max_age = yt.YTQuantity(max_age, 'Myr')
    
    #Create spatial patch
    box = cell_patch(ds, center, size, height)

    #Create array of masses & ages of all particles
    star_mass = box['particle_mass'].in_units('Msun')
    ages = box['creation_time']

    #Select new star particles in patch
    new_stars = ages > 0
    new_star_time = ages[new_stars]
    new_star_ages = (ds.current_time - new_star_time).in_units('Myr')


    young = new_star_ages < max_age
    stars = new_star_ages[young]

    #Get masses of new stars
    new_star_mass = star_mass[new_stars]
    young_star_mass = new_star_mass[young]

    #Find total mass of new stars, then surface density
    total_star_mass = np.sum(young_star_mass)
    sfrDen = (total_star_mass/((size**2)*max_age)).in_units('Msun/yr/kpc**2')

    return sfrDen


def sigmaCalc(ds, step, age, size, height, low, high):
    #Give in output file names
    sfFile = input('Enter SFR file name: ')
    gFile = input('Enter Gas file name: ')
    #Read in limits 
    low_x, low_y, low_z = low
    high_x, high_y, high_z = high
    #Create empty lists to hold gas SD and SFR SD values
    gd = []
    sfrd = []

    for x in np.arange(low_x, high_x, step): #find x of center for each patch
        for y in np.arange(low_y, high_y, step): #find y of center
            for z in np.arange(low_z, high_z, step):
                print('x: ', x)
                print('y: ', y)
                print('z: ', z)
                SF = SFDensity(ds, (x, y, z), size, height, age) #Find SFRSD
                G = GasDensity(ds, (x, y, z), size, height) #Find GSD
                print('SFRSD: ', SF)
                print('GSD: ', G)
                sfrd.append(SF) #Add to arrays
                gd.append(G) #Add to arrays
                
                sfsd = np.asarray(sfrd)
                gsd = np.asarray(gd)
                np.savetxt(sfFile, sfsd)
                np.savetxt(gFile, gsd)

    return sfsd, gsd



def expectedSFR_calculator(ds, center, size, height):

    #Create patch of cells
    box = cell_patch(ds, center, size, height)
    #Make inputs into YT Quantites
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')
    #Find cells in patch
    cell_sfr = box['expectedSFR'].in_units('Msun/yr/kpc**3')

    v_cells = box['cell_volume']

        #Create Temperature Cut
#    cell_temp = ad['temperature'].in_units('K')
#    cell_temp_patch = cell_temp[box]
#    cold_gas = cell_temp_patch < 10**4
    #Add up mass of cells in patch under temp cut
#    cold_gas_mass = patch_cells[cold_gas]
    cell_sfr = cell_sfr*v_cells
    total_sfr = np.sum(cell_sfr)
    density = (total_sfr/size**2).in_units('Msun/yr/kpc**2')
    return density

def expSFRSD(ds, step, size, height, low, high):
    #Give in output file names
    expsfFile = input('Enter expected SFR file name: ')
#    gFile = input('Enter GSD file name: ')
    #Read in limits 
    low_x, low_y, low_z = low
    high_x, high_y, high_z = high
    #Create empty lists to hold gas SD and SFR SD values

    sfrd = []

    for x in np.arange(low_x, high_x, step): #find x of center for each patch
        for y in np.arange(low_y, high_y, step): #find y of center
            for z in np.arange(low_z, high_z, step):
                print('x: ', x)
                print('y: ', y)
                SF = expectedSFR_calculator(ds, (x, y, z), size, height) #Find SFRSD

                print('SFRSD: ', SF)

                sfrd.append(SF) #Add to arrays

    
    sfsd = np.asarray(sfrd)

    np.savetxt(expsfFile, sfsd)

    return sfsd


def sigmaPlotter(g, sfr, esf, age_in_Myr, size_in_kpc, limit):
    
    age_in_yr = age_in_Myr*(10**6)
    upper_limit = np.zeros(len(sfr))
    for i in range(len(sfr)):
        if sfr[i] == 0:
            upper_limit[i] = 50./age_in_yr/(size_in_kpc**2)


    plt.loglog(g, esf, 'o', label = 'Analytic', color = 'b')            
    plt.loglog(g, sfr, 's', color = 'r', label = '{age_in_Myr} Myr'.format(age_in_Myr=age_in_Myr))

    if limit == True:
        plt.loglog(g, upper_limit, 'v', color = 'r')


    plt.ylabel('$\dot{\Sigma}_{SF}$ ($M_\odot  kpc^{-2} yr^{-1}$)')
    plt.xlabel('$\Sigma (M_\odot  pc^{-2})$')
#    plt.xlim((0, 50))
    plt.legend(loc = 'upper left')
    prompt = input('save? ')
    if prompt == ('yes') or prompt == ('y'):
        fname = input('Enter filename: ')
        plt.savefig(fname)
    else:
        pass


    plt.show()


#define constants                                                                                                                                                                                                                                                                          
G = yt.YTQuantity(6.67259e-8, 'cm**3/(g*s**2)')
cw = yt.YTQuantity(8.0e5, 'cm/s')
fw = 0.5
zd = 0.33
tCNM = yt.YTQuantity(243., 'K').in_cgs()
Eff = 0.01
SigmaSF_0 = yt.YTQuantity(2.5e-3, 'Msun/(Myr*pc**2)').in_cgs()
fc = 1.
kB = yt.YTQuantity(1.380658e-16, 'erg/K').in_cgs()
alpha = 5.

def Krumholz_model(ds, center, size, height):
    #Define Patch                                                                                                                                                                                                                                                                          
    box = cell_patch(ds, center, size, height)

    #Get values from simulation                                                                                                                                                                                                                                                            
    SigmaG = GasDensity(ds, center, size, height, temp = 'all')
    print('SigmaG: ', SigmaG)
    size = yt.YTQuantity(size, 'pc').in_cgs()
    height = yt.YTQuantity(height, 'pc').in_cgs()
    SigmaG_0 = (SigmaG/yt.YTQuantity(1., 'Msun/pc**2')).in_units('')
    Z = 0.1
    G0_guess = (yt.YTQuantity(10**-4, 'Msun/(yr*kpc**2)')/SigmaSF_0).in_units('')
    rhoSD = (np.sum(box['particle_mass'].in_units('g'))/((size**2)*(2*height))).in_cgs()


    def nonlinear_fH2(fH2, G0):

        Sigma_H2 = (fH2*SigmaG).in_units('Msun/pc**2').in_cgs()
        Sigma_HI = ((1-fH2)*SigmaG).in_units('Msun/pc**2').in_cgs()
        RH2 = Sigma_H2/Sigma_HI


        #Plug values into system of equations                                                                                                                                                                                                                                              
        nCNM_min = 31*G0*(1+3.1*(Z)**0.365)**-1

        nCNM_2p = yt.YTQuantity(23*G0*((1+3.1*Z**0.365)/4.1)**-1, 'cm**-3')

        Pth = ((np.pi*G*Sigma_HI**2)/(4.*alpha))*(1+2*RH2+np.sqrt(((1+2*RH2)**2)+(32*zd*alpha*fw*rhoSD*cw**2)/(np.pi*G*Sigma_HI**2)))

        nCNM_hydro = Pth/(1.1*kB*tCNM)

        nCNM = max([nCNM_hydro.in_units('cm**-3'), nCNM_2p.in_units('cm**-3')])

        n1 = nCNM/yt.YTQuantity(10., 'cm**-3')

        Chi = (7.2*G0/n1).in_units('')

        Tc = (0.066*fc*Z*SigmaG_0).in_units('')

        S = np.log(1 + 0.6*Chi + 0.01*Chi**2)/(0.6*Tc)

        if S < 2.:
            fH2_inferred = 1 - (3/4)*S/(1+0.25*S)
        else:
            fH2_inferred = 0

        return fH2 - fH2_inferred



    def nonlinear_G0(G0):
        fH2_found = opt.brentq(nonlinear_fH2, 0, 1, args = (G0))

        tff = yt.YTQuantity(31*SigmaG_0**(-1/4), 'Myr')

        SigmaSF = fH2_found*Eff*SigmaG/tff

        return SigmaSF/SigmaSF_0 - G0

    a = 100
    try:
        G0_found = opt.brentq(nonlinear_G0, 10**-5, a)
    except:
        G0_found = 0.




    

    return (G0_found * SigmaSF_0).in_units('Msun/yr/kpc**2')


def Krumholz_data(ds, step, size, height, low, high):
    #Give in output file names
    KFile = input('Enter file name: ')

    #Read in limits 
    low_x, low_y, low_z = low
    high_x, high_y, high_z = high
    #Create empty lists to hold gas SD and SFR SD values

    KSD = []

    for x in np.arange(low_x, high_x, step): #find x of center for each patch
        for y in np.arange(low_y, high_y, step): #find y of center
            for z in np.arange(low_z, high_z, step):
                print('x: ', x)
                print('y: ', y)
                SF = Krumholz_model(ds, (x, y, z), size, height) #Find SFRSD

                print('SFRSD: ', SF)

                KSD.append(SF) #Add to arrays

    
    Ksd = np.asarray(KSD)

    np.savetxt(KFile, Ksd)

    return Ksd
