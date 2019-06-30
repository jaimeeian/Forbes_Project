import numpy as np
import matplotlib.pyplot as plt
import yt
from yt.data_objects.particle_filters import add_particle_filter
import scipy.optimize as opt
import pdb
import scipy.interpolate
from sfr_calculator import formed_star

add_particle_filter("formed_star", function=formed_star, filtered_type='all',
                    requires=["creation_time"])

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

def SFDensity(ds, center, size, height, cutoff):

    #define center, establish YT Quantities
    x, y, z = center
    size = yt.YTQuantity(size,'pc')
    height = yt.YTQuantity(height, 'pc')
    #max_age = yt.YTQuantity(max_age, 'Myr')
    cutoff = yt.YTQuantity(cutoff, 'Myr')
    max_age = (ds.current_time.in_units('Myr') - cutoff).in_units('Myr')

    #Create spatial patch
    box = cell_patch(ds, center, size, height)

    #Create array of masses & ages of all particles
    star_mass = box['particle_mass'].in_units('Msun')
    creation_time = box['creation_time']

    #Select new star particles in patch
    new_stars = creation_time > 0
    new_star_time = creation_time[new_stars].in_units('Myr')
    new_star_ages = (ds.current_time.in_units('Myr') - new_star_time).in_units('Myr')
    birth_years = creation_time[new_stars].in_units('Myr')

    #Create cutoff after SFR is stabilized (~200MYR)
#    cutoff = yt.YTQuantity(200, 'Myr')
#    cut_mask = new_star_ages < cutoff
#    cut_ages = new_star_ages[cut_mask]

    #young = new_star_ages < max_age
    #stars = new_star_ages[young]
    #stars = new_star_ages[young]
    young = birth_years > cutoff
    stars = birth_years[young]
    


    #Get masses of new stars
    new_star_mass = star_mass[new_stars]
#    new_star_mass = new_star_mass[cut_mask]
    young_star_mass = new_star_mass[young]

    #Find total mass of new stars, then surface density
    total_star_mass = np.sum(young_star_mass)
    sfrDen = (total_star_mass/((size**2)*max_age)).in_units('Msun/yr/kpc**2')

    return sfrDen


def sigmaCalc(ds, step, size, height, age, low, high):
    #print stats
    print('step size: ', step)
    print('size: ', size)
    print('cutoff: ', age)
    print('range: ', low, 'to', high)

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
                G = GasDensity(ds, (x, y, z), size, height, temp='all') #Find GSD
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
t_SF = yt.YTQuantity(2e9, 'yr')
Pth_0 = kB*yt.YTQuantity(3000, 'K*cm**-3')

def Ostriker_model(ds, center, size, height):
    box = cell_patch(ds, center, size, height)
  
    SigmaG = GasDensity(ds, center, size, height, temp = 'all')
    print('SigmaG: ', SigmaG)
    size = yt.YTQuantity(size, 'pc').in_cgs()
    height = yt.YTQuantity(height, 'pc').in_cgs()
    SigmaG_0 = (SigmaG/yt.YTQuantity(1., 'Msun/pc**2')).in_units('')
    Z = 0.1
    rhoSD = (np.sum(box['particle_mass'].in_units('g'))/((size**2)*(2*height))).in_cgs()
    print('rhoSD: ', rhoSD)

    Sigma_h = (2*alpha*Pth_0)/(np.pi*G*SigmaSF_0*t_SF)

    S = (8.*zd*alpha*fw*rhoSD*cw**2)/(np.pi*G*SigmaG**2)

    omega = SigmaG/Sigma_h
    print('Sigma_h: ', Sigma_h)
    print('s: ', S)
    print('omega: ', omega)

    def Gdiff_frac(x):
        return omega*(1.+(1.+(1./(omega*(1.-x)))+(S/(1.-x)**2.))**(1./2.))-(1/x)

    x = opt.brentq(Gdiff_frac, 10**-5, 1)
    print('x: ', x)
    Sigma_diff = x*SigmaG
    Sigma_GBC = (1-x)*SigmaG
    print('Sigma_diff', Sigma_diff)
    print('Sigma_GBC', Sigma_GBC)

        #diffuse-gas thermal equilibrium
    
#    phi_d = (1./4.)*(1+3*(Z*SigmaG/SigmaG_0)**0.4)
    
 #   Pth = (1-phi_d)*(Pth_0/SigmaSF_0)*(Sigma_GBC/t_SF)

  #  Pth = ((np.pi*G*Sigma_diff**2)/(4*alpha))*(1+2*(Sigma_GBC/Sigma_dff)+((1+2*(Sigma_GBC/Sigma_diff))**2+((32*zd*alpha*fw*cw**2)/(np.pi*G))*(rhoSD/(Sigma_diff**2)))**(1/2))

    SigmaSF = (Sigma_GBC/t_SF)
    
    print('tSF: ', t_SF)
    print('SigmaSFR: ', SigmaSF)

    return SigmaSF.in_units('Msun/yr/kpc**2')
"""

def Ostriker_model(ds, center, size, height):
    box = cell_patch(ds, center, size, height)
  


    SigmaG = GasDensity(ds, center, size, height, temp = 'all')
    print('SigmaG: ', SigmaG)
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')
    SigmaG_0 = (SigmaG/yt.YTQuantity(1., 'Msun/pc**2')).in_units('')
    Z = 0.1

    rhoSD = (np.sum(box['particle_mass'].in_units('g'))/((size**2)*(2*height))).in_cgs()

    SigmaSF_low = yt.YTQuantity(3e-10, 'Msun/pc**2/yr')*((SigmaG/yt.YTQuantity(10., 'Msun*pc**-2'))*(1.+3.*((Z*SigmaG)/(yt.YTQuantity(10., 'Msun*pc**-2')))**0.4)*((2./alpha)*(SigmaG/yt.YTQuantity(10., 'Msun*pc**-2'))+np.sqrt(50*fw/alpha)*np.sqrt(rhoSD/yt.YTQuantity(0.1, 'Msun*pc**-3'))))

    SigmaSF = ((t_SF/SigmaG)+(1./SigmaSF_low))**-1.

    return SigmaSF
"""

def Ostriker_data(ds, step, size, height, low, high):
    #Give in output file names                                                                              
    OFile = input('Enter file name: ')

    #Read in limits                                                                                         
    low_x, low_y, low_z = low
    high_x, high_y, high_z = high
    #Create empty lists to hold gas SD and SFR SD values                                                    

    OSD = []


    for x in np.arange(low_x, high_x, step): #find x of center for each patch                               
        for y in np.arange(low_y, high_y, step): #find y of center                                          
            for z in np.arange(low_z, high_z, step):
                print('x: ', x)
                print('y: ', y)
                SF = Ostriker_model(ds, (x, y, z), size, height) #Find SFRSD                                

                print('SFRSD: ', SF)

                OSD.append(SF) #Add to arrays                                                               


    Osd = np.asarray(OSD)

    np.savetxt(OFile, Osd)

    return Osd



def _rgas(field, data):
    ''' Velocity of gas towards or away from the center of the galaxy, in cylindrical units, i.e. v dot r, where r is the cylindrical radius sqrt(x^2+y^2). '''
    x = data['x'] - data.ds.arr(0.5, 'code_length')
    y = data['y'] - data.ds.arr(0.5, 'code_length')
    r = np.sqrt(x*x + y*y)
    return r

def _rParticle(field, data):
    ''' Velocity of gas towards or away from the center of the galaxy, in cylindrical units, i.e. v dot r, where r is the cylindrical radius sqrt(x^2+y^2). '''
    x = data['particle_position_x'] - data.ds.arr(0.5, 'code_length')
    y = data['particle_position_y'] - data.ds.arr(0.5, 'code_length')
    r = np.sqrt(x*x + y*y)
    return r

def _vrParticle(field, data):
    ''' Velocity of gas towards or away from the center of the galaxy, in cylindrical units, i.e. v dot r, where r is the cylindrical radius sqrt(x^2+y^2). '''
    x = data['particle_position_x'] - data.ds.arr(0.5, 'code_length')
    y = data['particle_position_y'] - data.ds.arr(0.5, 'code_length')
    r = np.sqrt(x*x + y*y)
    vx = data['particle_velocity_x']
    vy = data['particle_velocity_y']
    return vx*x/r + vy*y/r 

def _vthetaParticle(field,data):
    ''' Velocity of gas around the center of the galaxy, in cylindrical units, i.e. v dot theta hat, where theta is the angular position of the cell in the x-y plane with respect to the x-axis. '''
    x = data['particle_position_x'] - data.ds.arr(0.5, 'code_length')
    y = data['particle_position_y'] - data.ds.arr(0.5, 'code_length')
    r = np.sqrt(x*x + y*y)
    vx = data['particle_velocity_x']
    vy = data['particle_velocity_y']
    return vy*x/r - vx*y/r


def _vr(field, data):
    ''' Velocity of gas towards or away from the center of the galaxy, in cylindrical units, i.e. v dot r, where r is the cylindrical radius sqrt(x^2+y^2). '''
    x = data['x'] - data.ds.arr(0.5, 'code_length')
    y = data['y'] - data.ds.arr(0.5, 'code_length')
    r = np.sqrt(x*x + y*y)
    vx = data['x-velocity']
    vy = data['y-velocity']
    return vx*x/r + vy*y/r 

def _vtheta(field,data):
    ''' Velocity of gas around the center of the galaxy, in cylindrical units, i.e. v dot theta hat, where theta is the angular position of the cell in the x-y plane with respect to the x-axis. '''
    x = data['x'] - data.ds.arr(0.5, 'code_length')
    y = data['y'] - data.ds.arr(0.5, 'code_length')
    r = np.sqrt(x*x + y*y)
    vx = data['x-velocity']
    vy = data['y-velocity']
    return vy*x/r - vx*y/r

def _expectedSFR(field, data):
    # The density must exceed the Jeans density on the finest resolution level, i.e. stars can only form in cells that 
    # obey this condition:
    cond = data['density'] > data.ds.gamma*np.pi*yt.units.kboltz*data['temperature']/(16.0 * yt.units.G * 1.1 * yt.units.mass_hydrogen* data['dx']*data['dx'])
    # all other cells are have zero SFR density. Because of the way YTArrays work, we need to define an array of
    # zeros with the correct dimensions, i.e. mass/volume/time:
    ret = yt.YTArray(np.zeros(np.shape(data['density'])), 'g/cm**3/s')
    # cells that meet the condition have a SFR density given by 0.01 * density/freefall time.
    ret[cond] = 0.01 * data['density'][cond] / np.sqrt(3.0*np.pi/(32.0 * yt.units.G* data['density'][cond])).in_units('s')
    return ret

yt.add_field( "r-velocity", function=_vr, units='km/s')
yt.add_field( "r", function=_rgas, units='kpc')
yt.add_field( "particle_r", function=_rParticle, units='kpc')
yt.add_field( "theta-velocity", function=_vtheta, units='km/s')
yt.add_field( "particle_velocity_r", function=_vrParticle, units='km/s')
yt.add_field( "particle_velocity_theta", function=_vthetaParticle, units='km/s')
yt.add_field( "expectedSFR", function=_expectedSFR, units='msun/yr/pc**3')


class massProfile:
    def __init__(self, center, ds):
        self.cached=False
        self.center = center
        self.f = None
        self.ds = ds
    def massInterior(self, r):
        if not self.cached:
            self.cache()
        # the interpolation function os for logr-logM
        print('r: ', r.in_units('kpc'))
        return yt.YTQuantity(np.power(10.0, self.f(np.log10(r.in_units('kpc')))), 'msun')
    def vcirc(self,r):
        return np.sqrt(yt.units.G*self.massInterior(r)/r)
    def beta(self, r):
        # logarithmic rotation curve slope dlnv/dlnr
        if not self.cached:
            self.cache()
        return (np.log(self.vcirc(r*1.01))-np.log(self.vcirc(r*0.99)))/(np.log(r*1.01)-np.log(r*0.99))


    def cache(self):
        print(self.center)
        self.cached=True
        ad = self.ds.all_data()
        rsquaredParticles = np.power(ad['particle_position_x']-self.center[0],2) + np.power(ad['particle_position_y']-self.center[1],2) + np.power(ad['particle_position_z']-self.center[2],2)
        rsquaredCells= np.power(ad['x']-self.center[0],2) + np.power(ad['y']-self.center[1],2) + np.power(ad['z']-self.center[2],2)
        npts = 100
        rs = yt.YTArray( np.power(10.0, np.linspace(-2, 0, npts)), 'kpc')
        ms = yt.YTArray( np.zeros(npts), 'msun')
        for i in range(npts):
            theseParticles = rsquaredParticles<rs[i]*rs[i]
            theseCells = rsquaredCells<rs[i]*rs[i]
            massInParticles = np.sum( ad['particle_mass'][theseParticles] )
            massInCells = np.sum( ad['cell_mass'][theseCells] )
            ms[i] = (massInParticles+massInCells).in_units('msun')
            print('rs: ', rs)
            print('ms: ', ms)
            print('log10(rs): ', np.log10(rs))
            print('log10(ms): ', np.log10(ms))
            print('rs limits: ', max(np.log10(rs)), min(np.log10(rs)))
            print('ms limits: ', max(np.log10(ms)), min(np.log10(ms)))

        self.f = scipy.interpolate.interp1d( np.log10(rs), np.log10(ms), kind='quadratic', fill_value="extrapolate")



def Qpatch(data, width, theMassProfile ):
    ''' Estimate multi-component Toomre Q within a patch of a disk. The inputs are:
            data: a YT data object, or at least a dictionary containing all the necessary data fields. No further filtering should be required of this data object, so this should be a YTRegion or similar object, not all_data().
            width: a YTQuantity containing the width of the (assumed-to-be-square) patch of the disk
            theMassProfile: an instance of massProfile
            '''
    # the mass-weighted average cylindrical radius of this patch
    r = np.average( data['r'], weights=data['cell_mass'] ).in_units('kpc')

    # the angular velocity at that radius

    Omega = (np.sqrt(yt.units.G*theMassProfile.massInterior(r.in_units('kpc'))/(r*r*r))).in_units('s**-1')

    # the epicyclic frequency
    kappa = (np.sqrt(2.0*(theMassProfile.beta(r.in_units('kpc'))+1.0)) * Omega).in_units('s**-1')

    # mass-weighted average velocity in the r direction within this patch
    meanVelocity = (np.average(data['r-velocity'], weights=data['cell_mass'])).in_units('km/s')
    # mass-weighted velocity dispersion in the r direction including the contribution of the sound speed
    sigma = (np.sqrt( np.sum(( np.power(data['sound_speed'],2.0) + np.power(data['r-velocity']-meanVelocity,2))*data['cell_mass'] )  / np.sum(data['cell_mass']) )).in_units('km/s')
    # the column density
    col = (np.sum( data['cell_mass'])/(width*width)).in_units('Msun/pc**2')
    # the gas-component of Q.
    Qg = (kappa*sigma/(np.pi*yt.units.G*col)).in_units('')


    # pick out stars [as opposed to dark matter]
    stars = data['particle_type']==2
    # average radial velocity of the stars
    avgr = (np.average(data['particle_velocity_r'][stars], weights=data['particle_mass'][stars])).in_units('km/s')
    # velocity dispersion of the stars in the radial direction
    sigmast = (np.sqrt( np.sum( np.power(data['particle_velocity_r'][stars]-avgr,2.0)*data['particle_mass'][stars] ) / np.sum(data['particle_mass'][stars]) )).in_units('km/s')
    # average z-velocity of the stars
    avgz = (np.average(data['particle_velocity_z'][stars], weights=data['particle_mass'][stars])).in_units('km/s')
    # velocity dispersion of the stars in the z-direction
    sigmazz = (np.sqrt( np.sum( np.power(data['particle_velocity_z'][stars]-avgz,2.0)*data['particle_mass'][stars] ) / np.sum(data['particle_mass'][stars]) )).in_units('km/s')
    # column density of stars
    colst = (np.sum( data['particle_mass'][stars] ) / (width*width)).in_units('msun/pc**2')
    # stellar part of Q
    Qst = (kappa*sigmast/(np.pi*yt.units.G*colst)).in_units('')

    # Formula for Q estimated by Romeo & Wiegert (2011)
    W = (2.0*sigmast*sigma/(sigmast**2+sigma**2)).in_units('') # weighting
    Tg = 1.5 # thickness correction for gas
    Tst = (0.8 + sigmazz/sigmast).in_units('') # thickness correction for stars
    if Qst*Tst > Qg*Tg:
        QRW = 1.0/ ( W/Qst + 1.0/Qg)
    else:
        QRW = 1.0/ ( 1.0/Qst + W/Qg)
    return QRW

def FG_model(ds, center, size, height, mP):
    try:
        box = cell_patch(ds, center, size, height)
#    gal_center = ds.arr([0.5, 0.5, 0.6], 'code_length')
#    mP = massProfile(gal_center, ds)
        size = yt.YTQuantity(size, 'pc')
        Q = Qpatch(box, size, mP)

        phi = 1./Q
        F = 1.
        retMomentum = yt.YTQuantity(3000., 'km/s')
        SigmaG = GasDensity(ds,center, size, height, temp = 'all')

        SigmaSF = (2.*np.sqrt(2.)*np.pi*G*Q*phi/F)*((retMomentum)**-1)*(SigmaG**2)


 #    SigmaSF = yt.YTQuantity(13., 'Msun/yr/kpc**2')*(Q*phi/F)*((retMomentum/yt.YTQuantity(3000., 'km/s'))**-1)*(SigmaG/yt.YTQuantity(10**3, 'Msun/pc**2'))**2


        print('SigmaG: ', SigmaG)
        print('SigmaSFR: ', SigmaSF)
        return SigmaSF.in_units('Msun/yr/kpc**2')

    except:
        pass

def FG_data(ds, step, size, height, low, high, mP):
    #Give in output file names                                                                              
    FGFile = input('Enter file name: ')

    #Read in limits                                                                                         
    low_x, low_y, low_z = low
    high_x, high_y, high_z = high
    #Create empty lists to hold gas SD and SFR SD values                                                    

    FGSD = []


    for x in np.arange(low_x, high_x, step): #find x of center for each patch                               
        for y in np.arange(low_y, high_y, step): #find y of center                                          
            for z in np.arange(low_z, high_z, step):
                print('x: ', x)
                print('y: ', y)
                SF = FG_model(ds, (x, y, z), size, height, mP) #Find SFRSD                                

                print('SFRSD: ', SF)

                FGSD.append(SF) #Add to arrays                                                               


    FGsd = np.asarray(FGSD)

    try:
        np.savetxt(FGFile, FGsd)
    except:
        np.savetxt(FGFile, FGsd, fmt='%s')


    return FGsd



def SFDensity_age_cut(ds, center, size, height, age_range):

    #define center, establish YT Quantities
    x, y, z = center
    size = yt.YTQuantity(size,'pc')
    height = yt.YTQuantity(height, 'pc')
    min_age, max_age = age_range
    min_age = yt.YTQuantity(min_age, 'Myr')
    max_age = yt.YTQuantity(max_age, 'Myr')

    #Create spatial patch
    box = cell_patch(ds, center, size, height)

    #Select new star particles in patch
    formed_star_birth = box['formed_star', 'creation_time'].in_units('Myr')
    formed_star_ages = (ds.current_time.in_units('Myr') - formed_star_birth).in_units('Myr')
    formed_star_mass = box['formed_star', 'particle_mass']

    #Select stars between min and max age
    low_filter = formed_star_ages > min_age
    high_filter = formed_star_ages < max_age
    filtered_mass = formed_star_mass[low_filter & high_filter]
    
    #Find total mass of stars, then SD
    total_star_mass = np.sum(filtered_mass)
    sfrDen = (total_star_mass/((size**2)*max_age)).in_units('Msun/yr/kpc**2')

    return sfrDen


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


def sigmaCalc_age(ds, step, size, height, age_range, low, high, age_bins = 100):
    #print stats
    print('step size: ', step)
    print('size: ', size)
    print('age range: ', age_range[0], 'to', age_range[1])
    print('range: ', low, 'to', high)

    #Give in output file names
    sfFile = input('Enter SFR file name: ')
    gFile = input('Enter Gas file name: ')
    #Read in limits 
    low_x, low_y, low_z = low
    high_x, high_y, high_z = high
    #Create empty lists to hold gas SD and SFR SD values
    gd = []
    sfrd = []

    for age in np.arange(age_range[0], age_range[1], age_bins):
        print('age',age)
        for x in np.arange(low_x, high_x, step): #find x of center for each patch
            print(x)
            for y in np.arange(low_y, high_y, step): #find y of center
                print(y)
                for z in np.arange(low_z, high_z, step):
                    print(z)
                    if age == age_range[0]:
                        print('x: ', x)
                        print('y: ', y)
                        print('z: ', z)
                
                        G = GasDensity(ds, (x, y, z), size, height, temp='all') #Find GSD
                        print('GSD: ', G)
                        gd.append(G) #Add to arrays
                        gsd = np.asarray(gd)
                        np.savetxt(gFile+'.txt', gsd)
                    else:
                        pass

                    SF = SFDensity_age_cut(ds, (x, y, z), size, height, (age_range[0], age)) #Find SFRSD
                    print('SFRSD: ', SF)
                    sfrd.append(SF) #Add to arrays
                    sfsd = np.asarray(sfrd)
                    np.savetxt(sfFile+'age'+str(age)+'.txt', sfsd)
        sfrd.clear()
        del sfsd

#    return sfsd, gsd
