import numpy as np
import matplotlib.pyplot as plt
import yt


def cell_patch(data, x, y, size, height):

    ad = data.all_data()
    #covert inputs into YT Quantities
    x = yt.YTQuantity(x, 'pc')
    y = yt.YTQuantity(y, 'pc')
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')
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
    return cell_box


def particle_patch(data, x, y, size, height):
    
    ad = data.all_data()
    
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



def GasDensity(ds, x, y, size, height):
    ad = ds.all_data()
    #Create patch of cells
    box = cell_patch(ds, x, y, size, height)
    #Make inputs into YT Quantites
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')
    #Find cells in patch
    cell_mass = ad['cell_mass'].in_units('Msun')
    patch_cells = cell_mass[box]
    #Create Temperature Cut
    cell_temp = ad['temperature'].in_units('K')
    cell_temp_patch = cell_temp[box]
    cold_gas = cell_temp_patch < 10**4
    #Add up mass of cells in patch under temp cut
    cold_gas_mass = patch_cells[cold_gas]
    total_mass = np.sum(cold_gas_mass)
    #Calculate density (Divide by area)
    density = total_mass/(size**2)
    return density

def SFDensity(ds, center, size, height, max_age):
    ad = ds.all_data()
    #define center, establish YT Quantities
    x, y = center
    size = yt.YTQuantity(size,'pc')
    height = yt.YTQuantity(height, 'pc')
    max_age = yt.YTQuantity(max_age, 'Myr')
    
    #Create array of masses & ages of all particles
    particle_mass = ad['particle_mass'].in_units('Msun')
    ages = ad['creation_time']
    #Create spatial patch
    box = particle_patch(ds, x, y, size, height)
    #Select new star particles in patch
    patched_ages = ages[box]
    new_stars = patched_ages > 0
    new_star_ages = patched_ages[new_stars]
    new_star_ages = ds.current_time - new_star_ages
    Myr_ages = new_star_ages.in_units('Myr')

    young = Myr_ages < max_age
    stars = Myr_ages[young]
    #Get masses of new stars
    star_mass = particle_mass[box]
    new_star_mass = star_mass[new_stars]
    young_star_mass = new_star_mass[young]
    #Find total mass of new stars, then surface density
    patched_stars = np.sum(young_star_mass)
    sfrDen = patched_stars/((size**2)*max_age)
    sfrDen = sfrDen.in_units('Msun/yr/kpc**2')
    return sfrDen


def sigmaCalc(ds, step, age, spatial_range):
    #Give in output file names
    sfFile = input('Enter SFR file name: ')
    gFile = input('Enter Gas file name: ')
    #Read in limits 
    low, high = spatial_range
    #Create empty lists to hold gas SD and SFR SD values
    gd = []
    sfrd = []

    for x in np.arange(low, high, step): #find x of center for each patch
        for y in np.arange(low, high, step): #find y of center
            print('x: ', x)
            print('y: ', y)
            SF = SFDensity(ds, (x, y), 100, 500, age) #Find SFRSD
            G = GasDensity(ds, x, y, 100, 500) #Find GSD
            print('SFRSD: ', SF)
            print('GSD: ', G)
            sfrd.append(SF) #Add to arrays
            gd.append(G) #Add to arrays
    
            sfsd = np.asarray(sfrd)
            gsd = np.asarray(gd)
            np.savetxt(sfFile, sfsd)
            np.savetxt(gFile, gsd)

    return sfsd, gsd



def expectedSFR_calculator(ds, x, y, size, height):
    ad = ds.all_data()
    #Create patch of cells
    box = cell_patch(ds, x, y, size, height)
    #Make inputs into YT Quantites
    size = yt.YTQuantity(size, 'pc')
    height = yt.YTQuantity(height, 'pc')
    #Find cells in patch
    cell_sfr = ad['expectedSFR'].in_units('Msun/yr/kpc**3')
    patch_cells = cell_sfr[box]
    v_cells = ad['cell_volume']
    v_cells_patched = v_cells[box]
        #Create Temperature Cut
#    cell_temp = ad['temperature'].in_units('K')
#    cell_temp_patch = cell_temp[box]
#    cold_gas = cell_temp_patch < 10**4
    #Add up mass of cells in patch under temp cut
#    cold_gas_mass = patch_cells[cold_gas]
    cell_sfr = patch_cells*v_cells_patched
    total_sfr = np.sum(cell_sfr)
    density = (total_sfr/size**2).in_units('Msun/yr/kpc**2')
    return density

def expSFRSD(ds, step, spatial_range):
    #Give in output file names
    expsfFile = input('Enter expected SFR file name: ')
#    gFile = input('Enter GSD file name: ')
    #Read in limits 
    low, high = spatial_range
    #Create empty lists to hold gas SD and SFR SD values

    sfrd = []
#   gd = []
    for x in np.arange(low, high, step): #find x of center for each patch
        for y in np.arange(low, high, step): #find y of center
            print('x: ', x)
            print('y: ', y)
            SF = expectedSFR_calculator(ds, x, y, 100, 500) #Find SFRSD
#            G = GasDensity(ds, x, y, 100, 500) #Find GasSD
            print('SFRSD: ', SF)
#            print('GSD: ', G)
            sfrd.append(SF) #Add to arrays
#            gd.append(G)
    
    sfsd = np.asarray(sfrd)
#    gsd = np.asarray(gd)
    np.savetxt(expsfFile, sfsd)
#    np.savetxt(gFile, gsd)
    return sfsd


def sigmaPlotter(g, sfr, esf, age_in_yr, size_in_kpc):

    upper_limit = np.zeros(len(sfr))
    for i in range(len(sfr)):
        if sfr[i] == 0:
            upper_limit[i] = 50./age_in_yr/(size_in_kpc**2)
            
    plt.loglog(g, sfr, 's', color = 'r', label = 'Simulation')
    plt.loglog(g, esf, 'o', label = 'Analytic')
    plt.loglog(g, upper_limit, 'v')
    plt.ylabel('$\dot{\Sigma}_{SF}$ ($M_\odot  kpc^{-2} yr^{-1}$)')
    plt.xlabel('$\Sigma (M_\odot  pc^{-2})$')
    
    plt.legend(loc = 'upper left')
    plt.show()
