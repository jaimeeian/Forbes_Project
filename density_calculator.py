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



