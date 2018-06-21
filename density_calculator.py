import numpy as np
import matplotlib.pyplot as plt
import yt


def cell_patch(data, x, y, size):

    ad = data.all_data()
        
    #Find the cells in the patch
    x_cells = ad['x'] - data.arr(0.5, 'code_length')

    y_cells = ad['y'] - data.arr(0.5, 'code_length')

#    z_cells = ad['z'] - data.arr(0.5, 'code_length')

    x_cells_pc = x_cells.in_units('pc')
    y_cells_pc = y_cells.in_units('pc')
#    z_cells_pc = z_cells.in_units('pc')
    xbounds_cells = np.logical_and(x_cells_pc > x - size, x_cells_pc < x + size)
    ybounds_cells = np.logical_and(y_cells_pc > y - size, y_cells_pc < y + size)
    cell_box = np.logical_and(xbounds_cells, ybounds_cells)

    print('total cells in x mask: ', np.sum(xbounds_cells))
    print('total cells in y mask: ', np.sum(ybounds_cells))
    print('total cells: ', np.sum(cell_box))
    return cell_box


def particle_patch(data, x, y, size):
    
    ad = data.all_data()

    #Find all the particles in the patch
    x_particles = ad['particle_position_x'] - data.arr(0.5, 'code_length')

    y_particles = ad['particle_position_y'] - data.arr(0.5, 'code_length')

    x_particles_pc = x_particles.in_units('pc')
    y_particles_pc = y_particles.in_units('pc')

    xbounds_particles = np.logical_and(x_particles_pc > x - size, x_particles_pc < x + size)
    ybounds_particles = np.logical_and(y_particles_pc > y - size, y_particles_pc < y + size)
    particle_box = np.logical_and(xbounds_particles, ybounds_particles)

    print('total particles in x mask: ', np.sum(xbounds_particles))
    print('total particles in y mask: ', np.sum(ybounds_particles))
    print('total particles: ', np.sum(particle_box))
    return particle_box



def GasDensity(ds, x, y, size):
    ad = ds.all_data()

    box = cell_patch(ds, x, y, size)

    cell_mass = ad['cell_mass']
    patch_cells = cell_mass[box]
    total_mass = np.sum(patch_cells)
    density = total_mass/(size**2)
    return density

def SFDensity(ds, center, size, max_age):
    ad = ds.all_data()
    x, y = center


    """
    new_stars = ad['creation_time'] >0
    xpos = ad['particle_position_x']
    ypos = ad['particle_position_y']
    zpos = ad['particle_position_z']
    new_stars_x = xpos[new_stars]
    new_stars_y = ypos[new_stars]
    new_stars_z = zpos[new_stars]
    """
    particle_mass = ad['particle_mass'].in_units('Msun')
    box = particle_patch(ds, x, y, size)
    
    ages = ad['creation_time']
    patched_ages = ages[box]
    new_stars = patched_ages > 0
    new_star_ages = patched_ages[new_stars]
    Myr_ages = new_star_ages.in_units('Myr')
    young = Myr_ages < max_age
#    print('total new stars: ', new_star_ages.shape)
 #   print('total stars in box: ', ages.shape)
  #  print('young stars mask: ', young.shape)
    stars = Myr_ages[young]
    star_mass = particle_mass[box]
    star_mass = star_mass[new_stars]
    star_mass = star_mass[young]
    sfrDen = star_mass/((size**2)*max_age)
    
    return sfrDen


def sigmaPlotter(ds, step, age):
    sfFile = input('Enter SFR file name: ')
    gFile = input('Enter Gas file name: ')
    g = np.array([])
    sf = np.array([])
    for x in np.arange(-300, 300, step):
        for y in np.arange(-300, 300, step):
            SF = SFDensity(ds, (x, y), 100, age)
            G = GasDensity(ds, x, y, 100)
            np.append(sf, SF)
            np.append(g, G)
    
    np.savetxt(sf, sfFile)
    np.savetxt(g, gFile)
    return sf, g
