import numpy as np
import matplotlib.pyplot as plt
import yt

def patch(data, x, y):
#    ds = yt.load(data)
    ad = data.all_data()
    
#    x = x - 0.5
#   y = y - 0.5

    #xpos = ad['particle_position_y] - data.arr(0.5, 'code_length')
#   ypos = ad['particle_position_y'] - data.arr(0.5, 'code_length')
    
#    slc = yt.SlicePlot(data, 'z', 'Density', center=[x, y, 0.6], width=(10,'kpc') )

#    ax = list(slc.plots.values())[0].axes

#    ax.scatter(xpos, ypos)
    
#    dist = np.sqrt(x**2 + y**2)
    
    x_cells = ad['x'] - data.arr(0.5, 'code_length')
#    print(xs.min())
    y_cells = ad['y'] - data.arr(0.5, 'code_length')
#    z_cells = ad['z'] - data.arr(0.6, 'code_length')
    x_cells_pc = x_cells.in_units('pc')
    y_cells_pc = y_cells.in_units('pc')
#    z_cells_pc = z_cells.in_units('pc')
    xbounds = np.logical_and(x_cells_pc > x - 100, x_cells_pc < x + 100)
    ybounds = np.logical_and(y_cells_pc > y - 100, y_cells_pc < y + 100)
    box = np.logical_and(xbounds, ybounds)
    print('total cells in x mask: ', np.sum(xbounds))
    print('total cells in y mask: ', np.sum(ybounds))
    print('overlap: ', np.sum(box))
    return box
#    ybox = np.logical_and(ybound < y + 100, ybound > y - 100)
 #   box = np.logical_and(xbox, ybox)
#    return box
