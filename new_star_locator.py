import numpy as np

ds = yt.load('DD1150/DD1150')
slc = yt.SlicePlot(ds, 'z', 'Density', center=[0.5,0.5,0.6], width=(10,'kpc') )
ax = list(slc.plots.values())[0].axes
ax.scatter(recent_stars_x.in_units('kpc'), recent_stars_y.in_units('kpc'), c = 'red',s = 0.5); slc.save('new_stars.png')
