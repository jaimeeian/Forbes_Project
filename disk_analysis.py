import yt
import numpy as np
import pdb


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
yt.add_field( "theta-velocity", function=_vtheta, units='km/s')
yt.add_field( "expectedSFR", function=_expectedSFR, units='msun/yr/pc**3')


def define_patch( data, location, width, height, max_temp):
    ''' Pick out cells in dataset "data" centered on the (x,y,z) tuple "location", with a rectangular prism of width "width", height "height", and only return cells with temperatures less than max_temp. Returns a mask of True/False values with the same shape as all gas-based fields in "data".'''
    x = (data['x'] - data.ds.arr(0.5, 'code_length')).in_units('kpc')
    y = (data['y'] - data.ds.arr(0.5, 'code_length')).in_units('kpc')
    z = (data['z'] - data.ds.arr(0.6, 'code_length')).in_units('kpc')
    loc = (location[0].in_units('kpc'), location[1].in_units('kpc'), location[2].in_units('kpc'))
    wid = width.in_units('kpc')
    hei = height.in_units('kpc')
    mask = np.logical_and(x > loc[0]-wid/2, x<loc[0]+wid/2)
    mask = np.logical_and(mask, y > loc[1]-wid/2)
    mask = np.logical_and(mask, y < loc[1]+wid/2)
    mask = np.logical_and(mask, z > loc[2]-hei/2)
    mask = np.logical_and(mask, z < loc[2]+hei/2)
    mask = np.logical_and(mask, data['temperature']<max_temp)
    return mask

def test_patch(data):
    ''' Define a patch in the middle of the galaxy, and add up all the star formation inside, estimated based on the gas-only quantities. Just testing if the patch function and the expectedSFR field are working as expected.'''
    a = define_patch( data, (yt.YTQuantity(0.1,'kpc'),yt.YTQuantity(0.1,'kpc'),yt.YTQuantity(0.0,'kpc')), yt.YTQuantity(3, 'kpc'), yt.YTQuantity(1.0, 'kpc'), yt.YTQuantity(1.0e4,'K'))
    print(np.sum((data['expectedSFR'][a] * data['cell_volume'][a]).in_units('msun/yr')))

if __name__=='__main__':
    ds = yt.load('DD1150/DD1150')
    ad = ds.all_data()
    test_patch(ad)
