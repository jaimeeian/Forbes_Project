import numpy as np

def findSFR(ds, bin):
    ad = ds.all_data()
    new_stars = ad['creation_time']>0
    sf = ad['creation_time'][new_stars]
    sfr = np.array([])
    for i in np.arange(0, len(sf), bin):
        age_bin = (sf > i) & (sf < i+bin)
        stars = sf[age_bin]
        np.append(sfr, len(stars))
    
    return sfr
