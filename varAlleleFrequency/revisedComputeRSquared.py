#!/usr/bin/env python

def computersquared():
    from scipy import stats
    import numpy as np
    x = np.random.random(10)
    y = np.random.random(10)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

    return r_value**2

if __name__ == '__main__':
    r2 = computersquared()
    print r2
