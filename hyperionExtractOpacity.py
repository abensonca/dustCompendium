# -*- coding: iso-8859-15 -*-
import sys
import numpy as np
import h5py
import h5py.highlevel
from hyperion.dust import SphericalDust

# Get arguments.
if ( len(sys.argv) != 2 ):
    sys.exit("Usage: hyperionExtractOpacity.py <dustFile>")
dustFile = sys.argv[ 1]  # Dust file name.

# Read dust data and determine the opacity to extinction at 0.55�m (V band).
dust  = SphericalDust()
dust.read(dustFile)
kappaV = dust.optical_properties.interp_chi_wav(0.55)

print kappaV
