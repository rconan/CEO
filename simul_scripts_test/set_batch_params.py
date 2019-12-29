import sys
import numpy as np
import math
#import ceo
#import scipy.io as sio
import os.path
import datetime
#import commands


# In[2]:
mag = float(sys.argv[1])
angle = float(sys.argv[2])
np.savez('batch_params',mag=mag, alpha_ps=angle*60.)
print "Saving param files with: mag %s, angle %s"%(mag, angle)

