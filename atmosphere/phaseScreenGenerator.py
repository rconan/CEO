import ceo
import json
from subprocess import call
import os
from datetime import datetime
import numpy as np

SEED = int(datetime.now().timestamp())

#for SEED in [int(datetime.now().timestamp()),]:
for k,L0 in enumerate(np.arange(5,95,5,dtype=np.float)):
    #L0   = 30.0
    filepath = "/home/ubuntu/DATA/"
    filename = filepath+"gmtAtmosphereL0%d_%d"%(L0,SEED)
    atm_prms = {'r0':0.15,'L0':L0,'L':26.0,'NXY_PUPIL':1300,'fov':20*ceo.constants.ARCMIN2RAD,
               'duration':15,'N_DURATION':80,'filename':filename+'.bin',
                'SEED':SEED}
    gpu_id = k%8
    with open(filename+'.json', 'w') as outfile:
        print(">>> "+outfile.name)
        json.dump(atm_prms, outfile, sort_keys = False, indent = 4, ensure_ascii=False)

    ceo.setDevice(gpu_id)
    atm = ceo.GmtAtmosphere(**atm_prms)

    #os.system("zip -0 -j "+filename+".zip "+filename+".json "+filename+".bin")
    #os.system("aws s3 cp "+filename+".zip s3://gmto.rconan/PhaseScreens/")
    #os.system("rm "+filename+".zip "+filename+".json "+filename+".bin")
