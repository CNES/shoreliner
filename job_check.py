import os
import matplotlib.pyplot as plt
import numpy as np
import yaml
from yaml.loader import SafeLoader

file_inputs = 'config.yaml'
inputs = yaml.load(open(os.path.join('./',file_inputs),'rb'),Loader = SafeLoader)

index_name = inputs['WaterlineIndex']
path = os.getcwd()

JOBS = False

#print(os.listdir())
c=0
N = []

for i in os.listdir():
  if i[:4] == 'poly':
    print(' ')
    print('----- '+ i + ' -----')
    tmp = os.path.join(path,i)
    if index_name in os.listdir(tmp):
      nb = len(os.listdir(os.path.join(tmp, index_name)))
      N.append(nb)
      print(str(nb)+' '+index_name+' in the folder '+i)
      if nb<=200:
        print('Launch of GEE image acquisition for ROI '+i)
        c+=1
        if JOBS:
          os.system('sbatch ./'+i+'/job.slurm')
    else:
      print('No '+index_name+' in the folder '+i)
      print('Launch of GEE image acquisition for ROI '+i)
      c+=1
      if JOBS:
        os.system('sbatch ./'+i+'/job.slurm')

print('In total, '+str(c)+' acquisition jobs have been launched')
plt.hist(N, np.arange(0,2500,10))
plt.savefig('hist_jobs.png')