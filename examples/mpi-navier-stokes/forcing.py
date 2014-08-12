
import numpy as np
import h5py

forcing = {}
for attr in [ 'kx', 'ky', 'kz', 'ax', 'ay', 'az', 'px', 'py', 'pz', 'pt', 'ft' ]:
    forcing[attr] = []

for kx in range(4):
    for ky in range(4):
        for kz in range(4):
            if kx + ky + kz > 4:
                continue

            r = np.random.rand(8)

            forcing['kx'].append(kx)
            forcing['ky'].append(ky)
            forcing['kz'].append(kz)
            forcing['ax'].append(r[0])
            forcing['ay'].append(r[1])
            forcing['az'].append(r[2])
            forcing['px'].append(2*np.pi*r[3])
            forcing['py'].append(2*np.pi*r[4])
            forcing['pz'].append(2*np.pi*r[5])
            forcing['ft'].append((1 + r[6])*np.pi)
            forcing['pt'].append(2*np.pi*r[7])

h5 = h5py.File('forcing.h5', 'w')
for attr in [ 'kx', 'ky', 'kz', 'ax', 'ay', 'az', 'px', 'py', 'pz', 'pt', 'ft' ]:
    h5.create_dataset(attr, data=np.asarray(forcing[attr], np.float64))
h5.close()
