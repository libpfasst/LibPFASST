import h5py
import getProblemDefinition
filename = "test.h5"
problem = getProblemDefinition.get()
file = h5py.File(filename,'w')
for key1 in problem.keys():
    group1 = file.create_group(key1)
    for key2 in problem[key1].keys():
        group2 = group1.create_group(key2)
        for key3 in problem[key1][key2].keys():
            dset     = group2.create_dataset(key3, problem[key1][key2][key3][2], problem[key1][key2][key3][0])
            dset[...] = problem[key1][key2][key3][1]
file.close()
