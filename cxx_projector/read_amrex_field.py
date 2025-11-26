import numpy as np
import yt 
import sys 
import h5py

#simple script to read an array from an amrex plotfile, save it as a text file to be loaded by another program
#workaround for reading the array from the raw data in the plotfile itself
#the order of command line arguments is: plotfile_name field_name savefile_name 

pltfile_name = sys.argv[1]
field_name = sys.argv[2]
save_name = sys.argv[3]

ds = yt.load(pltfile_name) 
ad = ds.all_data()

print(field_name)
field_names = ["x", "y", "z", "dx", "dy", "dz", field_name]
fields = []

for this_name in field_names: 
    #Put lengths in code units since that's how stuff in the projector is defined
    if this_name == "x" or this_name == "dx": 
        field = ad[this_name] / (ds.domain_right_edge[0] - ds.domain_left_edge[0])
        print("x/dx: ", field.min(), field.max())

    elif this_name == "y" or this_name == "dy": 
        field = ad[this_name] / (ds.domain_right_edge[1] - ds.domain_left_edge[1]) 
        print("y/dy: ", field.min(), field.max())

    elif this_name == "z" or this_name == "dz": 
        field = ad[this_name] / (ds.domain_right_edge[2] - ds.domain_left_edge[2]) 
        print("z/dz: ", field.min(), field.max())
    else: 
        field = ad[this_name]
        print("fld: ", field.min(), field.max())

    if this_name == "x" or this_name == "y" or this_name == "z": 
        field += 0.5
        print("post add ", field.min(), field.max())
    fields.append(field) 
with h5py.File(save_name + ".h5", "w") as hf: 
    for i in range(len(field_names)): 
        hf.create_dataset(field_names[i], data=fields[i])
