import unyt
import swiftsimio as sw
import numpy as np
import h5py
import tracemalloc
dir="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/flamingo_0077.hdf5"

soap_dir="../../../mnt/su3ctm/ludlow/Flamingo/L1000N0900/HYDRO_FIDUCIAL/SOAP-HBT/halo_properties_0077.hdf5"
tracemalloc.start()
#dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/colibre_with_SOAP_membership_0127.hdf5"
#soap_dir="../../../mnt/su3-pro/colibre/L0012N0094/THERMAL_AGN/SOAP/halo_properties_0127.hdf5"
data_h=sw.load(soap_dir)
input_ids=np.array(data_h.input_halos.halo_catalogue_index)
host_ids=np.array(data_h.soap.host_halo_index)
halo_ids=np.arange(0,len(host_ids),1)
mass=np.array(data_h.bound_subhalo.total_mass)   
mask=(mass>10000)*(host_ids==-1)
#main_ids=halo_ids[mask]
#mask=mask+np.isin(host_ids,main_ids)
sample_input_ids=input_ids[mask]
data_h=[]
input_ids=[]
host_ids=[]
halo_ids=[]
mass=[]
main_ids=[]
print("Done\n")
f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/halos.hdf5', 'a')
g=f["halos"]
g.create_dataset("input_ids",data=input_ids)
f.close()
'''
#mask = sw.mask(dir,spatial_only=False)
#boxsize = mask.metadata.boxsize
#load_region = [[0.0 * b, 0.5 * b] for b in boxsize]
#mask.constrain_spatial(load_region)
#mask.constrain_mask("gas", "group_nr_bound", -0.9*unyt.dimensionless,np.inf*unyt.dimensionless)
#mask.constrain_mask("stars", "group_nr_bound", -0.9*unyt.dimensionless,np.inf*unyt.dimensionless)
#mask.constrain_mask("dark_matter", "group_nr_bound", -0.9*unyt.dimensionless,np.inf*unyt.dimensionless)
#  mask.constrain_mask("stars", "group_nr_bound", -1.1*unyt.dimensionless, -0.9*unyt.dimensionless)
#  mask.constrain_mask("dark_matter", "group_nr_bound", -1.1*unyt.dimensionless, -0.9*unyt.dimensionless)
#data_h=sw.load(soap_dir)
#input_ids=np.array(data_h.input_halos.halo_catalogue_index)
#host_ids=np.array(data_h.soap.host_halo_index)
#data=sw.load(dir)#, mask=mask)
#print("start")
#print(tracemalloc.get_traced_memory())

#print(len(data.dark_matter.coordinates[0]))
#member_dm=np.array(data.dark_matter.group_nr_bound,dtype=np.int32)
#print(len(member_dm[member_dm!=-1]))
#print(tracemalloc.get_traced_memory())
#pmask=np.isin(member_dm,sample_input_ids)
#print(tracemalloc.get_traced_memory())
#member_dm=member_dm[pmask]
#print(len(member_dm))
#print(tracemalloc.get_traced_memory())
#Coord_dm=np.array(data.dark_matter.coordinates,dtype=np.float32)[pmask]
#potential_dm=np.array(data.dark_matter.potentials,dtype=np.float32)
#print("dm done")
#print(tracemalloc.get_traced_memory())
#member_star=np.array(data.stars.group_nr_bound,dtype=np.int32)
#pmask=np.isin(member_star,sample_input_ids)
#member_star=member_star[pmask]
#print(len(member_star))
#Coord_star=np.array(data.stars.coordinates,dtype=np.float32)[pmask]
#lum_z=np.array(data.stars.luminosities.GAMA_z,dtype=np.float32)[pmask]
#potential_star=np.array(data.stars.potentials,dtype=np.float32)
print("star done")
print(tracemalloc.get_traced_memory())
member_gas=np.array(data.gas.group_nr_bound,dtype=np.int32)

pmask=np.isin(member_gas,sample_input_ids)
#member_gas=member_gas[pmask]
#Coord_gas=np.array(data.gas.coordinates,dtype=np.float32)[pmask]
#potential_gas=np.array(data.gas.potentials,dtype=np.float32)

xray_lum_high=np.array(data.gas.xray_luminosities.erosita_high,dtype=np.float32)[pmask]
print("xlum!=0",len(xray_lum_high[xray_lum_high!=0]))
xray_lum_low=np.array(data.gas.xray_luminosities.erosita_low,dtype=np.float32)[pmask]
#a_lastagn=np.array(data.gas.last_agnfeedback_scale_factors,dtype=np.float32)[pmask]
print("gas done")
print(tracemalloc.get_traced_memory())
tracemalloc.stop()
print("Done\n")

f = h5py.File('/home/jyang/data/Flamingo/L1000N0900/cluster_particles.hdf5', 'a')
#dm = f.create_group("PartType1")
#dm.create_dataset("Coordinates", data=Coord_dm)
#dm.create_dataset("potentials", data=potential_dm)
#dm.create_dataset("member",data=member_dm)
#g= f.create_group("PartType0")
g=f["PartType0"]
del g["xray_lum_erosita_high"]
del g["xray_lum_erosita_low"]
#g.create_dataset("Coordinates", data=Coord_gas)
#g.create_dataset("potentials", data=potential_gas)
#g.create_dataset("xray_lum_erosita_high", data=xray_lum_high)
#g.create_dataset("xray_lum_erosita_low", data=xray_lum_low)
#g.create_dataset("a_agn", data=a_lastagn)
#g.create_dataset("member",data=member_gas)
#s= f.create_group("PartType2")
#s.create_dataset("Coordinates", data=Coord_star)
#s.create_dataset("potentials", data=potential_star)
#s.create_dataset("member",data=member_star)
#s.create_dataset("lum_z",data=lum_z)
f.close()
'''