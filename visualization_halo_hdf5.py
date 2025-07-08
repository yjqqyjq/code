#visualization hdf5
import numpy as np
import unyt
import functions as fn
import h5py
path="/Users/24756376/data/Flamingo/L1000N0900/"
id=10
main_id=fn.halo_ids[fn.halo_ids<=0]
subarg=np.argwhere((fn.halo_ids>id)*(fn.halo_ids<id+1))
mainarg=np.argwhere((fn.halo_ids==-id))
center=fn.centers[fn.halo_ids<=0][id]
Coord_dm=[]
Coord_g=[]
Coord_s=[]
member_dm=[]
member_g=[]
member_s=[]
for i in range(0,len(subarg)):
  
#print(center)
  particle=fn.load_particles(path,fn.halo_ids[subarg[i]],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
  Coord_dm.append(particle[0][0])
  member_dm.append(np.ones(len(particle[0][0]))*(i+1)*0.1+0.1)
  Coord_g.append(particle[1][0])
  member_g.append(np.ones(len(particle[1][0]))*(i+1)*0.1+0.1)
  Coord_s.append(particle[2][0])
  member_s.append(np.ones(len(particle[2][0]))*(i+1)*0.1+0.1)
main_p=fn.load_particles(path,fn.halo_ids[mainarg],dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":[],"gas":[],"stars":[]},mode="halo")
Coord_dm.append(main_p[0][0])
Coord_g.append(main_p[1][0])
Coord_s.append(main_p[2][0])
member_dm.append(np.ones(len((main_p[0][0])))*0.1)
member_g.append(np.ones(len((main_p[1][0])))*0.1)
member_s.append(np.ones(len((main_p[2][0])))*0.1)
Coord_dm=np.concatenate(Coord_dm)
Coord_g=np.concatenate(Coord_g)
Coord_s=np.concatenate(Coord_s)
member_dm=np.concatenate(member_dm)
member_g=np.concatenate(member_g)
member_s=np.concatenate(member_s)


f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(int(id))+'.hdf5', 'w')
dm = f.create_group("PartType1")
dm.create_dataset("Coordinates", data=Coord_dm-center)
dm.create_dataset("member",data=member_dm)
g = f.create_group("PartType0")
g.create_dataset("Coordinates", data=Coord_g-center)
g.create_dataset("member",data=member_g)
st = f.create_group("PartType2")
st.create_dataset("Coordinates", data=Coord_s-center)
st.create_dataset("member",data=member_s)

f.close()
#print(Coord_dm)
