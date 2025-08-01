#visualization hdf5
import numpy as np
import unyt
import functions as fn
import h5py
path="/Users/24756376/data/Flamingo/L1000N0900/"
id=-2795

mode=1
dm_ext=[]
g_ext=["xray_lum_erosita_low","T"]
s_ext=[]
main_id=fn.halo_ids[fn.halo_ids<=0]
mainarg=np.argwhere((fn.halo_ids==-id))
if dm_ext!=[]:
  extra_dm=[]
if g_ext!=[]:
  extra_g=[]
if s_ext!=[]:
  extra_s=[]
if mode=="sub":
  subarg=np.argwhere((fn.halo_ids>id)*(fn.halo_ids<id+1)+(fn.halo_ids==-id))
  center=fn.centers[fn.halo_ids<=0][id]
  N_dm_sub=fn.N_dm[(fn.halo_ids>id)*(fn.halo_ids<id+1)+(fn.halo_ids==-id)]
  N_g_sub=fn.N_g[(fn.halo_ids>id)*(fn.halo_ids<id+1)+(fn.halo_ids==-id)]
  N_s_sub=fn.N_s[(fn.halo_ids>id)*(fn.halo_ids<id+1)+(fn.halo_ids==-id)]
  member_dm=[]
  member_g=[]
  member_s=[]
  center_sub=fn.centers[(fn.halo_ids>id)*(fn.halo_ids<id+1)+(fn.halo_ids==-id)]

  for i in range(0,len(subarg)):
  
#print(center)
    
    member_dm.append(np.ones(N_dm_sub[i])*(i+1))
 
    member_g.append(np.ones(N_g_sub[i])*(i+1))
   
    member_s.append(np.ones(N_s_sub[i])*(i+1)) 

  member_dm=np.concatenate(member_dm,dtype=np.float16)
  member_g=np.concatenate(member_g,dtype=np.float16)
  member_s=np.concatenate(member_s,dtype=np.float16)

if mode==1:
  particle=fn.load_particles(path,id,dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext},mode="halo")
  center=fn.centers[fn.halo_ids==id]
  print("1")
  
if mode==2:
  particle=fn.load_particles(path,id,dm=1,g=1,s=1,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext},mode="cluster")
  center=center=fn.centers[fn.halo_ids<=0][int(id)]
  print("2")

#extra entry
Coord_dm=particle[0][0]
Coord_g=particle[1][0]

Coord_s=particle[2][0]
if dm_ext!=[]:    
      for j in range(0,len(dm_ext)):       
        extra_dm.append(particle[0][j+1])
if g_ext!=[]:    
      for j in range(0,len(g_ext)):    
        extra_g.append(particle[1][j+1])
       
if s_ext!=[]:    
      for j in range(0,len(s_ext)):
         extra_s.append(particle[2][j+1])

if mode =="halo":
  f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(f"{np.float32(id):.2f}")+'.hdf5', 'w')
else:
  f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(int(id))+'.hdf5', 'w')
dm = f.create_group("PartType1")
dm.create_dataset("Coordinates", data=Coord_dm-center)
if dm_ext!=[]:
  for i in range(0,len(dm_ext)):
    dm.create_dataset(dm_ext[i],data=extra_dm[i])
g = f.create_group("PartType0")
g.create_dataset("Coordinates", data=Coord_g-center)
if g_ext!=[]:
  for i in range(0,len(g_ext)):
    
    g.create_dataset(g_ext[i],data=extra_g[i])
st = f.create_group("PartType2")
st.create_dataset("Coordinates", data=Coord_s-center)
if s_ext!=[]:
  for i in range(0,len(s_ext)):
    st.create_dataset(s_ext[i],data=extra_s[i])
if mode=="sub":
  dm.create_dataset("member",data=member_dm)
  g.create_dataset("member",data=member_g)
  st.create_dataset("member",data=member_s)
f.close()