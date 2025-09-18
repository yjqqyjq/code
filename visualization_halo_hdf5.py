#visualization hdf5
import numpy as np
import unyt
import functions as fn
import h5py
from tqdm import tqdm
path="/Users/24756376/data/Flamingo/L1000N0900/"

#[1815 5673 6402  836 2113 9249 9065 4356 2852 4957]
#[2.08984375 1.49609375 1.44921875 2.11328125 1.859375   1.32421875
# 1.27148438 1.5859375  1.75585938 1.56640625]

def vis(id):
  dms=1
  gs=1
  ss=0
  mode="halor"
  radius=1#in r200
  dm_ext=["Velocities"]
  g_ext=["Velocities"]#["xray_lum_erosita_low","T"]
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
    subarg=np.argwhere((fn.halo_ids>-id)*(fn.halo_ids<-id+1)+(fn.halo_ids==id))
    center=fn.centers[fn.halo_ids<=0][-id]
    N_dm_sub=fn.N_dm[(fn.halo_ids>-id)*(fn.halo_ids<-id+1)+(fn.halo_ids==id)]
    N_g_sub=fn.N_g[(fn.halo_ids>-id)*(fn.halo_ids<-id+1)+(fn.halo_ids==id)]
    N_s_sub=fn.N_s[(fn.halo_ids>-id)*(fn.halo_ids<-id+1)+(fn.halo_ids==id)]
    member_dm=[]
    member_g=[]
    member_s=[]
    center_sub=fn.centers[(fn.halo_ids>-id)*(fn.halo_ids<-id+1)+(fn.halo_ids==id)]

    for i in range(0,len(subarg)):
  
#print(center)
    
      member_dm.append(np.ones(N_dm_sub[i])*(i+1))
 
      member_g.append(np.ones(N_g_sub[i])*(i+1))
   
      member_s.append(np.ones(N_s_sub[i])*(i+1)) 

    member_dm=np.concatenate(member_dm,dtype=np.float16)
    member_g=np.concatenate(member_g,dtype=np.float16)
    member_s=np.concatenate(member_s,dtype=np.float16)

  if mode=="halo" or "halor":
    particle=fn.load_particles(path,id,dm=dms,g=gs,s=ss,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext},mode="halo")
    particles=[]
    if mode=="halor":
      
      if dms==1:
        data=[]
        for i in range(0,len(particle[0])):
         
          r=particle[0][0][:,0]**2+particle[0][0][:,1]**2+particle[0][0][:,2]**2
          r=np.sqrt(r)
          data.append(particle[0][i][r<fn.r100[fn.halo_ids<=0][-id]])
        particles.append(data)
      if gs==1:
        data=[]
    
        for i in range(0,len(particle[1])):
         
          r=particle[1][0][:,0]**2+particle[1][0][:,1]**2+particle[1][0][:,2]**2
          r=np.sqrt(r)
          data.append(particle[1][i][r<fn.r100[fn.halo_ids<=0][-id]])
        particles.append(data)
      if ss==1:
        data=[]
      
        for i in range(0,len(particle[2])):
         
          r=particle[2][0][:,0]**2+particle[2][0][:,1]**2+particle[2][0][:,2]**2
          r=np.sqrt(r)
          data.append(particle[2][i][r<fn.r100[fn.halo_ids<=0][-id]])
        particles.append(data)
      particle=particles
    print("1")
  elif mode=="region":
    particle=fn.load_regions(path,id,radius*fn.r100[fn.halo_ids<=0][-id],dm=dms,g=gs,s=ss,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext},mode="all")
  elif mode=="unbound":
    particle=fn.load_regions(path,id,radius*fn.r100[fn.halo_ids<=0][-id],dm=dms,g=gs,s=ss,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext},mode="unbound")
  else:
    particle=fn.load_particles(path,id,dm=dms,g=gs,s=ss,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext},mode="cluster")
  
    print("1")

#extra entry
  if dms==1:
    Coord_dm=particle[0][0]
    if dm_ext!=[]:    
        for j in range(0,len(dm_ext)):   
          if dm_ext[j]=="Velocities":    
            extra_dm.append(np.sum(particle[0][j+1]*Coord_dm,axis=1))#radial velocity
          else:
            extra_dm.append(particle[0][j+1])
  if gs==1:
    Coord_g=particle[1][0]
    if g_ext!=[]:    
        for j in range(0,len(g_ext)):    
          if g_ext[j]=="Velocities":    
            extra_g.append(np.sum(particle[1][j+1]*Coord_g,axis=1))#radial velocity
          else:
            extra_g.append(particle[1][j+1])
  if ss==1:
    Coord_s=particle[2][0]


       
    if s_ext!=[]:    
     
        for j in range(0,len(s_ext)):
          if s_ext[j]=="Velocities":    
            extra_s.append(np.sum(particle[2][j+1]*Coord_s,axis=1))#radial velocity
          else:
            extra_s.append(particle[2][j+1])

  if mode =="halo":
    f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(f"{np.float32(id):.2f}")+'.hdf5', 'w')
  elif mode=="region":
    f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(int(id))+'_'+str(np.float16(radius))+'_r100.hdf5', 'w') 
  elif mode=="unbound":
    f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(int(id))+'_'+str(np.float16(radius))+'_r100_unbound.hdf5', 'w')
  else:
    f = h5py.File('/Users/24756376/data/Flamingo/L1000N0900/halos/'+str(int(id))+'.hdf5', 'w')
  if dms==1:
    dm = f.create_group("PartType1")
    dm.create_dataset("Coordinates", data=Coord_dm)
    if dm_ext!=[]:
      for i in range(0,len(dm_ext)):
        dm.create_dataset(dm_ext[i],data=extra_dm[i])
  if gs==1:
    g = f.create_group("PartType0")
    g.create_dataset("Coordinates", data=Coord_g)
    if g_ext!=[]:
      for i in range(0,len(g_ext)):
    
        g.create_dataset(g_ext[i],data=extra_g[i])
  if ss==1:
    st = f.create_group("PartType2")
    st.create_dataset("Coordinates", data=Coord_s)
    if s_ext!=[]:
      for i in range(0,len(s_ext)):
        st.create_dataset(s_ext[i],data=extra_s[i])
  if mode=="sub":
    if dms==1:
      dm.create_dataset("member",data=member_dm)
    if gs==1:
      g.create_dataset("member",data=member_g)
    if ss==1:
      st.create_dataset("member",data=member_s)
  f.close()

for i in tqdm(range(0,101)):
  vis(-i)