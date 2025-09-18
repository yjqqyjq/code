#check density
import numpy as np
import unyt
import functions as fn
import h5py
from tqdm import tqdm
path="/Users/24756376/data/Flamingo/L1000N0900/"
dm_ext=[]
g_ext=[]#["xray_lum_erosita_low","T"]
s_ext=[]
r200=fn.r200[fn.r200>0][0:10000]
Rho=np.zeros(10000)
for i in tqdm(range(0,10000)):
    particle=fn.load_regions(path,-i,r200[i],dm=1,g=1,s=0,coordinate=1,extra_entry={"dm":dm_ext,"gas":g_ext,"stars":s_ext})
    M_dm=len(particle[0][0])*45.2#m in 10^9 Msun
    M_g=len(particle[1][0])*8.56
    M200=M_dm+M_g
    Rho[i]=M200#critical density in Msun/Mpc^3
print(Rho/4*3/r200**3/127.1/200/np.pi)