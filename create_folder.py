import os 
import time

# Details for SSH connection to COSMA
##### CHANGE THIS
for i in range(27,78):
  file="flamingo_00"+str(i).zfill(2)
  print(file)
  os.system(f"mkdir /cosma8/data/do012/dc-yang9/data/Flamingo/L1000N1800/'{file}'/")

