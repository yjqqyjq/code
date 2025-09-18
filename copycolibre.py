import os 
import time

# Details for SSH connection to COSMA
##### CHANGE THIS
uname="dc-yang9"
password="1qaz2wsx#EDC"
fpath_cosma='/cosma8/data/dp004/flamingo/Runs/L2800N5040/HYDRO_FIDUCIAL'
fpath_hyades='/mnt/su3-pro/L2800N5040'


# Snapshots with SOAP
snapshots_all=[78] ##### CHANGE THIS

# Create directories if they don't exist
if not os.path.exists(f'{fpath_hyades}/'):
    try:
        os.makedirs(f'{fpath_hyades}/')
    except:
        pass

# Loop over snapshots
for snap in snapshots_all:
    # Download particle data files
    fpath=f'{fpath_cosma}/snapshots/flamingo_{str(snap).zfill(4)}' ##### CHANGE THIS
    print('Downloading particle data files:',fpath)
    os.system(f"/home/rwright/code/modulefiles/sshpass-1.10/bin/sshpass -p '{password}' rsync -vPrh -e 'ssh -c aes256-ctr' {uname}@login8.cosma.dur.ac.uk://{fpath} {fpath_hyades}/")
    time.sleep(5)

