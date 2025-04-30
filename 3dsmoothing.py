from scipy.ndimage import gaussian_filter
import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
folder = "/home/jyang/plot/Flamingo/L0200N0720/halo_2d_lum_rotate/66_kde0/"
img = Image.open(folder+" 0.7208, 0.9334, 0.9379 .png")#.convert("L")
img_array = np.array(img)
img_smooth = gaussian_filter(img_array, sigma=1,mode='constant')
plt.imshow(img_smooth)

