import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

df = pd.read_csv("../../test_terrain_data.csv")

# make target grid
lon_lin = np.linspace(df.longitude.min(), df.longitude.max(), 100)
lat_lin = np.linspace(df.latitude.min(), df.latitude.max(), 100)
X, Y = np.meshgrid(lon_lin, lat_lin)

# interpolate elevations to grid
Z = griddata(
    (df['longitude'], df['latitude']),
    df['elevation'],
    (X, Y),
    method='cubic'
)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X, Y, Z)
plt.show()
