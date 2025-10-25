import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

df = pd.read_csv("../../test_terrain_data.csv")

# reshape to 2D grid
lat = np.sort(df['latitude'].unique())
lon = np.sort(df['longitude'].unique())
Z = df.pivot(index='latitude', columns='longitude', values='elevation').values

X, Y = np.meshgrid(lon, lat)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X, Y, Z)      # <-- surface!
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")
ax.set_zlabel("Elevation (m)")
plt.show()
