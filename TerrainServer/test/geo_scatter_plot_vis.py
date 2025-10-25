from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv("../../test_terrain_data.csv")

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(df['longitude'], df['latitude'], df['elevation'])
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')
ax.set_zlabel('Elevation (m)')
plt.show()
