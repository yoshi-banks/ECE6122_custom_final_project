import folium
import pandas as pd

df = pd.read_csv("../../test_terrain_data.csv")

# Center on first point
m = folium.Map(location=[df.latitude.mean(), df.longitude.mean()], zoom_start=12)

for _, row in df.iterrows():
    folium.CircleMarker(location=[row.latitude, row.longitude],
                        radius=4,
                        popup=f"Elev: {row.elevation}m").add_to(m)

m.save("map.html")
