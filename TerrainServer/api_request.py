import requests

url = "https://api.opentopodata.org/v1/srtm90m"
data = {
    "locations": "-43.5,172.5|27.6,1.98",
     "interpolation": "cubic",
}
response = requests.post(url, json=data)
