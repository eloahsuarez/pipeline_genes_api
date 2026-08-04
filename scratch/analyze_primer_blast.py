import requests

url = "https://www.ncbi.nlm.nih.gov/tools/primer-blast/"
response = requests.get(url)
print(response.status_code)
html = response.text
with open("primer_blast.html", "w", encoding="utf-8") as f:
    f.write(html)
print("Saved to primer_blast.html")
