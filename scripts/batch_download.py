import requests
import time
import sys
import os

def download_uniprot_batches(query, output_dir, base_name, batch_size=500):
    # Ensure the directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": query, "format": "fasta", "size": batch_size}
    
    url = base_url
    batch_count = 1
    
    while url:
        response = requests.get(url, params=params if batch_count == 1 else None)
        
        if response.status_code == 200:
            filename = os.path.join(output_dir, f"{base_name}_part_{batch_count}.fasta")
            with open(filename, "w") as f:
                f.write(response.text)
            print(f"Saved: {filename}")
            
            if "Link" in response.headers:
                links = response.headers["Link"].split(",")
                next_link = [l for l in links if 'rel="next"' in l]
                if next_link:
                    url = next_link[0].split(";")[0].strip("<>")
                    batch_count += 1
                    params = None
                    time.sleep(1)
                else:
                    url = None
            else:
                url = None
        else:
            print(f"Error: {response.status_code}")
            break

if __name__ == "__main__":
    # Usage: python download_uniprot.py "query" "output_directory" "filename_prefix"
    query_arg = sys.argv[1]
    dir_arg = sys.argv[2]
    prefix_arg = sys.argv[3]
    download_uniprot_batches(query_arg, dir_arg, prefix_arg)