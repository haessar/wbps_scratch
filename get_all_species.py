#! /usr/bin/env python3
from bs4 import BeautifulSoup
import os
import os.path
import requests
import shutil
import sys

base_url = "https://ftp.ebi.ac.uk/pub/databases/wormbase/parasite/releases/current/species/"


def download_file(url, datadir):
    local_filename = url.split('/')[-1]
    with requests.get(url, stream=True) as r:
        with open(os.path.join(datadir, local_filename), 'wb') as f:
            shutil.copyfileobj(r.raw, f)
    return local_filename


def soup_resp(url):
    resp = requests.get(url)
    if resp.status_code == 200:
        return BeautifulSoup(resp.text)
    else:
        raise Exception("failed to request url {} with status code {}".format(url, resp.status_code))


if __name__ == "__main__":
    output_dir = sys.argv[1]
    for r1 in soup_resp(base_url).find("table").find_all("tr"):
        a1 = r1.find('a')
        if a1 and "_" in a1.text:
            species = a1.text
            for r2 in soup_resp(base_url + species).find("table").find_all("tr"):
                a2 = r2.find('a')
                if a2 and r2.find('img').attrs.get('src') == "/icons/folder.gif":
                    acc = a2.text
                    for r3 in soup_resp(base_url + species + acc).find('table').find_all("tr"):
                        a3 = r3.find('a')
                        if a3:
                            if a3.text.endswith("genomic.fa.gz") or a3.text.endswith("annotations.gff3.gz"):
                                if a3.text not in os.listdir(output_dir):
                                    print(a3.text)
                                    download_file(base_url + species + acc + a3.text, datadir=output_dir)
