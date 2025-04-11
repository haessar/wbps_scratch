#! /usr/bin/env python3
import argparse
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
    parser = argparse.ArgumentParser()
    parser.add_argument('output_dir')
    parser.add_argument('--species-list', '-s', type=str, default=None)
    parser.add_argument('--softmasked', "-m", action="store_true")
    args = parser.parse_args()

    species_list = []
    if args.species_list:
        with open(args.species_list, "r") as f:
            species_list = f.read().splitlines()
    FA_SUFFIX = "genomic_softmasked.fa.gz" if args.softmasked else "genomic.fa.gz"
    for r1 in soup_resp(base_url).find("table").find_all("tr"):
        a1 = r1.find('a')
        if a1 and "_" in a1.text:
            species = a1.text
            for r2 in soup_resp(base_url + species).find("table").find_all("tr"):
                a2 = r2.find('a')
                if a2 and r2.find('img').attrs.get('src') == "/icons/folder.gif":
                    acc = a2.text
                    if species_list and not any([s for s in species_list if s.startswith(species.strip("/")) and s.endswith(acc.strip("/").lower())]):
                        continue
                    for r3 in soup_resp(base_url + species + acc).find('table').find_all("tr"):
                        a3 = r3.find('a')
                        if a3:
                            if a3.text.endswith(FA_SUFFIX) or a3.text.endswith("annotations.gff3.gz"):
                                if a3.text not in os.listdir(args.output_dir):
                                    print(a3.text)
                                    download_file(base_url + species + acc + a3.text, datadir=args.output_dir)
