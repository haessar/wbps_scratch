#! /usr/bin/env python3
import codecs
import json
import sys
import urllib.request

url = "https://alphafold.ebi.ac.uk/api/prediction/{}"

if __name__ == "__main__":
    acc = sys.argv[1]
    response = urllib.request.urlopen(url.format(acc))
    if response.status == 200:
        reader = codecs.getreader("utf-8")
        obj = json.load(reader(response))
        print(">" + obj[0]["gene"])
        print(obj[0]["uniprotSequence"])
