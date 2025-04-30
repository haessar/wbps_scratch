#! /usr/bin/env python3
import codecs
import json
import sys
import urllib.request

URL = "https://alphafold.ebi.ac.uk/api/prediction/{}"

if __name__ == "__main__":
    acc = sys.argv[1]
    response = urllib.request.urlopen(URL.format(acc))
    if response.status == 200:
        reader = codecs.getreader("utf-8")
        obj = json.load(reader(response))
        print(">" + acc)
        print(obj[0]["uniprotSequence"])
