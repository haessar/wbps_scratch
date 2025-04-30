from collections import defaultdict
import json
import re

import gffutils
import requests

INTERPRO_API = "https://www.ebi.ac.uk/interpro/api/entry/interpro/{}"


def get_pfams_for_feature_accession(f: gffutils.Feature, ipa_map: defaultdict[list]):
    pfams = []
    if f.featuretype == "mRNA":
        if "info" in f.attributes:
            for ipa in re.findall(r"InterPro accession:(IPR[0-9]*)\s", dict(f.attributes)["info"][0]):
                if ipa not in ipa_map:
                    resp = requests.get(INTERPRO_API.format(ipa), timeout=300)
                    try:
                        pfams = list(json.loads(resp.content)['metadata']['member_databases'].get('pfam', {}).keys())
                        ipa_map[ipa].extend(pfams)
                    except KeyError:
                        pass
                    except json.decoder.JSONDecodeError:
                        print(f)
                        raise
                else:
                    pfams = ipa_map[ipa]
    return pfams
