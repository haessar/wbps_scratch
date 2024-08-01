import re


def extract_accessions_from_transcript(tran):
    try:
        info = tran.attributes.get("info")[0]
    except TypeError:
        print(f"No InterPro accessions found for transcript {tran.id}")
        return
    print(info)
    for acc in info.split("\n"):
        yield re.search(r'IPR\d+', acc).group(), acc.split("description:")[1]
