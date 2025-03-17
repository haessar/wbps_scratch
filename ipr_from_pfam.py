import json

from bs4 import BeautifulSoup
from requests_html import HTMLSession


if __name__ == "__main__":
    URL = "https://www.ebi.ac.uk/interpro/entry/pfam/?search={}#table"
    TE_PFAM_PATH = "data/transposable_elements/transposable_element_pfams.json"
    session = HTMLSession()
    with open(TE_PFAM_PATH, "r") as f:
        te_pfams = json.load(f)
    for pfam in te_pfams.keys():
        resp = session.get(URL.format(pfam))
        resp.html.render()
        soup = BeautifulSoup(resp.content, 'html.parser')
        soup.find_all(lambda tag: tag.name == "a" and "IPR" in tag.text)
