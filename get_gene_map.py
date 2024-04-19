#! /usr/bin/env python3
import re
import sys

regex = r"gene:([a-zA-Z0-9\_\-]+).*\s([=ckmnjeosxiypru])\s.*gene:([a-zA-Z0-9\_\-]+)"

if __name__ == "__main__":
    fn = sys.argv[1]
    map_length = 0
    map = {}
    with open(fn, "r") as f:
        for line in f.read().splitlines():
            try:
                match = re.match(regex, line).groups()
                print(match[0] + "\t" + match[2] + "\t" + match[1])
                map[match[0]] = match[1]
                map_length += 1
            except AttributeError:
                pass
