#! /usr/bin/env python3
import re
import sys

PATTERN = r"gene:([a-zA-Z0-9\_\-]+).*\s([=ckmnjeosxiypru])\s.*gene:([a-zA-Z0-9\_\-]+)"

if __name__ == "__main__":
    fn = sys.argv[1]
    MAP_LENGTH = 0
    MAP = {}
    with open(fn, "r") as f:
        for line in f.read().splitlines():
            try:
                match = re.match(PATTERN, line).groups()
                print(match[0] + "\t" + match[2] + "\t" + match[1])
                MAP[match[0]] = match[1]
                MAP_LENGTH += 1
            except AttributeError:
                pass
