#!/usr/bin/env python

import sys

pn = sys.argv[1]
tim = sys.argv[2]
out = sys.argv[3]

with open(pn, 'r') as pn_f, open(tim, 'r') as tim_f, open(out, 'w') as out_f:
    toas_done = 0
    pns = pn_f.readlines()
    for line in tim_f.readlines():
        if "FORMAT" in line.strip() or "MODE" in line.strip() or toas_done >= len(pns) or line.strip()[0] == 'C':
            print(line.strip(), file=out_f)
        else:
            print(f"{line.strip()} -pn {pns[toas_done].strip()}", file=out_f)
            toas_done += 1

