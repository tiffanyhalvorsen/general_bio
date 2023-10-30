#!/usr/bin/env python3

import sys
import re

fastq = sys.argv[1]

reads = 0
lines = 0
with open(fastq, "r") as file:
	for line in file:
		lines += 1
		for match in re.findall(r"^@.+", line):
			reads += 1

print(lines, reads)
		
