#! /usr/bin/env python
import math
import sys

_, process_orgs, total_orgs, available_cpus = sys.argv


if __name__ == "__main__":
    alloc_cpus = int(math.floor(float(process_orgs) / float(total_orgs) * int(available_cpus)))
    print(max(1, alloc_cpus))
