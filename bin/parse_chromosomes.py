#! /usr/bin/env python
from collections import Counter
from itertools import zip_longest
import os
import re
import sys

patterns = ["[0-9][0-9]", "[0-9]",]


def common_prefixes(li):
    threshold = len(li)
    prefix = []
    prefixes = []
    for chars in zip_longest(*li, fillvalue=''):
        char, count = Counter(chars).most_common(1)[0]
        if count == 1:
            break
        elif count < threshold:
            if prefix:
                prefixes.append((''.join(prefix), threshold))
            threshold = count
        prefix.append(char)
    if prefix:
        prefixes.append((''.join(prefix), threshold))
    return prefixes


def parse_chromosomes(gff):
    species, _ = os.path.splitext(os.path.basename(gff))
    regex = ""
    with open(gff, 'rt') as f:
        data = f.readlines()
    regions = [list(filter(None, l.split(" ")))[1] for l in data if l.startswith("##sequence-region")]
    chromosomes = [r for r in regions if "_" in r]
    if chromosomes:
        with_suffix = list(filter(re.compile(r".*_v\d+$").match, chromosomes))
        if len(with_suffix) > 1:
            prefix, freq = common_prefixes(with_suffix)[0]
        elif len(chromosomes) > 1:
            prefix, freq = common_prefixes(chromosomes)[0]
        else:
            prefix, freq = chromosomes[0].split("_")[0] + "_", 1
        if prefix and not (freq == 1 and len(regions) > 1):
            suffix = re.search(r"_v\d+", with_suffix[0]).group() if any(with_suffix) else ""
            r = re.compile(r"{}.*{}$".format(prefix, suffix))
            filt_chromosomes = list(filter(r.match, chromosomes))
            max_match = 0
            best_pattern = ""
            for pattern in patterns:
                r = re.compile(r"{}{}{}$".format(prefix, pattern, suffix))
                len_match = len(list(filter(r.match, filt_chromosomes)))
                if len_match > max_match:
                    max_match = len_match
                    best_pattern = pattern
            if not best_pattern:
                best_pattern = ".*"
            regex = "{}({}){}".format(prefix, best_pattern, suffix)
    return "\t".join(filter(None, (species, regex)))


if __name__ == "__main__":
    gff = sys.argv[1]
    print(parse_chromosomes(gff))
