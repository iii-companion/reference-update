#! /usr/bin/env python3
from collections import Counter
from itertools import groupby, zip_longest
from operator import itemgetter
import os
import re
import sys

patterns = [
    "([0-9][0-9])",
    "([0-9][0-9]){0,1}(MIT|API){0,1}",
    "([0-9])",
    "([MDCLXVI]+[a-z]?)"
]

class ChromosomeParser:
    def __init__(self, gff) -> None:
        self.species = ''
        self.prefix = ''
        self.pattern = ''
        self.suffix = ''
        self.sequence_regions = []
        self.chr_candidates = []
        self._parse(gff)
    
    @property
    def regex(self) -> str:
        return "{}{}{}".format(self.prefix, self.pattern, self.suffix)
    
    def __repr__(self) -> str:
        return "\t".join(filter(None, (self.species, self.regex)))
    
    def _parse(self, gff) -> None:
        self.species, _ = os.path.splitext(os.path.basename(gff))
        with open(gff, 'rt') as f:
            data = f.readlines()
        self.sequence_regions = [list(filter(None, l.split(" ")))[1] for l in data if l.startswith("##sequence-region")]
        self.chr_candidates = [r for r in self.sequence_regions if {".", "_"}.intersection(r)]
    
    def generate_regex(self) -> None:
        self.has_suffix()
        prefix_pool = self._iter_prefix_pool()
        while self.pattern == "":
            max_match = 0
            try:
                prefix, freq = next(prefix_pool)
            except StopIteration:
                print(self.species)
                return
                # self.pattern = ".*"
                # self.prefix, _, _ = next(self._iter_prefix_pool())
                # raise
            if prefix and not (freq == 1 and len(self.sequence_regions) > 1):
                r = re.compile(r"{}.*{}$".format(prefix, self.suffix), re.IGNORECASE)
                filt_candidates = list(filter(r.match, self.chr_candidates))
                for pattern in patterns:
                    r = re.compile(r"{}{}{}$".format(prefix, pattern, self.suffix), re.IGNORECASE)
                    len_match = len(list(filter(r.match, filt_candidates)))
                    if len_match > max_match:
                        max_match = len_match
                        self.prefix = prefix
                        self.pattern = pattern
        
    def has_suffix(self) -> bool:
        """
        Assume any sequence_region with a version number "_v#" is a likely candidate for a chromosome.
        Filter candidates accordingly.
        """
        with_suffix = list(filter(re.compile(r".*_v\d+$").match, self.chr_candidates))
        if any(with_suffix):
            self.chr_candidates = with_suffix
            self.suffix = re.search(r"_v\d+", with_suffix[0]).group()
            return True
        return False
    
    def _iter_prefix_pool(self):
        yield from common_prefixes(self.chr_candidates)
        # yield from [
        #     common_prefixes(self.chr_candidates),
        #     (self.chr_candidates[0].split("_")[0] + "_", 1),
        # ]
    
    def correct_prefix(self):
        try:
            r = re.compile(r"{}$".format(self.regex), re.IGNORECASE)
            first_match = list(filter(r.match, self.chr_candidates))[0]
            prefix = re.split('_|\.', first_match)[0]
            separator = first_match.split(prefix)[1][0]
            self.prefix = prefix + separator
        except:
            print(self.species)

def common_prefixes(li):
    prefixes = []
    for first_letter, prefix_batch in groupby(sorted(li), key=itemgetter(0)):
        threshold = len(li)
        prefix = []
        for chars in zip_longest(*map(str.lower, prefix_batch), fillvalue=''):
            char, count = Counter(chars).most_common(1)[0]
            if count == 1:
                break
            elif count < threshold:
                if prefix:
                    prefixes.append((''.join(prefix), threshold))
                    break
                threshold = count
            prefix.append(char)    
    return sorted(prefixes, key=itemgetter(1))


if __name__ == "__main__":
    gff = sys.argv[1]
    parser = ChromosomeParser(gff)
    parser.generate_regex()
    parser.correct_prefix()
    print(parser)
