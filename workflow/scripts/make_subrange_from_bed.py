#!/usr/bin/env python

from __future__ import print_function

import sys

max_range_width = int(sys.argv[1])

for line in sys.stdin:
        id, range_start, range_end = line.rstrip("\n").split("\t")
        range_end = int(range_end)
        assert range_end > 0, range_end

        subrange_start = 0
        while subrange_start < range_end:
                subrange_end = min(range_end, subrange_start + max_range_width)
                print("%s\t%s\t%s" % (id, subrange_start, subrange_end))
                subrange_start = subrange_end

