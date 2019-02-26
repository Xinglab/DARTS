#!/usr/bin/env python

"""this is a simple script to concatenate space-separated
list into a comma-separated string
Zijun Zhang
10.4.2017
"""

import sys

x_list = sys.argv[1:]
print ','.join(x_list)
print '' 