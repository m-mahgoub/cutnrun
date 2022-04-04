#!/usr/bin/env python

import sys

outfile = open (sys.argv[2], 'w')
outfile.write ('sample_id,control' + '\n')

with open (sys.argv[1], 'r') as f:
    for line in f.readlines():
        line_as_list = line.rstrip().split(',')
        if line_as_list[0] != 'sample_id' and line_as_list[-1] != '' :
            sample_input_line = ','.join ([line_as_list[0], line_as_list[-1]]) + '\n'
            outfile.write(sample_input_line)
outfile.close()