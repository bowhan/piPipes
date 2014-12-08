#! /usr/bin/env python

# piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
# Copyright (C) 2014  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

'''
this string extract sequence for a specific organism from fasta file
it requires the 3-letter abbreviation of the organism to be repsent in the first 3 letters
of the header, like
>cel-let-7-5p MIMAT0000001 Caenorhabditis elegans let-7-5p
TGAGGTAGTAGGTTGTATAGTT
'''

def sep_sequence_from_fa(file, str):
    match = False
    for line in open(file):
        if line[0] == '>':
            if line[1:4] == str:
                print line,
                match = True
            else:
                match = False
        elif match:
            print line.replace('U', 'T'),

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 3:
        print "usage: %s input.fa dme" % sys.argv[0]
        sys.exit(1)
    sep_sequence_from_fa(sys.argv[1], sys.argv[2])