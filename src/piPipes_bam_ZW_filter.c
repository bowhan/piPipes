/*
	# piPipes, a set of pipelines for PIWI-interacting RNA (piRNA) and transposon analysis
	# Copyright (C) 2014-2016  Bo Han, Wei Wang, Zhiping Weng, Phillip Zamore
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
*/

#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"

void usage (const char* name) {
	fprintf (stderr, "this script is part of piPipes package and used to filter the csem bam output based on the ZW:f tag.\n\nusage: %s %s %s[%.1f]\n", name, "input.bam", "threshold", 0.5);
}

uint8_t* check_bam_aux_get(const bam1_t* aln, const char* tag, char type) {
	uint8_t* p = bam_aux_get(aln, tag);
	if (p) {
		if (*p == type) return p;
		else {
			fprintf(stderr, "%s field of type '%c', expected '%c'\n", tag, *p, type);
			exit (1);
		}
	} else {
		fprintf(stderr, "can't find %s field\n", tag);
		exit (1);
	}
	return NULL;
}

int main(int argc, char** argv)
{
	double threshold = 0.5;
	if (argc < 2) {
		usage (argv[0]);
		return 1;
	}
	if (argc > 2) {
		threshold = atof (argv[2]);
	}
	samFile* in = sam_open(argv[1], "r");
	bam_hdr_t* header = sam_hdr_read(in);
	bam1_t* aln = bam_init1();
	uint8_t* p;
	int exit_code = 0;
	htsFile* out = hts_open("-", "wb");
	if (out == NULL) {
		fprintf(stderr, "Error opening stdout\n");
		exit_code = 1;
		return exit_code;
	}
	sam_hdr_write(out, header);
	while (sam_read1(in, header, aln) >= 0) {
		if (p = check_bam_aux_get(aln, "ZW", 'f')) {
			if (bam_aux2f(p) > threshold) {
				if (sam_write1(out, header, aln)<0) {
					fprintf(stderr, "error writing output\n");
					exit_code = 1;
					break;
				}
			}
		}
	}
	bam_destroy1(aln);
	bam_hdr_destroy(header);
	sam_close(out);
	sam_close(in);
	return exit_code;
}
