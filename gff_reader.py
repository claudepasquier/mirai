#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""Copyright 2015-2019 UNS-CNRS

   This file is part of MiRAI

   MiRAI is free software: you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   MiRAI is distributed in the hope that it will be useful, but WITHOUT
   ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
   License for more details.

   You should have received a copy of the GNU General Public License
   along with MiRAI.  If not, see <http://www.gnu.org/licenses/>


   Author: Claude Pasquier (I3S Laboratory, UCA)
   Contact: claude.pasquier@unice.fr 

   Created on 26 mai 2015

   utilities for reading GFF files
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"

def get_neighbors(filename, distance, mirna_type):
    """Construct a dict mapping between miRNAs and their neighbors

    filename is the location of the GFF file
    distance is the minimum distance between a miRNA and its neighbor
    mirna_type is the type of miRNA :
        'miRNA_primary_transcript' for stem-loop
        'miRNA' for mature form
    """
    mirna_to_neighbor = {}
    memo_loc = []
    with open(filename) as infile:
        for line in infile:
            if not line:
                continue
            if line.startswith('#'):
                continue
            line = line.split('\t')
            if line[2] != mirna_type:
                continue
            chromosome = line[0]
            start_pos = line[3]
            end_pos = line[4]
            attrib = line[8].strip().split(';')
            for att in attrib:
                key_value = att.split('=')
                if key_value[0] == 'Name':
                    mir_id = key_value[1]
                    break
            if memo_loc:
                ctr = 0
                while ctr < len(memo_loc):
                    if memo_loc[ctr][1] != chromosome:
                        memo_loc.pop(ctr)
                    elif memo_loc[ctr][3] < int(start_pos) - distance:
                        memo_loc.pop(ctr)
                    else:
                        if mir_id not in mirna_to_neighbor:
                            mirna_to_neighbor[mir_id] = set([])
                        mirna_to_neighbor[mir_id].add(memo_loc[ctr][0])
                        if memo_loc[ctr][0] not in mirna_to_neighbor:
                            mirna_to_neighbor[memo_loc[ctr][0]] = set()
                        mirna_to_neighbor[memo_loc[ctr][0]].add(mir_id)
                        ctr += 1
            memo_loc.append((mir_id, chromosome, int(start_pos), int(end_pos)))
    return mirna_to_neighbor

def get_mirna_to_neighbor(filename, distance):
    """Construct a dict mapping miRNAs to neighbor miRNAs

    filename is the location of the GFF file
    distance is the minimum distance between a miRNA and its neighbor
    """
    return get_neighbors(filename, distance, 'miRNA_primary_transcript')

def get_mature_to_neighbor(filename, distance):
    """Construct a dict mapping mature miRNAs to neighbor miRNAs

    filename is the location of the GFF file
    distance is the minimum distance between a miRNA and its neighbor
    """
    return get_neighbors(filename, distance, 'miRNA')

if __name__ == '__main__':
    pass
