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

   this program retrieve abstacts corresponding to the PMID used in mirBase
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"

import urllib.request
import os.path
import shutil
from mirdata import MirData
from mirai_parameter import MiRAIParam

def main():
    """main function"""
    feat = MirData.read_miRNA_features(set(['RX']))
    pmids = set()
    for values in feat.values():
        for value in values:
            try:
                pmids.add(int(value))
            except ValueError:
                pass

    query = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&rettype=abstract&retmode=xml&id=" + ','.join([str(x) for x in pmids])
    fileout, headers = urllib.request.urlretrieve(query)
    shutil.move(fileout, os.path.join(MiRAIParam.DATADIR, 'generated', 'pubmed', 'abstracts_of_mirbase.xml'))
#     print("the result can be found in the file: ", fileout)

if __name__ == '__main__':
    main()
