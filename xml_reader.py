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

   this file contains functions related to the handling of MeSH XML
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"

import os.path
import xml.etree.ElementTree as ET

def read_pmid_to_abstract(filename):
    """read a Pubmed XML file and return a dictionary mapping PMID to abstract"""
    pmid2abstract = {}
    if os.path.getsize(filename) <= 1:
        return pmid2abstract
    tree = ET.parse(filename)
    root = tree.getroot()
    for desc in root.iter('PubmedArticle'):
        pmid = desc.find('MedlineCitation/PMID').text
        abstract = ""
        for term in desc.iter('AbstractText'):
            if term.text:
                abstract += ' ' + term.text
        pmid2abstract[pmid] = abstract
    return pmid2abstract
