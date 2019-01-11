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

   this file contains classes and methods related to miRNA information
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"
__all__ = ['MeSHData']

import logging
import os.path
import re
import xml.etree.ElementTree as ET
import json
import sys
from mirai_parameter import MiRAIParam



DATADIR = MiRAIParam.DATADIR
XMLFILE = os.path.normpath(os.path.join(DATADIR, 'www.nlm.nih.gov', 'mesh', 'desc2015.xml'))
XMLSUPFILE = os.path.normpath(os.path.join(DATADIR, 'www.nlm.nih.gov', 'mesh', 'supp2015.xml'))
CURRENTDIR = os.path.dirname(sys.argv[0])
TMPDIR = os.path.normpath(os.path.join(DATADIR, "generated"))

JSONFILE = os.path.join(TMPDIR, "meshData.json")

class MeSHData:
    """group MeSH related data"""

    term2id = {}
    id2term = {}
    tree2id = {}
    id2tree = {}

    alias2mesh = {}
    mesh2alias = {}

    @staticmethod
    def reverse_dicts():
        """reverse dictionaries"""
        MeSHData.id2term = {value: key for key, value in MeSHData.term2id.items()}
        for key, value in MeSHData.tree2id.items():
            try:
                MeSHData.id2tree[value].add(key)
            except KeyError:
                MeSHData.id2tree[value] = set([key])
        for key, value in MeSHData.alias2mesh.items():
            try:
                MeSHData.mesh2alias[value].add(key)
            except KeyError:
                MeSHData.mesh2alias[value] = set([key])

    @staticmethod
    def read_data():
        """read Mesh data from json or XML file"""
        if os.path.exists(JSONFILE):
            logging.info("## reading MeSH data from %s", JSONFILE)
            with open(JSONFILE, 'r') as infile:
                json_data = json.load(infile)
                MeSHData.term2id = json_data[0]
                MeSHData.tree2id = json_data[1]
                MeSHData.alias2mesh = json_data[2]
                MeSHData.reverse_dicts()
        else:
            logging.info("## reading MeSH data from " + XMLFILE)
            tree = ET.parse(XMLFILE)
            root = tree.getroot()
            for desc in root.iter('DescriptorRecord'):
                name = desc.find('DescriptorName/String').text
                mesh_id = desc.find('DescriptorUI').text
#                 MeSHData.id2term[mesh_id] = name
                MeSHData.term2id[name] = mesh_id
                for tree in desc.iter('TreeNumber'):
                    tree_id = tree.text
                    assert tree_id not in MeSHData.tree2id
                    MeSHData.tree2id[tree_id] = mesh_id
#                     try:
#                         MeSHData.id2tree[mesh_id].add(tree_id)
#                     except KeyError:
#                         MeSHData.id2tree[mesh_id] = set([tree_id])
                MeSHData.alias2mesh[name.lower()] = name
                for term in desc.iter('Term'):
                    alias = term.findtext('String')
                    MeSHData.alias2mesh[alias.lower()] = name
            MeSHData.read_supplement_data()
            MeSHData.reverse_dicts()
            with open(JSONFILE, 'w') as outfile:
                json.dump([MeSHData.term2id, MeSHData.tree2id, MeSHData.alias2mesh], outfile, indent=2)

    @staticmethod
    def read_supplement_data():
        """read supplement Mesh data from XML file"""
        tree = ET.parse(XMLSUPFILE)
        root = tree.getroot()
        for desc in root.iter('SupplementalRecord'):
            mesh_mapping = []
            for term in desc.iter('DescriptorName'):
                mesh_mapping.append(term.findtext('String'))
            mesh_mapping = "|".join(mesh_mapping)
            aliases = []
            aliases.append(desc.find('SupplementalRecordName/String').text)
            for term in desc.iter('ConceptName'):
                aliases.append(term.findtext('String'))
            for term in desc.iter('Term'):
                aliases.append(term.findtext('String'))

            for alias in aliases:
                MeSHData.alias2mesh[alias.lower()] = mesh_mapping

    @staticmethod
    def get_id(term=None, tree=None):
        """retrieve the MeSH id from the term or a mesh tree"""
        if term:
            return MeSHData.term2id[term]
        else:
            return MeSHData.tree2id[tree]

    @staticmethod
    def get_term(mesh_id=None, tree=None):
        """retrieve the MeSH term from the term ID or a mesh tree"""
        if mesh_id:
            return MeSHData.id2term[mesh_id]
        else:
            return MeSHData.id2term[MeSHData.tree2id[tree]]

    @staticmethod
    def get_tree(mesh_id=None, term=None):
        """retrieve the MeSH trees from a term ID or a mesh term"""
        if mesh_id:
            return MeSHData.id2tree[mesh_id]
        else:
            return MeSHData.id2tree[MeSHData.term2id[term]]

    @staticmethod
    def convert2mesh(arg):
        """return the mesh term corresponding to the arg"""
        # remove expression inside brackets
        modified_arg = re.sub(r"\(.*?\)", "", arg)
        # strip line
        modified_arg = modified_arg.strip()
        # capitalize each words
        modified_arg = modified_arg.lower()
        if modified_arg in MeSHData.alias2mesh:
            return MeSHData.alias2mesh.get(modified_arg)
        return MeSHData.alias2mesh.get(arg, "not_mapped:" + arg)

    @staticmethod
    def convert2meshid(arg):
        """return the mesh id corresponding to the arg"""
        if arg in MeSHData.term2id:
            return MeSHData.term2id.get(arg)
        return arg

    @staticmethod
    def get_childs(value):
        """retrieve the child terms of the parameter in mesh"""
        if value not in MeSHData.term2id:
            return set()
        mesh_ids = MeSHData.get_tree(term=value)

#         if isinstance(mesh_ids, str):
#             mesh_ids = set([mesh_ids])
        # retrieve all descendant of the term
        childs = set()
        for tree, meshid in MeSHData.tree2id.items():
            for mesh_id in mesh_ids:
                if tree != mesh_id and tree.startswith(mesh_id):
                    childs.add(MeSHData.get_term(mesh_id=meshid))
        return childs
