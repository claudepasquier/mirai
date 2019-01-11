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
__all__ = ['MirData']

import os.path
import pandas as pd
import common_utils
import gff_reader
from mirai_parameter import MiRAIParam

DATADIR = MiRAIParam.DATADIR

class MirData:
    """ group miRNA related data """

    _miRNAList = None
    _mirAliases = None
    _mature2miRNA = None
    _miRNA2mature = None
    _family2mir = None
    _mir2family = None
    _dead_miRNA = None
    _mature2target = None
    _mature2neighbor = None
    _miRNA2neighbor = None
    _ID2AC = None
    _AC2ID = None

    @staticmethod
    def read_aliases():
        MirData._mirAliases = {}
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'aliases.txt')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                line = line.split('\t')
                mir_ac = line[0].strip()
                mir_ids = line[1].strip(';').split(';')
                try:
                    MirData._mirAliases[mir_ac].update(set(mir_ids))
                except KeyError:
                    MirData._mirAliases[mir_ac] = set(mir_ids)
   
    @staticmethod
    def read_miRNAList():
        """Read the list of miRNAs"""
        MirData._miRNAList = set()
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miRNA.dat')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('ID'):
                    miRNA = line.split()[1]
                if line.startswith('DE   Homo sapiens'):
                    MirData._miRNAList.add(miRNA)
    
    @staticmethod
    def read_mature2mir():
        MirData._mature2miRNA = {}    
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miRNA.dat')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('ID'):
                    miRNA = line.split()[1]
                if line.startswith('FT                   /product="'):
                    product = line.split('"')[1]
                    try:
                        MirData._mature2miRNA[product].add(miRNA)
                    except KeyError:
                        MirData._mature2miRNA[product] = set([miRNA])

    @staticmethod
    def read_ID2AC():
        MirData._ID2AC = {}    
        MirData._AC2ID = {}    
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miRNA.dat')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('ID'):
                    miRNA_ID = line.split()[1]
                if line.startswith('AC'):
                    miRNA_AC = line.split()[1][:-1]
                    MirData._ID2AC[miRNA_ID] = miRNA_AC   
                    MirData._AC2ID[miRNA_AC] = miRNA_ID   
                if line.startswith('FT                   /accession="'):
                    mature_AC = line.split('"')[1]
                if line.startswith('FT                   /product="'):
                    product = line.split('"')[1]
                    MirData._ID2AC[product] = mature_AC   
                    MirData._AC2ID[mature_AC] = product   


    @staticmethod
    def read_mir2mature():
        MirData._miRNA2mature = {}    
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miRNA.dat')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('ID'):
                    miRNA = line.split()[1]
                if line.startswith('FT                   /product="'):
                    mature = line.split('"')[1]
                    try:
                        MirData._miRNA2mature[miRNA].add(mature)
                    except KeyError:
                        MirData._miRNA2mature[miRNA] = set([mature])
                        
    @staticmethod
    def read_miRNA_features(features_id):
        file_in = open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miRNA.dat')), 'r')
    #     file_in = open('miRNA.dat', 'r')
        assoc_map = {}
        for line in file_in:
            line = line.strip()
            if not line:
                continue
            if line.startswith('ID'):
                miRNA = line.split()[1]
                assoc_map[miRNA] = []
            if line[:2] in features_id:
                features = line[5:].strip('";.').lower().split()
                assoc_map[miRNA].extend(features)
        file_in.close()
        return assoc_map

    @staticmethod
    def read_mature2neighbors():
        filename = os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'hsa.gff3'))
        MirData._mature2neighbor = gff_reader.get_mature_to_neighbor(filename, 50000)    
    
    @staticmethod
    def read_miRNA2neighbors():
        filename = os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'hsa.gff3'))
        MirData._miRNA2neighbor = gff_reader.get_mirna_to_neighbor(filename, 50000)

    @staticmethod
    def read_mature2target():
        MirData._mature2target = {}
        file = os.path.normpath(os.path.join(DATADIR, 'mirtarbase.mbc.nctu.edu.tw', 'hsa_MTI.csv'))
        pd_mir2targets = pd.read_csv(file, delimiter='\t',  index_col='miRNA', usecols=[1,3], header=0)
        MirData._mature2target = common_utils.create_assoc_map(pd_mir2targets, 'Target Gene')


    @staticmethod
    def read_mirna_families():
        MirData._family2mir = {}    
        MirData._mir2family = {}    
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miFam.dat')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('ID'):
                    family = line.split()[1]
                if line.startswith('MI'):
                    miRNA = line.split()[2]
                    try:
                        MirData._family2mir[family].add(miRNA)
                    except KeyError:
                        MirData._family2mir[family] = set([miRNA])
                    if miRNA in MirData._mir2family:
                        raise Exception("problem with miRNA families")
                    MirData._mir2family[miRNA] = family
    
    @staticmethod
    def read_dead_miRNA():
        MirData._dead_miRNA = set()    
        with open(os.path.normpath(os.path.join(DATADIR, 'www.mirbase.org', 'miRNA.dead')), 'r') as infile:
            for line in infile:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('ID'):
                    miRNA = line.split()[1]
                    MirData._dead_miRNA.add(miRNA)

    @staticmethod
    def add_mir_type(arg):
        """add the type of the argument"""
        arg = arg.strip()
        if MirData.is_mature(arg):
            return 'M:' + arg
        if MirData.is_premirna(arg):
            return 'm:' + arg
        if MirData.is_family(arg):
            return 'f:' + arg
        if MirData.is_dead(arg):
            return 'd:' + arg
        # do a search in the aliases
        for mir in MirData.get_aliases(arg):
            if MirData.is_mature(mir):
                return 'M:' + mir
            if MirData.is_premirna(mir):
                return 'm:' + mir
            if MirData.is_family(mir):
                return 'f:' + mir
            if MirData.is_dead(mir):
                return 'd:' + mir
    
        # as some mature forms are sometimes written in lowercase 
        # do a search by using the standard name of the mature form
        if 'mir' in arg:
            arg = arg.replace('mir', 'miR')
            return MirData.add_mir_type(arg)
        return 'u:unknow'
    
    @staticmethod
    def get_aliases(mirId):
        """return a set of aliases of the parameter id"""
        if MirData._mirAliases == None:
            MirData.read_aliases()
        aliases = MirData._mirAliases.get(MirData.mirID2AC(mirId), set())
        return aliases - set(mirId)

    @staticmethod
    def get_all_spellings(mirId):
        """return a set of aliases and spelling forms of the parameter id"""
        spellings = set()
        spellings.add(mirId) 
        spellings.update(MirData.get_aliases(mirId))
#         if mirId == "hsa-mir-499a":
#             spellings.remove("hsa-mir-499b")
        if mirId == "hsa-mir-499b":
            spellings.remove("hsa-mir-499a")
        if mirId == "hsa-mir-29b-1":
            spellings.remove("hsa-mir-29b-2")
        if mirId == "hsa-mir-29b-2":
            spellings.remove("hsa-mir-29b-1")
        short_names = [t.replace('hsa-', '') for t in spellings]
        spellings.update(short_names)
        spellings.update(t.replace('mir-', 'miRNA-') for t in short_names)
        spellings.update(t.replace('mir-', 'miRNA ') for t in short_names)
        spellings.update(t.replace('mir-', 'miRNA') for t in short_names)
        spellings.update(t.replace('mir-', 'microRNA-') for t in short_names)
        spellings.update(t.replace('mir-', 'microRNA ') for t in short_names)
        spellings.update(t.replace('mir-', 'microRNA') for t in short_names)
        spellings.update(t.replace('miR-', 'miRNA-') for t in short_names)
        spellings.update(t.replace('miR-', 'miRNA ') for t in short_names)
        spellings.update(t.replace('miR-', 'miRNA') for t in short_names)
        spellings.update(t.replace('miR-', 'microRNA-') for t in short_names)
        spellings.update(t.replace('miR-', 'microRNA ') for t in short_names)
        spellings.update(t.replace('miR-', 'microRNA') for t in short_names)
        return spellings
    
    @staticmethod
    def mature2mir(mirId):
        """return a set of stem-loops corresponding to a mature form"""
        if MirData._mature2miRNA == None:
            MirData.read_mature2mir()
        return MirData._mature2miRNA.get(mirId, set())
    
    @staticmethod
    def mir2mature(mirId):
        """return a set of stem-loops corresponding to a mature form"""
        if MirData._miRNA2mature == None:
            MirData.read_mir2mature()
        return MirData._miRNA2mature.get(mirId, set())

    @staticmethod
    def mirID2AC(mirId):
        """return the accession number corresponding to a mirId"""
        if MirData._ID2AC == None:
            MirData.read_ID2AC()
        return MirData._ID2AC.get(mirId, "")

    @staticmethod
    def mirAC2ID(mirAC):
        """return the mirID corresponding to an accession number"""
        if MirData._AC2ID == None:
            MirData.read_ID2AC()
        return MirData._AC2ID.get(mirAC, "")
    
    @staticmethod
    def family2mir(famId):
        """return a set of stem-loops belonging to a family"""
        if MirData._family2mir == None:
            MirData.read_mirna_families()
        if ':' in famId:
            famId = famId[2:]
        return MirData._family2mir.get(famId, set())

    @staticmethod
    def mir2family(mirId):
        """return the family of a pre-mir"""
        if MirData._mir2family == None:
            MirData.read_mirna_families()
        if ':' in mirId:
            mirId = mirId[2:]
        return MirData._mir2family.get(mirId, None)
    
    @staticmethod
    def is_mature(mirId):
        """return True if mirId identifies a mature form"""
        if mirId.startswith('M:'):
            return True
        if ':' in mirId:
            return False
        if MirData._mature2miRNA == None:
            MirData.read_mature2mir()
        return mirId in MirData._mature2miRNA
    
    @staticmethod
    def is_premirna(mirId):
        """return True if mirId identifies a trem-loop"""
        if mirId.startswith('m:'):
            return True
        if ':' in mirId:
            return False
        if MirData._miRNA2mature == None:
            MirData.read_mir2mature()
        return mirId in MirData._miRNA2mature
            
    @staticmethod
    def is_family(mirId):
        """return True if mirId identifies a family"""
        if mirId.startswith('f:'):
            return True
        if ':' in mirId:
            return False
        if MirData._family2mir == None:
            MirData.read_mirna_families()
        return mirId in MirData._family2mir

    @staticmethod
    def is_dead(mirId):
        """return True if mirId identifies a family"""
        if mirId.startswith('d:'):
            return True
        if ':' in mirId:
            return False
        if MirData._dead_miRNA == None:
            MirData.read_dead_miRNA()
        return mirId in MirData._dead_miRNA

    @staticmethod
    def get_premir_to_mature():
        """return the mapping from premir to mature"""
        if MirData._miRNA2mature == None:
            MirData.read_mir2mature()
        return MirData._miRNA2mature

    @staticmethod
    def get_mature_to_target():
        """return the mapping from mature to target"""
        if MirData._mature2target == None:
            MirData.read_mature2target()
        return MirData._mature2target

    @staticmethod
    def get_mir_to_neighbor():
        """return the mapping from miRNA to theirs neighbors"""
        if MirData._miRNA2neighbor == None:
            MirData.read_miRNA2neighbors()
        return MirData._miRNA2neighbor

    @staticmethod
    def get_mir_to_family():
        """return the mapping from miRNA to their family"""
        if MirData._mir2family == None:
            MirData.read_mirna_families()()
        return MirData._mir2family
    
    @staticmethod
    def get_family2mir():
        """return the mapping from family to member miRNA"""
        if MirData._family2mir == None:
            MirData.read_mirna_families()()
        return MirData._family2mir

    @staticmethod
    def get_mature_to_neighbor():
        """return the mapping from mature miRNA to theirs neighbors"""
        if MirData._mature2neighbor == None:
            MirData.read_mature2neighbors()
        return MirData._mature2neighbor
    
    @staticmethod
    def get_miRNAs():
        """get the list of all pre miRNAs"""
        if MirData._miRNAList == None:
            MirData.read_miRNAList()
        return MirData._miRNAList
        
    
    @staticmethod
    def to_mature(arg):
        """return the set of mature ids corresponding to the argument"""
        if ':' not in arg:
            arg = MirData.add_mir_type(arg)
        mirType = arg[0]
        mirId = arg[2:]
        if mirType == 'M':
            return mirId
        if mirType == 'm':
            return MirData.mir2mature(mirId)
        if mirType == 'f':
            matures = set()
            for mir in MirData.family2mir(mirId):
                matures.update(MirData.mir2mature(mir))
            return matures
        if mirType == 'd':
            return None
        if mirType == 'u':
            return None
        
    @staticmethod
    def to_premir(arg):
        """return the set of step-loop ids corresponding to the argument"""
        if ':' not in arg:
            arg = MirData.add_mir_type(arg)
        mirType = arg[0]
        mirId = arg[2:]
        if mirType == 'M':
            return MirData.mature2mir(mirId)
        if mirType == 'm':
            return mirId
        if mirType == 'f':
            return MirData.family2mir(mirId)
        if mirType == 'd':
            return None
        if mirType == 'u':
            return None

if __name__ == '__main__':
    pass
