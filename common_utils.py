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

   general functions
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"

import pandas as pd

def convert_pipe_separated_to_list(arg):
    """convert string of pipe separated terms to list
    
    do nothing if the argument is a list
    """
    if '|' in arg:
        return arg.split('|')
    else:
        return arg
    
def combine_dicts(dict1, dict2):
    """combine two dicts of key ->list by appending the lists associated with identical keys"""
    return dict((k, dict1.get(k,[]) + dict2.get(k,[])) for k in dict1.keys() | dict2.keys())


def combine_dictset(dict1, dict2):
    """combine two dicts of key -> set by merging the sets associated with identical keys"""
    return dict((k, dict1.get(k,set()) | dict2.get(k,set())) for k in dict1.keys() | dict2.keys())

def create_assoc_map(data_frame, column):
    """
    create an dict key -> set of values from a dataframe
    
    if a cell contains the pipe char, then a set is created
    """
    assoc = {}
    for key, value in data_frame.T.items():
        if pd.isnull(value[column]):
            continue
        try:
            if '|' in key: # extending the key
                keys = key.split('|')
            else:
                keys = [key]
        except TypeError:
            values = [value[column]]
        try:
            if '|' in value[column]: # extending the values
                values = value[column].split('|')
            else:
                values = [value[column]]
        except TypeError:
            values = [value[column]]
        # dispatching keys and values
        for k1 in keys:
            for v1 in values:
                if k1 not in assoc:
                    assoc[k1] = set([])
                assoc[k1].add(v1)
        values = set([])
        if isinstance(value[column], set):
            values.update(v[column])
        else:
            values.add(value[column])
    return assoc

def set_default(obj):
    """ default function for json encoder that transform set into list"""
    if isinstance(obj, set):
        return list(obj)
    raise TypeError

def get_stats_on_assocmap(assoc):
    diseases = set()
    count = 0
    for values in assoc.values():
        diseases.update(values)
        count += len(values)
    return("the matrix contains {0} rows, {1} columns and {2} associations".format(len(assoc), len(diseases), count))

def int2bin(number):
    """convert from positive integer to list of binary bits, msb at index 0"""
    if number:
        bits = []
        while number:
            number, remainder = divmod(number, 2)
            bits.insert(0, remainder)
        return bits
    else: return [0]

def bin2gray(bits):
    """convert from a vector of bits to gray code"""
    return bits[:1] + [i ^ ishift for i, ishift in zip(bits[:-1], bits[1:])]

def gray2bin(bits):
    """convert from Gray code to a vector of bits"""
    vbit = [bits[0]]
    for nextb in bits[1:]:
        vbit.append(vbit[-1] ^ nextb)
    return vbit

def bin2int(bits):
    "convert from binary bits, msb at index 0 to integer"""
    i = 0
    for bit in bits:
        i = i * 2 + bit
    return i
