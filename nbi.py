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

   this file contains functions related to NBI algorithm
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"

import numpy as np
import gensim.matutils

def create_matrix(corpus, nb_columns, val2index):
    """create a numpy matrix from a gensim corpus"""
#     print(len(corpus), nb_columns)
    M = np.zeros((len(corpus), nb_columns))
    for num_row in range(len(corpus)):
        for t in corpus[num_row]:
            M[num_row, val2index[t[0]]] = t[1]
    return M

def do_nbi_old(corpus, nb_columns):
    """perform NBI on the corpus"""
    # remove empty rows in corpus
    ignored_rows = set()
    new_corpus = []
    for num_row in range(len(corpus)):
        if not corpus[num_row]:
            ignored_rows.add(num_row)
            continue
        else:
            new_corpus.append(corpus[num_row])
    corpus = new_corpus
    
    val2index = {}
    index2val = {}
    index = 0
    for l in corpus:
        for v in l:
            if v[0] not in val2index:
                val2index[v[0]] = index
                index2val[index] = v[0]
                index += 1
    nb_columns = index
    M = create_matrix(corpus, nb_columns, val2index)
    M = M.astype(float, order='C')
    F1 = nbi_iterate(M, 0.95)

    new_corpus = []
    for row in range(F1.shape[0]):
        corpus_row = []
        for column in range(F1.shape[1]):
            val = F1[row, column]
            if val > 0:
                corpus_row.append((index2val[column],val))
        new_corpus.append(corpus_row)
    nb_row_in_corpus = 0
    num_row = 0
    corpus = []
    try:
        while True:
            if nb_row_in_corpus in ignored_rows:
                corpus.append([])
                nb_row_in_corpus += 1
                continue
            corpus.append(new_corpus[num_row])
            nb_row_in_corpus += 1
            num_row +=1
    except IndexError:
        pass
    return corpus

def do_nbi(corpus, nb_columns):
    """perform NBI on the corpus by using gensim conversion functions"""
    max_value = 0
    for l in corpus:
        for v in l:
            if v[0] > max_value:
                max_value = v[0]
    nb_columns = max_value+1
    M = gensim.matutils.corpus2dense(corpus, nb_columns)

    sum_rows = np.zeros(M.shape[0])
    sum_columns = np.zeros(M.shape[1])
    for row in range(M.shape[0]):
        for column in range(M.shape[1]):
            val = M[row, column]
            if val > 0:
                sum_rows[row] += val
                sum_columns[column] += val
                
    R = np.diag(sum_rows)
    H = np.diag(sum_columns)
    Wmm = np.dot(np.dot(M,np.linalg.inv(H)).T,np.dot(np.linalg.inv(R),M))
    F1 = np.dot(M, Wmm)
    
    new_corpus = gensim.matutils.Dense2Corpus(F1)
    return new_corpus


def nbi_iterate(M, alpha):
    """core of NBI calculation"""
    sum_rows = np.zeros(M.shape[0])
    sum_columns = np.zeros(M.shape[1])
    for row in range(M.shape[0]):
        for column in range(M.shape[1]):
            val = M[row, column]
            if val > 0:
                sum_rows[row] += val
                sum_columns[column] += val
    R = np.diag(sum_rows)
    H = np.diag(sum_columns)
    Wmm = np.dot(np.dot(M,np.linalg.inv(H)).T,np.dot(np.linalg.inv(R),M))
#     Wnn = np.dot(np.dot(np.linalg.inv(R),M),np.dot(M,np.linalg.inv(H)).T)

    F1 = alpha * np.dot(M, Wmm) + (1 - alpha) * M
#     F1a = np.dot(M.T,Wnn).T
    return F1

if __name__ == '__main__':
    F = np.array([[1,0,0,1,0,0],
                   [0,1,0,0,0,0],
                   [1,0,0,0,1,0],
                   [0,0,1,0,0,1],
                   [0,1,0,0,1,0]])
    F = np.array([[1,0,0],
                 [1,0,1],
                 [0,1,1]])
      
    F = np.array([[1,0,0,0],
                  [1,1,0,0],
                  [0,1,1,0],
                  [0,1,0,1],
                  [0,0,1,1]])


#     M = np.random.rand(500,1000)
    alph = 1
    for ctr in range(1):
        F = nbi_iterate(F, alph)
    #print(F)
#     print("finished")
#     F2 = nbi_iterate(F1)
    
#     M = np.array([[1,0,1],[1,1,0],[0,0,1]])


#     sim = np.array([[1,0.2,0],[0.2,1,0.8],[0,0.8,1]])
#     M = np.dot(F0, sim)
#     print(M)
#     m = M.max(axis=1)
#     m = m[:, np.newaxis]
#     print(M / m)
#     dfd
#     
#     
#     F0 = np.array([[1,0,1],[1,1,0],[0,0,1]])
# 
#     R=np.diag([2,2,1])
#     H=np.diag([2,1,2])
# 
#     Wmm = np.dot(np.dot(F0,np.linalg.inv(H)).T,np.dot(np.linalg.inv(R),F0))
#     Wnn = np.dot(np.dot(np.linalg.inv(R),F0),np.dot(F0,np.linalg.inv(H)).T)
# 
#     F1 = np.dot(F0, Wmm)
#     print(F0)
#     print(F1)
#     
#     # let say tha columns 0 an 1 are similars
#     sim_columns = np.array([[1,0.8,0],[0.8,1,0],[0,0,0]])
#     print(sim_columns)
#     print(F0.dot(np.linalg.inv(sim_columns)))
