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

   this file contains the classe storing program parameters
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"
__all__ = ['MiRAIParam']

import os.path
from datetime import datetime

class MiRAIParam:
    """group program parameters"""

    DATADIR = os.path.normpath("data")#os.path.expandvars("data")
    TMPDIR = os.path.normpath(os.path.join(DATADIR, "generated"))
    OUTPUTDIR = os.path.normpath("output")

    TIMESTAMP = datetime.now().strftime("%Y-%m-%d:%H:%M:%S")

#     PROJECTION_METHODS = set(['nbi-train', 'tfidf-train', 'nbi-test', 'tfidf-test']) # possible values: 'nbi-train', 'tfidf-train', nbi-test', 'tfidf-test'
    PROJECTION_MATRIX1_NBI = True
    PROJECTION_MATRIX2_TFIDF = True
    PROJECTION_MATRIX3_NBI = True

    USE_PRODUCT_DATA = True
    USE_MIRBASE_ABSTRACT_DATA = True
    USE_PUBMED_ABSTRACT_DATA = True
    USE_DESCRIPTION_DATA = False
    USE_TARGET_DATA = True
    USE_NEIGHBOR_DATA = True
    USE_FAMILY_DATA = True
    USE_DISEASE_SIM_ON_DATA = True
    USE_DISEASE_SIM_ON_QUERY = False

    DISEASE_SIMILARITY_CUTOFF = [0.21, 0.7, 0.95]
    DISEASE_SIMILARITY_CUT1 = 0.21
    DISEASE_SIMILARITY_CUT2 = 0.7
    DISEASE_SIMILARITY_CUT3 = 0.95

    LSI_DIMENSION = 400

    FPR_CUTOFF = 0.85 # used to find dubious associations
    TPR_CUTOFF = 0.85 # used to predict new associations

    EXPAND_DISEASES = False
    # group miRNA that share similar mature forms
    # only used if USE_MATURE_FORM = False
    GROUP_WHEN_SIMILAR_MATURE = True
    REMOVE_OUTLIERS = False
    REMOVE_NOT_TESTED = False

    PERFORM_EVALUATION = True
    EVALUATION_FOR_ALL_TERMS = True
    EVALUATION_MIN_TERM_ASSOC = 20

    PERFORM_PREDICTION = False

    PROCESS_MIR2DISEASE = True
    MIR2DISEASE_FILES = {'hmdd': os.path.normpath(os.path.join(DATADIR, 'miRNA-assoc', 'hmdd_v2.txt')),
                         'mircancer': os.path.normpath(os.path.join(DATADIR, 'miRNA-assoc', 'miRCancerJune2015.txt')),
                         'mir2disease': os.path.normpath(os.path.join(DATADIR, 'miRNA-assoc', 'mir2disease.txt')),
                         'phenomir': os.path.normpath(os.path.join(DATADIR, 'miRNA-assoc', 'phenomir-2.0.txt'))}
    MIR2DISEASE_SOURCE = 'hmdd'

    MODEL_FILE = os.path.join(TMPDIR, 'model_raw')
    
    def __init__(self):
        self.DATADIR = MiRAIParam.DATADIR
        self.TMPDIR = MiRAIParam.TMPDIR
        self.OUTPUTDIR = MiRAIParam.OUTPUTDIR
        self.TIMESTAMP = MiRAIParam.TIMESTAMP
        self.PROJECTION_MATRIX1_NBI = MiRAIParam.PROJECTION_MATRIX1_NBI
        self.PROJECTION_MATRIX2_TFIDF = MiRAIParam.PROJECTION_MATRIX2_TFIDF
        self.PROJECTION_MATRIX3_NBI = MiRAIParam.PROJECTION_MATRIX3_NBI
        self.USE_PRODUCT_DATA = MiRAIParam.USE_PRODUCT_DATA
        self.USE_MIRBASE_ABSTRACT_DATA = MiRAIParam.USE_MIRBASE_ABSTRACT_DATA
        self.USE_PUBMED_ABSTRACT_DATA = MiRAIParam.USE_PUBMED_ABSTRACT_DATA
        self.USE_DESCRIPTION_DATA = MiRAIParam.USE_DESCRIPTION_DATA
        self.USE_TARGET_DATA = MiRAIParam.USE_TARGET_DATA
        self.USE_NEIGHBOR_DATA = MiRAIParam.USE_NEIGHBOR_DATA
        self.USE_FAMILY_DATA = MiRAIParam.USE_FAMILY_DATA
        self.USE_DISEASE_SIM_ON_DATA = MiRAIParam.USE_DISEASE_SIM_ON_DATA
        self.USE_DISEASE_SIM_ON_QUERY = MiRAIParam.USE_DISEASE_SIM_ON_QUERY
        self.DISEASE_SIMILARITY_CUTOFF = MiRAIParam.DISEASE_SIMILARITY_CUTOFF
        self.LSI_DIMENSION = MiRAIParam.LSI_DIMENSION
        self.FPR_CUTOFF = MiRAIParam.FPR_CUTOFF
        self.TPR_CUTOFF = MiRAIParam.TPR_CUTOFF
        self.EXPAND_DISEASES = MiRAIParam.EXPAND_DISEASES
        self.GROUP_WHEN_SIMILAR_MATURE = MiRAIParam.GROUP_WHEN_SIMILAR_MATURE
        self.REMOVE_OUTLIERS = MiRAIParam.REMOVE_OUTLIERS
        self.REMOVE_NOT_TESTED = MiRAIParam.REMOVE_NOT_TESTED
        self.PERFORM_EVALUATION = MiRAIParam.PERFORM_EVALUATION
        self.EVALUATION_FOR_ALL_TERMS = MiRAIParam.EVALUATION_FOR_ALL_TERMS
        self.EVALUATION_MIN_TERM_ASSOC = MiRAIParam.EVALUATION_MIN_TERM_ASSOC
        self.PERFORM_PREDICTION = MiRAIParam.PERFORM_PREDICTION
        self.PROCESS_MIR2DISEASE = MiRAIParam.PROCESS_MIR2DISEASE
        self.MIR2DISEASE_FILES = MiRAIParam.MIR2DISEASE_FILES
        self.MIR2DISEASE_SOURCE = MiRAIParam.MIR2DISEASE_SOURCE
        self.MODEL_FILE = MiRAIParam.MODEL_FILE
