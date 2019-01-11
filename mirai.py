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

   entry point of the MiRAI software
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"

import sys
import logging
from gensim import corpora, models, similarities, matutils
import random
import nbi
from meshdata import MeSHData
from diseasedata import DiseaseData
import pandas as pd
import json
import yaml
import common_utils
import numpy as np
from sklearn import metrics
import scipy.spatial
from mirdata import MirData
import xml_reader
import os.path
import matplotlib
from mirai_parameter import MiRAIParam

matplotlib.use('SVG')
import matplotlib.pyplot as plt


def align_dicts(dict1, dict2):
    """align the two dict in order that they have the same keys"""
    for key in dict1:
        if key not in dict2:
            dict2[key] = []

def create_word_map(assoc, param):
    """create word map storing word frequencies"""
    word_map = {}
    if param.USE_DESCRIPTION_DATA:
        word_map = common_utils.combine_dicts(word_map, add_features(assoc, set(['CC'])))
    if param.USE_MIRBASE_ABSTRACT_DATA:
        word_map = common_utils.combine_dicts(word_map, add_abstract(assoc, param))
    if param.USE_PUBMED_ABSTRACT_DATA:
        word_map = common_utils.combine_dicts(word_map, add_abstract2(assoc, param))
    align_dicts(assoc, word_map)
    return word_map

def create_target_map(assoc, param):
    """create target map storing targets"""
    extension = {}
    if param.USE_TARGET_DATA:
        extension = common_utils.combine_dictset(extension, add_targets(assoc))
    align_dicts(assoc, extension)
    return extension

def create_neighbor_map(assoc, param):
    """create neighbor map storing neighbors"""
    extension = {}
    if param.USE_NEIGHBOR_DATA:
        extension = common_utils.combine_dictset(extension, add_neighbors(assoc))
    align_dicts(assoc, extension)
    return extension

def create_family_map(assoc, param):
    """create family map storing families"""
    extension = {}
    if param.USE_FAMILY_DATA:
        extension = common_utils.combine_dictset(extension, add_family(assoc))
    align_dicts(assoc, extension)
    return extension

def create_assoc_map_ext2_to_delete(assoc, param):
    """create assoc map storing associated data"""
    extension = {}
    if param.USE_PRODUCT_DATA:
        extension = common_utils.combine_dictset(extension, add_products(assoc))
    if param.USE_TARGET_DATA:
        extension = common_utils.combine_dictset(extension, add_targets(assoc))
    if param.USE_NEIGHBOR_DATA:
        extension = common_utils.combine_dictset(extension, add_neighbors(assoc))
    if param.USE_FAMILY_DATA:
        extension = common_utils.combine_dictset(extension, add_family(assoc))

    # align the two dict in order that they have the same keys
    for key in assoc:
        if key not in extension:
            extension[key] = set()
    return extension

def extend_disease(set_of_diseases):
    """extend diseases in a set with their ancestors"""

    to_add = set([])
    for value in set_of_diseases:
        if value not in MeSHData.term2id:
            continue
        mesh_ids = MeSHData.get_tree(term=value)
        if isinstance(mesh_ids, str):
            mesh_ids = set([mesh_ids])
        for mesh_id in mesh_ids:
            id_split = mesh_id.split('.')
            for indx in range(min(len(id_split)-2, 5)):
                ancestor_term = MeSHData.get_term(tree='.'.join(id_split[:-(indx+1)]))
                to_add.update(ancestor_term)
    set_of_diseases.update(to_add)

def extend_similarities(set_of_diseases, ext_simils):
    """extend diseases in a set with similar diseases stored in a dict"""
    to_add = set([])
    for value in set_of_diseases:
        if value not in MeSHData.term2id:
            continue
        term_id = MeSHData.get_id(term=value)
        to_add.update(ext_simils.get(term_id, set([])))
    set_of_diseases.update(to_add)

def add_products(assoc):
    """add products associated with miRNA in assoc dict"""
    extension = {}
    mir_to_product = MirData.get_premir_to_mature()
    for key in assoc:
        if MirData.is_premirna(key):
            if key in extension:
                extension[key].update(mir_to_product.get(key[2:], set()))
            else:
                extension[key] = mir_to_product.get(key[2:], set())
        elif MirData.is_mature(key):
            for mir in MirData.mature2mir(key[2:]):
                if key in extension:
                    extension[key].update(mir_to_product.get(mir, set()))
                else:
                    extension[key] = mir_to_product.get(mir, set())
        else:
            raise
    return extension

def add_targets(assoc):
    """add targets associated with miRNA in assoc dict"""
    extension = {}
    mir_to_targets = MirData.get_mature_to_target()
    for key in assoc:
        if MirData.is_premirna(key):
            for mir in MirData.mir2mature(key[2:]):
                if key in extension:
                    extension[key].update(mir_to_targets.get(mir, set()))
                else:
                    extension[key] = mir_to_targets.get(mir, set())
        elif MirData.is_mature(key):
            if key in extension:
                extension[key].update(mir_to_targets.get(key[2:], set()))
            else:
                extension[key] = mir_to_targets.get(key[2:], set())
        else:
            raise
    return extension

def add_neighbors(assoc):
    """add neighbors of miRNA in assoc dict"""
    extension = {}
    mir_to_neighbors = MirData.get_mir_to_neighbor()
    for key in assoc:
        if MirData.is_premirna(key):
            neighbors = mir_to_neighbors.get(key[2:], set())
            for neighbor in neighbors:
                if key in extension:
                    extension[key].add("cluster:" + neighbor)
                else:
                    extension[key] = set(["cluster:" + neighbor])
        elif MirData.is_mature(key):
            for mir in MirData.mature2mir(key[2:]):
                neighbors = mir_to_neighbors.get(mir, set())
                for neighbor in neighbors:
                    if key in extension:
                        extension[key].add("cluster:" + neighbor)
                    else:
                        extension[key] = set(["cluster:" + neighbor])
        else:
            raise Exception("miRNA not identified")
    return extension

def add_family(assoc):
    """add family of miRNA in assoc dict"""
    extension = {}
    mir_to_family = MirData.get_mir_to_family()
    for key in assoc:
        if MirData.is_premirna(key):
            family = mir_to_family.get(key[2:], '')
            if key in extension:
                extension[key].add("family:" + family)
            else:
                extension[key] = set(["family:" + family])
        elif MirData.is_mature(key):
            for mir in MirData.mature2mir(key[2:]):
                family = mir_to_family.get(mir, '')
                if key in extension:
                    extension[key].add("family:" + family)
                else:
                    extension[key] = set(["family:" + family])
        else:
            raise Exception("miRNA not identified")
    return extension

def add_features(assoc, features):
    """add features of miRNA in assoc dict"""
    extension = {}
    mir_to_features = MirData.read_miRNA_features(features)
    for key in assoc:
        if MirData.is_premirna(key):
            if key in extension:
                extension[key].extend(mir_to_features.get(key[2:], []))
            else:
                extension[key] = mir_to_features.get(key[2:], [])
        elif MirData.is_mature(key):
            for mir in MirData.mature2mir(key[2:]):
                if key in extension:
                    extension[key].extend(mir_to_features.get(mir, []))
                else:
                    extension[key] = mir_to_features.get(mir, [])
        else:
            raise
    return extension

def add_abstract(assoc, param):
    """add terms from abstracts of papers associated to miRNA in assoc dict"""
    mir_to_pmid = MirData.read_miRNA_features(set(['RX']))

    filename = os.path.normpath(os.path.join(param.DATADIR, 'generated', 'pubmed', 'abstracts_of_mirbase.xml'))
    pmid2abstract = xml_reader.read_pmid_to_abstract(filename)

    extension = {}
    for key in assoc:
        if MirData.is_premirna(key):
            pmids = mir_to_pmid.get(key[2:], [])
            for pmid in pmids:
                try:
                    int(pmid)
                except ValueError:
                    continue
                if key in extension:
                    extension[key].extend(pmid2abstract[pmid].split())
                else:
                    extension[key] = pmid2abstract[pmid].split()
        elif MirData.is_mature(key):
            for mir in MirData.mature2mir(key[2:]):
                pmids = mir_to_pmid.get(mir, [])
                for pmid in pmids:
                    try:
                        int(pmid)
                    except ValueError:
                        continue
                    if key in extension:
                        extension[key].extend(pmid2abstract[pmid].split())
                    else:
                        extension[key] = pmid2abstract[pmid].split()
        else:
            raise
    return extension

def add_abstract2(assoc, param):
    """add terms from all abstracts of papers associated to miRNA in assoc dict"""
    extension = {}
    for key in assoc:
        if MirData.is_premirna(key):
            filename = os.path.normpath(os.path.join(param.DATADIR, 'generated', 'pubmed', 'mirPapers', key[2:]+'.xml'))
            if not os.path.isfile(filename):
                continue
            #print(filename)
            pmid2abstract = xml_reader.read_pmid_to_abstract(filename)
            for _, abstract in pmid2abstract.items():
                if key in extension:
                    extension[key].extend(abstract.split())
                else:
                    extension[key] = abstract.split()
        elif MirData.is_mature(key):
            for mir in MirData.mature2mir(key[2:]):
                filename = os.path.normpath(os.path.join(param.DATADIR, 'generated', 'pubmed', 'mirPapers', mir[2:]+'.xml'))
                if not os.path.isfile(filename):
                    continue
                pmid2abstract = xml_reader.read_pmid_to_abstract(filename)
                for _, abstract in pmid2abstract.items():
                    if key in extension:
                        extension[key].extend(abstract.split())
                    else:
                        extension[key] = abstract.split()
        else:
            raise
    return extension

def lsi_search(search_term, dictionary2, index2, lsi2):
    """perform LSI search"""
    vec_bow = dictionary2.doc2bow(search_term)
    vec_lsi = lsi2[vec_bow] # convert the query to LSI space
    sims = index2[vec_lsi]
    return sorted(enumerate(sims), key=lambda item: -item[1])

def remove_rare_or_frequent(assoc):
    """remove miRNAs associated with many or few diseases"""
    disease_count = {}
    all_keys = list(assoc.keys())
    for value in assoc.values():
        for val in value:
            disease_count[val] = disease_count.get(val, 0) + 1
    disease_to_remove = set([])
    for key, value in disease_count.items():
        if value < len(all_keys)/100:
            disease_to_remove.add(key)
        if value > len(all_keys) * 0.99:
            disease_to_remove.add(key)
    for value in assoc.values():
        value -= disease_to_remove

    #remove miRNA associated with few diseasese
    for key in all_keys:
        if len(assoc[key]) < 1:
            del assoc[key]

def build_corpus(assoc, keys, global_dictionary, ext_simil, query_values=None, dummyline=None, nb_assoc_to_keep=0, param=None):
    """build the corpus used by LSI"""
    texts = []
    initial_values = set()
    ctr = 0
    for key in keys:
        values = assoc.get(key).copy()
        if query_values:
            if param.USE_DISEASE_SIM_ON_QUERY:
                queries = set()
                for term in query_values:
                    queries.add(term)
                    queries.update(MeSHData.get_childs(term))
                if values.intersection(queries):
                    ctr += 1
                if ctr > nb_assoc_to_keep:
                    values -= queries
            else:
                if values.intersection(query_values):
                    ctr += 1
                if ctr > nb_assoc_to_keep:
                    values -= query_values
        initial_values.update(values)
        if param.EXPAND_DISEASES:
            extend_disease(values)
        if param.USE_DISEASE_SIM_ON_DATA:
            extend_similarities(values, ext_simil)

        texts.append(sorted(values))
        if 'up' in values:
            texts[-1].append(' '.join(['up'+str(x) for x in range(10)]))
        if 'down' in values:
            texts[-1].append(" ".join(['down'+str(x) for x in range(10)]))
    if dummyline:
        texts.append(dummyline) # dummy line
    corpus = [global_dictionary.doc2bow(text) for text in texts]
    new_corpus = []
    for row in corpus:
        line = []
        for term in row:
            if global_dictionary.get(term[0]) in initial_values:
                line.append((term[0], term[1]))
            else:
                line.append((term[0], int(term[1])/2))
        new_corpus.append(line)
    return new_corpus

def build_corpus_ext(assoc, keys, global_dictionary, dummyline=None):
    """build the ext corpus used by LSI"""
    texts = []
    for key in keys:
        values = assoc[key].copy()
        texts.append(sorted(values))
    if dummyline:
        texts.append(dummyline) # dummy line
    corpus = [global_dictionary.doc2bow(text) for text in texts]
    return corpus

def build_corpus_ext2(assoc, keys, global_dictionary, dummyline=None):
    """build the ext2 corpus used by LSI"""
    texts = []
    for key in keys:
        values = assoc[key].copy()
        texts.append(sorted(values))
    if dummyline:
        texts.append(dummyline) # dummy line
    corpus = [global_dictionary.doc2bow(text) for text in texts]
    new_corpus = []
    for row in corpus:
        line = []
        for term in row:
            line.append((term[0], int(term[1])/4))
        new_corpus.append(line)
    return new_corpus

def extend_families(assoc):
    """replace families by the miRNAs in them"""
    to_delete = set()
    to_add = {}
    for key, value in assoc.items():
        if MirData.is_family(key):
            to_delete.add(key)
            mirs = MirData.family2mir(key)
            for mir in mirs:
                typed_mir = 'm:' + mir
                try:
                    to_add[typed_mir] = to_add[typed_mir] | value
                except KeyError:
                    to_add[typed_mir] = value
    for item in to_delete:
        del assoc[item]
    for key, value in to_add.items():
        try:
            assoc[key] = assoc[key] | value
        except KeyError:
            assoc[key] = value

def map_to_mature(assocs):
    """ map keys of the dict to mature forms"""
    for assoc in assocs:
        to_delete = set()
        to_add = {}
        for key, value in assoc.items():
            if MirData.is_premirna(key):
                to_delete.add(key)
                key = key[2:]
                mirs = MirData.mir2mature(key)
                for mir in mirs:
                    typed_mir = 'M:' + mir
                    try:
                        to_add[typed_mir] = to_add[typed_mir] | value
                    except TypeError:
                        to_add[typed_mir] = to_add[typed_mir] + value
                    except KeyError:
                        to_add[typed_mir] = value
        for item in to_delete:
            del assoc[item]
        for key, value in to_add.items():
            try:
                assoc[key] = assoc[key] | value
            except TypeError:
                assoc[key] = assoc[key] + value
            except KeyError:
                assoc[key] = value

def map_to_premir(assocs, param):
    """ map keys of the dict to premir form forms"""
    for assoc in assocs:
        to_delete = set()
        to_add = {}
        for key, value in assoc.items():
            if MirData.is_mature(key):
                to_delete.add(key)
                key = key[2:]
                mirs = MirData.mature2mir(key)
                for mir in mirs:
                    typed_mir = 'm:' + mir
                    try:
                        to_add[typed_mir] = to_add[typed_mir] | value
                    except TypeError:
                        to_add[typed_mir] = to_add[typed_mir] + value
                    except KeyError:
                        to_add[typed_mir] = value
        for item in to_delete:
            del assoc[item]
        for key, value in to_add.items():
            try:
                assoc[key] = assoc[key] | value
            except TypeError:
                assoc[key] = assoc[key] + value
            except KeyError:
                assoc[key] = value

        if param.GROUP_WHEN_SIMILAR_MATURE:
            to_delete = set()
            to_add = {}
            mirnas = assoc.keys()
            for key, value in assoc.items():
                togroup = set()
                matures = MirData.mir2mature(key)
                for mature in matures:
                    togroup.update(MirData.mature2mir(mature))
                togroup.intersection_update(mirnas)
                if len(togroup) == 1:
                    continue
                to_delete.update(togroup)
                group_name = 'm' + '/'.join(list(togroup))
                try:
                    to_add[group_name] = to_add[group_name] | set(value)
                except TypeError:
                    to_add[group_name] = to_add[group_name] + list(value)
                except KeyError:
                    to_add[group_name] = value
            for item in to_delete:
                del assoc[item]
            for key, value in to_add.items():
                assoc[key] = value

def get_similarities2(param):
    """transform for each miRNA, similarity data to a list of similar diseases"""
    with open(os.path.normpath(os.path.join(param.DATADIR, 'generated', 'disease-similarities.json')), 'r') as simil_file:
        simils = json.load(simil_file)
    ext_simil = {}
    for key, value in simils.items():
        if param.DISEASE_SIMILARITY_CUTOFF and value > param.DISEASE_SIMILARITY_CUTOFF[0]:
            pair = key.split('-')
            ext_simil[pair[0]] = ext_simil.get(pair[0], set([]))
            ext_simil[pair[1]] = ext_simil.get(pair[1], set([]))
        for cutoff in param.DISEASE_SIMILARITY_CUTOFF:
            if value > cutoff:
                ext_simil[pair[0]].add(pair[1]+'_'+str(cutoff))
                ext_simil[pair[1]].add(pair[0]+'_'+str(cutoff))
    return ext_simil

def create_dictionary(disease_map, word_map, target_map, neighbor_map, family_map, ext_simil, param):
    """create a dictionary from data stored in three matrixes"""
    all_text = []
    all_distinct_values = set()
    for key in disease_map.keys():
        values = disease_map.get(key).copy()
        if param.EXPAND_DISEASES:
            extend_disease(values)
        if param.USE_DISEASE_SIM_ON_DATA:
            extend_similarities(values, ext_simil)
        values.update(target_map.get(key, set()))
        values.update(neighbor_map.get(key, set()))
        values.update(family_map.get(key, set()))
        line = sorted(values)
        if key in word_map:
            line.extend(word_map[key])
        all_text.append(set(line))
        all_distinct_values.update(set(line))
    dictionary = corpora.Dictionary(all_text, prune_at=None)
#   dictionary.filter_extremes(no_below=5, no_above=0.8)
    return dictionary

def evaluate_results(disease_map, word_map, target_map, neighbor_map, family_map, param):
    """main loop for evaluating the method"""
    infirmed_associations = {}
    predicted_associations = {}
    results = {}

    if param.EVALUATION_FOR_ALL_TERMS:
        diseases = set([])
        nb_mirna_by_disease = {}
        for values in disease_map.values():
            diseases.update(values)
            for value in values:
                nb_mirna_by_disease[value] = nb_mirna_by_disease.get(value, 0) + 1
        diseases = sorted(diseases)
        set_of_targets = []
        for disease in diseases:
            if nb_mirna_by_disease[disease] >= param.EVALUATION_MIN_TERM_ASSOC:
                set_of_targets.append([disease])
#    set_of_targets = [["Lupus Nephritis"]]
    nb_assoc = [0 for _ in range(len(set_of_targets))]
    tested_target = set([x[0] for x in set_of_targets])
    for key, value in disease_map.items():
        for ctr in range(len(set_of_targets)):
            if set_of_targets[ctr][0] in value:
                nb_assoc[ctr] += 1
        if param.REMOVE_NOT_TESTED:
            value &= tested_target
#     set_of_targets=[["Breast Neoplasms"]]
    logging.info("associations %sim", nb_assoc)

    if param.REMOVE_OUTLIERS:
        remove_rare_or_frequent(disease_map)
    # remove rare diseases and frequent diseases
    nb_data = len(list(disease_map.keys()))
    nb_partitions = 5
    partition_keys = {}
    all_keys = set(disease_map.keys())
    logging.info("len all_keys= %d", len(all_keys))

    for ctr in range(nb_partitions):
        partition_keys[ctr] = []
        for _ in range(int(nb_data/nb_partitions)):
            key = random.choice(list(all_keys))
            partition_keys[ctr].append(key)
            all_keys.remove(key)
        logging.info(len(partition_keys[ctr]))
    logging.info(len(all_keys))
    logging.info(len(disease_map))

    ext_simil = {}
    if param.USE_DISEASE_SIM_ON_DATA:
        ext_simil = get_similarities2(param)

    assocmap_ext = word_map
    assocmap_ext2 = target_map
    assocmap_ext3 = neighbor_map
    assocmap_ext4 = family_map

    global_dictionary = create_dictionary(disease_map, word_map, target_map, neighbor_map, family_map, ext_simil, param)

    av_score = {}
    nb_not_scored = {}
    for ctr in range(nb_partitions):
        testing_keys = partition_keys[ctr]
        training_keys = set(disease_map.keys()) - set(testing_keys)
        training_keys.update(all_keys)

        training_corpus = build_corpus(disease_map, training_keys, global_dictionary, ext_simil, param=param)
        training_corpus_ext = build_corpus_ext(assocmap_ext, training_keys, global_dictionary)
        training_corpus_ext2 = build_corpus_ext2(assocmap_ext2, training_keys, global_dictionary)
        training_corpus_ext3 = build_corpus_ext2(assocmap_ext3, training_keys, global_dictionary)
        training_corpus_ext4 = build_corpus_ext2(assocmap_ext4, training_keys, global_dictionary)

        training_dictionary = global_dictionary

        if param.PROJECTION_MATRIX1_NBI:
            training_corpus = nbi.do_nbi_old(training_corpus, len(training_dictionary.keys()))

        if param.PROJECTION_MATRIX2_TFIDF:
            tfidf = models.TfidfModel(training_corpus_ext, normalize=True)
            training_corpus_mod = tfidf[training_corpus_ext]
            training_corpus_ext = []
            for row in training_corpus_mod:
                training_corpus_ext.append(row)

        if param.PROJECTION_MATRIX3_NBI:
            training_corpus_ext2 = nbi.do_nbi_old(training_corpus_ext2, len(training_dictionary.keys()))

        for ctr in range(len(training_corpus)):
            dict1 = dict(training_corpus[ctr])
            dict1.update(dict(training_corpus_ext[ctr]))
            dict1.update(dict(training_corpus_ext2[ctr]))
            dict1.update(dict(training_corpus_ext3[ctr]))
            dict1.update(dict(training_corpus_ext4[ctr]))
            training_corpus[ctr] = list(dict1.items())

        lsi = models.LsiModel(training_corpus, id2word=training_dictionary, num_topics=param.LSI_DIMENSION) # initialize an LSI transformation
        average_auc = 0
        testing_corpus_ext = build_corpus_ext(assocmap_ext, testing_keys, global_dictionary)#all_distinct_values)
        testing_corpus_ext = [[] for x in testing_keys] # do not use litterature information
        testing_corpus_ext2 = build_corpus_ext2(assocmap_ext2, testing_keys, global_dictionary)
        testing_corpus_ext3 = build_corpus_ext2(assocmap_ext3, testing_keys, global_dictionary)
        testing_corpus_ext4 = build_corpus_ext2(assocmap_ext4, testing_keys, global_dictionary)

        if param.PROJECTION_MATRIX3_NBI:
            testing_corpus_ext2 = nbi.do_nbi_old(testing_corpus_ext2, len(training_dictionary.keys()))

        if param.PROJECTION_MATRIX2_TFIDF:
            tfidf = models.TfidfModel(testing_corpus_ext, normalize=True)
            testing_corpus_mod = tfidf[testing_corpus_ext]
            testing_corpus_ext = []
            for row in testing_corpus_mod:
                testing_corpus_ext.append(row)
        for env in set_of_targets:
            query_key = "#".join(sorted(env))

            confirmed = set()
            ignored = set()
            for key in testing_keys:
                query = set(env)
                if query.intersection(disease_map[key]):
                    confirmed.add(key)
            if not confirmed:
                logging.warning("%sim no_assoc", env[0])
                nb_not_scored[query_key] = nb_not_scored.get(query_key, 0) + 1
                continue

            testing_corpus = build_corpus(disease_map, testing_keys, global_dictionary, ext_simil, query_values=set(env), param=param)#all_distinct_values)

            if param.PROJECTION_MATRIX1_NBI:
                testing_corpus = nbi.do_nbi_old(testing_corpus, len(training_dictionary.keys()))

            for ctr in range(len(testing_corpus)):
                dict1 = dict(testing_corpus[ctr])
                dict1.update(dict(testing_corpus_ext2[ctr]))
                dict1.update(dict(testing_corpus_ext3[ctr]))
                dict1.update(dict(testing_corpus_ext4[ctr]))
                testing_corpus[ctr] = dict1.items()

#             index = similarities.MatrixSimilarity(lsi[testing_corpus])
            index = similarities.Similarity("./tmp", lsi[testing_corpus], num_features=param.LSI_DIMENSION)
            testing_dictionary = global_dictionary
            sims = lsi_search(env, testing_dictionary, index, lsi)

            nb_pos = len(confirmed)
            nb_neg = len(sims) - nb_pos
            true_pos = 0
            false_pos = 0
            true_neg = nb_neg
            false_neg = nb_pos
            ctr = 0
            true_neg = 0
            false_neg = 0
            true_pos = 0
            false_pos = 0
            ctr = 0
            for sim in sims:
                key = testing_keys[sim[0]]
                if key in ignored:
                    continue
                if key in confirmed:
                    false_neg += 1
                else:
                    true_neg += 1
                ctr += 1
                if ctr == 10000:
                    break
            ctr = 0
            tpr = []
            fpr = []
            positives = []
            for sim in sims:
                key = testing_keys[sim[0]]
                if key in ignored:
                    continue
                if key in confirmed:
                    true_pos += 1
                    false_neg -= 1
                    positives.append(1)
                else:
                    false_pos += 1
                    true_neg -= 1
                    positives.append(0)
                true_positive_rate = true_pos / (true_pos + false_neg)
                false_positive_rate = false_pos / (true_neg + false_pos)
                tpr.append(true_positive_rate)
                fpr.append(false_positive_rate)
                ctr += 1
                if ctr == 10000:
                    break

            average = metrics.auc(fpr, tpr)

            logging.info("positives=%d, %sim", nb_pos, positives)
#             for sim in sims[-1:-10:-1]:
#                 if testing_keys[sim[0]] in confirmed:
#                     logging.info("last confirmed=%sim", testing_keys[sim[0]])
            if average > 0.8: #only make a prediction if the AUC is good enough
                for ctr in range(len(sims)):
                    sim = sims[ctr]
                    if fpr[ctr] > param.FPR_CUTOFF and testing_keys[sim[0]] in confirmed:
                        try:
                            infirmed_associations[query_key].append(testing_keys[sim[0]])
                        except KeyError:
                            infirmed_associations[query_key] = [testing_keys[sim[0]]]
                    if tpr[ctr] < (1 - param.TPR_CUTOFF) and testing_keys[sim[0]] not in confirmed:
                        try:
                            predicted_associations[query_key].append(testing_keys[sim[0]])
                        except KeyError:
                            predicted_associations[query_key] = [testing_keys[sim[0]]]
            mirnas = [testing_keys[x[0]] for x in sims]
            encode = zip(mirnas, positives)
            try:
                results[query_key].append(encode)
            except KeyError:
                results[query_key] = [encode]
#             print(fpr)
#             print(tpr)
#             dd
            #plt.plot(fpr, tpr)
            #plt.savefig('myfig')
            logging.info("%sim average=%f", env, average)
            average_auc += average
            av_score[query_key] = av_score.get(query_key, 0) + average

        if set_of_targets:
            logging.info("average AUC = %f", average_auc/len(set_of_targets))
    global_score = 0
    results_file = os.path.join(param.OUTPUTDIR, 'results' + param.TIMESTAMP)
    with open(results_file + '.score', 'w') as res_file:
        res_file.write("query\tnb_assoc\tscore\n")
        for ctr in range(len(set_of_targets)):
            env_key = "#".join(sorted(set_of_targets[ctr]))
            print(env_key, "nb_assoc=", nb_assoc[ctr], "score=", av_score[env_key]/(nb_partitions - nb_not_scored.get(env_key, 0)))
            res_file.write(env_key + '\t' + str(nb_assoc[ctr]) + '\t' + str(av_score[env_key]/nb_partitions) + '\n')
            global_score += av_score[env_key]/nb_partitions
    with open(results_file + '-trace' + '.yaml', 'w') as res_file:
        yaml.dump(results, res_file, indent=2)
    with open(results_file + '-infirmed' + '.json', 'w') as res_file:
        json.dump(infirmed_associations, res_file, indent=2)
    with open(results_file + '-predicted' + '.json', 'w') as res_file:
        json.dump(predicted_associations, res_file, indent=2)

    logging.info("global score = %f", global_score/len(set_of_targets))
    return global_score/len(set_of_targets)
        
def predict_association(disease_map, word_map, target_map, neighbor_map, family_map, param):
    """main loop for predicting associations"""

    queries = set([])
    nb_assoc = {}
    for values in disease_map.values():
        queries.update(values)
        for value in values:
            nb_assoc[value] = nb_assoc.get(value, 0) +1
    queries = sorted([[x] for x in queries])

    if os.path.isfile(param.MODEL_FILE + '.lsi'):
        lsi = models.LsiModel.load(param.MODEL_FILE + '.lsi')
        dictionary = corpora.Dictionary.load(param.MODEL_FILE + '.dict')
        index = similarities.MatrixSimilarity.load(param.MODEL_FILE + '.index')

        training_corpus = corpora.MmCorpus(param.MODEL_FILE + '.mm')
        with open(param.MODEL_FILE + '.json', 'r') as keys_file:
            training_keys = json.load(keys_file)
    else:
        training_keys = list(disease_map.keys())

        if param.REMOVE_OUTLIERS:
            remove_rare_or_frequent(disease_map)

        ext_simil = {}
        if param.USE_DISEASE_SIM_ON_DATA:
            ext_simil = get_similarities2(param)

        # create a dictionary
        dictionary = create_dictionary(disease_map, word_map, target_map, neighbor_map, family_map, ext_simil, param)

        training_corpus = build_corpus(disease_map, training_keys, dictionary, ext_simil, param=param)
        training_corpus_ext = build_corpus_ext(word_map, training_keys, dictionary)
        training_corpus_ext2 = build_corpus_ext2(target_map, training_keys, dictionary)
        training_corpus_ext3 = build_corpus_ext2(neighbor_map, training_keys, dictionary)
        training_corpus_ext4 = build_corpus_ext2(family_map, training_keys, dictionary)

        if param.PROJECTION_MATRIX1_NBI:
            training_corpus = nbi.do_nbi_old(training_corpus, len(dictionary.keys()))

        if param.PROJECTION_MATRIX3_NBI:
            training_corpus_ext2 = nbi.do_nbi_old(training_corpus_ext2, len(dictionary.keys()))

        if param.PROJECTION_MATRIX2_TFIDF:
            tfidf = models.TfidfModel(training_corpus_ext, normalize=True)
            training_corpus_mod = tfidf[training_corpus_ext]
            training_corpus_ext = []
            for row in training_corpus_mod:
                training_corpus_ext.append(row)

        for ctr in range(len(training_corpus)):
            dict1 = dict(training_corpus[ctr])
            dict1.update(dict(training_corpus_ext[ctr]))
            dict1.update(dict(training_corpus_ext2[ctr]))
            dict1.update(dict(training_corpus_ext3[ctr]))
            dict1.update(dict(training_corpus_ext4[ctr]))
            training_corpus[ctr] = dict1.items()

        # initialize an LSI transformation
        lsi = models.LsiModel(training_corpus, id2word=dictionary, num_topics=param.LSI_DIMENSION)
        index = similarities.MatrixSimilarity(lsi[training_corpus])

        dictionary.save(param.MODEL_FILE + '.dict')
        index.save(param.MODEL_FILE + '.index')
        corpora.MmCorpus.serialize(param.MODEL_FILE + '.mm', training_corpus)
        lsi.save(param.MODEL_FILE + '.lsi') # same for tfidf, lda, ...
        with open(param.MODEL_FILE + '.json', 'w') as keys_file:
            json.dump(training_keys, keys_file, indent=2)

    cutoff = 0.8
    max_lookup = 100

    all_annotated = set()
    for key, value in disease_map.items():
        if value:
            all_annotated.add(key)

    predicted_associations = {}
    for diseases in queries:
        confirmed = set([])
        query = set(diseases)
        if param.USE_DISEASE_SIM_ON_QUERY:
            for disease in diseases:
                query.update(MeSHData.get_childs(disease))
                query.add(disease)
#         for key in training_keys:
#             try:
#                 if query.intersection(disease_map[key]):
#                     confirmed.add(key)
#             except KeyError:
#                 pass
        for key, value in disease_map.items():
#             print(key,value, query.intersection(value))
            if query.intersection(value):
                confirmed.add(key)

        sims = lsi_search(query, dictionary, index, lsi)

        # compute the index of the last sims > cutoff
        ctr = 0
        true_pos = 0
        false_pos = 0
        positive = len(confirmed)
        negative = len(sims) - len(confirmed)
        last_index = 0
        tpr = []
        fpr = []
        results = []
        for sim in sims:
            key = training_keys[sim[0]]

            if key in confirmed:
                true_pos += 1
                results.append(1)
            else:
                false_pos += 1
                results.append(0)

            if ctr > max_lookup:
                results = results[-max_lookup:]
            if sum(results) / len(results) >= cutoff:
#             if true_pos / (ctr + 1) >= cutoff:
                last_index = ctr
            ctr += 1
            true_positive_rate = true_pos / positive
            false_positive_rate = false_pos / negative
            tpr.append(true_positive_rate)
            fpr.append(false_positive_rate)
        average = metrics.auc(fpr, tpr)
        print(diseases, "AUC=", average, " confirmed=", len(confirmed))


        last_index += 10
        true_pos = 0

        predicted_accepted = []
        predicted_potential = []
        ctr = 0
        true_pos = 0
        false_pos = 0
        last_index = 1000
        new_associations = []
        for sim in sims:
            key = training_keys[sim[0]]
            if key in confirmed:
                true_pos += 1

                confirmed.remove(key)
                print(ctr, key, " confirmed ", sim[1], true_pos / positive, false_pos / negative, len(confirmed))
                predicted_accepted.extend(predicted_potential)
                predicted_potential = []
            else:
                false_pos += 1
                print(ctr, key, " ........ ", sim[1], true_pos / positive, false_pos / negative)
                predicted_potential.append((key, true_pos / positive))
                if false_pos / negative < 0.05 and true_pos / positive < 0.90:
                    new_associations.append(key)

            if last_index == ctr:
                break
            ctr += 1
        if new_associations:
            predicted_associations[diseases[0]] = new_associations
#     outfile = open('predict_mirna2disease.txt', 'w')
    with open(os.path.join(param.OUTPUTDIR, 'predicted_associations.json'), 'w') as prediction_file:
        json.dump(predicted_associations, prediction_file, indent=2)
#              outfile.write("{0}\t".format(diseases))
#             outfile.write("{0}\n".format(' '.join(sorted(new_associations))))

def add_empty_coordinates(bow, dim):
    """
    add coordinates with 0 values

    by default, Gensim does not keep coordinates with zero value
    this function add them in the bow
    """
    dict_bow = dict(bow)
    for coord in range(dim):
        if coord not in dict_bow:
            dict_bow[coord] = 0
    return list(dict_bow.items())

def lsi_search_term2doc(terms, dictionary, corpus, lsi, param):
    """reimplementation of lsi_search"""
    logging.debug(terms)
    if len(terms) == 1:
        vec_bow = dictionary.doc2bow(terms)
        gensim_coord = add_empty_coordinates(lsi[vec_bow], param.LSI_DIMENSION)
        base_coord = np.array([x[1] for x in gensim_coord])
    else:
        nb_coord = 0
        base_coord = np.array([])
        for term in terms:
            vec_bow = dictionary.doc2bow([term])
            if vec_bow:
                nb_coord += 1
                gensim_coord = add_empty_coordinates(lsi[vec_bow], param.LSI_DIMENSION)
                if base_coord.size == 0:
                    base_coord = np.array([x[1] for x in gensim_coord])
                else:
                    base_coord = base_coord + np.array([x[1] for x in gensim_coord])
        base_coord = base_coord / nb_coord
    sims = []
    coord = matutils.corpus2dense(lsi[corpus], param.LSI_DIMENSION)
    for ctr in range(coord.shape[1]):
        vec = coord[:, ctr]
        simil = 1 - scipy.spatial.distance.cosine(base_coord, vec)
        sims.append(simil)
    return sorted(enumerate(sims), key=lambda item: -item[1])

def lsi_search_doc2term(docs, dictionary, corpus, lsi, param):
    """inverse search: return list of terms similar to docs"""
    if len(docs) == 1:
        gensim_coord = add_empty_coordinates(lsi[corpus[docs[0]]], param.LSI_DIMENSION)
        base_coord = np.array([x[1] for x in gensim_coord])
    else:
        nb_coord = 0
        base_coord = np.array([])
        for doc in docs:
            nb_coord += 1
            gensim_coord = add_empty_coordinates(lsi[corpus[doc]], param.LSI_DIMENSION)
            if base_coord.size == 0:
                base_coord = np.array([x[1] for x in gensim_coord])
            else:
                base_coord = base_coord + np.array([x[1] for x in gensim_coord])
        base_coord = base_coord / nb_coord

    sims = []
    for gensim_id, word in dictionary.iteritems():
        vec_bow = dictionary.doc2bow([word])
        gensim_coord = add_empty_coordinates(lsi[vec_bow], param.LSI_DIMENSION)
        coord = [x[1] for x in gensim_coord]
        simil = 1 - scipy.spatial.distance.cosine(base_coord, coord)
        sims.append((gensim_id, simil))
    return sorted(sims, key=lambda item: -item[1])

def lsi_search_doc2term2(docs, dictionary, corpus, lsi):
    """return list of terms similar to docs
    version using gensim methods"""
    termcorpus = matutils.Dense2Corpus(lsi.projection.u.T)
    index = similarities.MatrixSimilarity(termcorpus)
    sims = index[lsi[corpus[docs[0]]]]
    return sorted(enumerate(sims), key=lambda item: -item[1])

def read_mir2disease(param):
    """read mir2disease data"""
    MeSHData.alias2mesh.update(DiseaseData.add_hmdd_aliases())
    mir2disease = pd.read_csv(param.MIR2DISEASE_FILES[param.MIR2DISEASE_SOURCE], delimiter='\t', usecols=[1, 2, 3], skipinitialspace=True, header=None, names=['miRNA', 'disease', 'pmid'])
    mir2disease = mir2disease.drop('pmid', axis=1)

    # keep only human miRNA
    mir2disease = mir2disease[mir2disease.miRNA.str[:3] == 'hsa']
    mir2disease['miRNA'] = mir2disease['miRNA'].apply(MirData.add_mir_type)

    # remove dead or unknown
    mir2disease = mir2disease[mir2disease.miRNA.str[:2] != 'd:']
    mir2disease = mir2disease[mir2disease.miRNA.str[:2] != 'u:']
    mir2disease['disease'] = mir2disease['disease'].map(MeSHData.convert2mesh)

    not_mapped = mir2disease[mir2disease.disease.str[:11] == 'not_mapped:']
    not_mapped = set(not_mapped.loc[:, 'disease'].unique())

    all_diseases = set(mir2disease.loc[:, 'disease'].unique())
    all_mirna = set(mir2disease.loc[:, 'miRNA'].unique())
#     print(all_mirna)
    mir2disease = mir2disease.set_index('miRNA')
    logging.info("number of miRNAs: %d", len(all_mirna))
    logging.info("number of diseases: %d", len(all_diseases))
    logging.info("not mapped diseases = %d out of %d: %s", len(not_mapped), len(all_diseases), not_mapped)
    return mir2disease

def read_mir2env(param):
    """read mir2env data"""
    mir2env = pd.read_csv(param.MIR2ENV_FILES[param.MIR2ENV_SOURCE], delimiter='\t', quotechar='\t', quoting=2, usecols=[1, 2, 5], skipinitialspace=True, header=0, names=['miRNA', 'product', 'env'])

    mir2env = mir2env[mir2env.miRNA.str[:3] == 'hsa'] # keep only humain miRNA

    mir2env['miRNA'] = mir2env['miRNA'].apply(MirData.add_mir_type)
    # remove dead or unknown
    mir2env = mir2env[mir2env.miRNA.str[:2] != 'd:']
    mir2env = mir2env[mir2env.miRNA.str[:2] != 'u:']

    all_env = set(mir2env.loc[:, 'env'].unique())

    mir2env['env'] = mir2env['env'].map(MeSHData.convert2mesh)
    not_mapped = mir2env[mir2env.env.str[:11] == 'not_mapped:']
    not_mapped = set(not_mapped.loc[:, 'env'].unique())

    all_mirna = set(mir2env.loc[:, 'miRNA'].unique())
    if param.USE_MATURE_FORM:
        mir2env = mir2env.set_index('product')
    else:
        mir2env = mir2env.set_index('miRNA')    
    logging.info("number of miRNAs: %d", len(all_mirna))
    logging.info("number of environments: %d", len(all_env))
    logging.info("not mapped environment = %d out of %d: %s", len(not_mapped), len(all_env), not_mapped)
    return mir2env

def read_mir2env_disease(param):
    """read mirendata to collect environments and diseases associated to miRNAs"""
    mir2env = pd.read_csv(param.MIR2ENV_FILES[param.MIR2ENV_SOURCE], delimiter='\t', quotechar='\t', quoting=2, usecols=[1, 2, 4, 5], skipinitialspace=True, header=0, names=['miRNA', 'product', 'disease', 'env'])

    mir2env = mir2env[mir2env.miRNA.str[:3] == 'hsa'] # keep only humain miRNA

    mir2env['miRNA'] = mir2env['miRNA'].apply(MirData.add_mir_type)
    # remove dead or unknown
    mir2env = mir2env[mir2env.miRNA.str[:2] != 'd:']
    mir2env = mir2env[mir2env.miRNA.str[:2] != 'u:']

    all_mirna = set(mir2env.loc[:, 'miRNA'].unique())
    logging.info("number of miRNAs: %d", len(all_mirna))

    all_env = set(mir2env.loc[:, 'env'].unique())
    mir2env['env'] = mir2env['env'].map(MeSHData.convert2mesh)
    not_mapped = mir2env[mir2env.env.str[:11] == 'not_mapped:']
    not_mapped = set(not_mapped.loc[:, 'env'].unique())
    logging.info("number of environments: %d", len(all_env))
    logging.info("not mapped environment = %d out of %d: %s", len(not_mapped), len(all_env), not_mapped)

    all_diseases = set(mir2env.loc[:, 'disease'].unique())
    mir2env['disease'] = mir2env['disease'].map(MeSHData.convert2mesh)
    not_mapped = mir2env[mir2env.disease.str[:11] == 'not_mapped:']
    not_mapped = set(not_mapped.loc[:, 'disease'].unique())
    logging.info("number of diseases: %d", len(all_diseases))
    logging.info("not mapped diseases = %d out of %d: %s", len(not_mapped), len(all_diseases), not_mapped)

    if param.USE_MATURE_FORM:
        mir2env = mir2env.set_index('product')
    else:
        mir2env = mir2env.set_index('miRNA')    
    return mir2env

def read_disease2env(param):
    """read disease2env data"""
    disease2env = pd.read_csv(param.DISEASE2ENV_FILES[param.DISEASE2ENV_SOURCE], compression='gzip', delimiter='\t', comment='#', usecols=[0, 3], skipinitialspace=True, header=0, names=['env', 'disease'])

    all_env = set(disease2env.loc[:, 'env'].unique())
    all_disease = set(disease2env.loc[:, 'disease'].unique())

    disease2env['env'] = disease2env['env'].map(MeSHData.convert2mesh)
    not_mapped = disease2env[disease2env.env.str[:11] == 'not_mapped:']
    not_mapped = set(not_mapped.loc[:, 'env'].unique())
    logging.info("not mapped environment = %d out of %d: %s", len(not_mapped), len(all_env), not_mapped)

    disease2env['disease'] = disease2env['disease'].map(MeSHData.convert2mesh)
    not_mapped = disease2env[disease2env.disease.str[:11] == 'not_mapped:']
    not_mapped = set(not_mapped.loc[:, 'disease'].unique())
    logging.info("not mapped diseases = %d out of %d: %s", len(not_mapped), len(all_disease), not_mapped)

    disease2env = disease2env.set_index('disease')
    return disease2env

def load_data(param):
    """ load data from files"""
    mir2disease = read_mir2disease(param)
    disease_map = common_utils.create_assoc_map(mir2disease, 'disease')
    extend_families(disease_map)
    logging.info('Number of miRNAs after family extension: %d', len(disease_map))

    word_map = create_word_map(set(disease_map.keys()), param)
    target_map = create_target_map(set(disease_map.keys()), param)
    neighbor_map = create_neighbor_map(set(disease_map.keys()), param)
    family_map = create_family_map(set(disease_map.keys()), param)

    dicts = (disease_map, word_map, target_map, neighbor_map, family_map)
    map_to_premir(dicts, param)
    logging.debug("disease_map:   %s", common_utils.get_stats_on_assocmap(disease_map))
    logging.debug("word_map: %s", common_utils.get_stats_on_assocmap(word_map))
    logging.debug("target_map: %s", common_utils.get_stats_on_assocmap(target_map))
    logging.debug("neighbor_map: %s", common_utils.get_stats_on_assocmap(neighbor_map))
    logging.debug("family_map: %s", common_utils.get_stats_on_assocmap(family_map))
    


    return disease_map, word_map, target_map, neighbor_map, family_map

def main():
    """main entry point"""
    param = MiRAIParam()
    MeSHData.read_data()
    disease_map, word_map, target_map, neighbor_map, family_map = load_data(param)
    
    if param.PERFORM_EVALUATION:
#        random.seed(0)
        score = 0
        nb_run = 1
        for _ in range(nb_run):
            score += evaluate_results(disease_map, word_map, target_map, neighbor_map, family_map, param)
        logging.info('final score=%f', score/nb_run)

    if param.PERFORM_PREDICTION:
        predict_association(disease_map, word_map, target_map, neighbor_map, family_map, param)

if __name__ == '__main__':
    logging.basicConfig(format='%(asctime)s : %(levelname)s : %(message)s', level=logging.INFO)
    if len(sys.argv) > 1:
        MiRAIParam.TIMESTAMP = sys.argv[1]
    main()

