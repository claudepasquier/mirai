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

   this file contains classes and methods related to disease information
"""

__author__ = "Claude Pasquier"
__copyright__ = "2015-2017 UNS-CNRS"
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Claude Pasquier"
__email__ = "claude.pasquier@unice.fr"
__status__ = "Development"
__all__ = ['DiseaseData']

class DiseaseData:
    """ group disease related data """

    @staticmethod
    def add_mircancer_aliases():
        """add data for mirCancer database"""
        alias2mesh = {}
        alias2mesh['diffuse large B-cell lymphoma'] = 'Lymphoma, Large B-Cell, Diffuse'
        alias2mesh['acute leukemia'] = 'Leukemia'
        alias2mesh['malignant mesothelioma'] = 'Mesothelioma'
        alias2mesh['primary cns lymphomas'] = 'Central Nervous System Neoplasms'
        alias2mesh['medullary thyroid carcinoma'] = 'Thyroid Neoplasms'
        alias2mesh['clear cell renal cell cancer'] = 'Carcinoma, Renal Cell'
        alias2mesh['monocytic leukemia'] = 'Leukemia, Monocytic, Acute'
        alias2mesh['splenic marginal zone lymphoma'] = 'Lymphoma, B-Cell'
        alias2mesh['mesenchymal cancer'] = 'Epithelial-Mesenchymal Transition'
        alias2mesh['follicular cancer'] = 'Thyroid Neoplasms'
        alias2mesh['renal clear cell carcinoma'] = 'Carcinoma, Renal Cell'
        alias2mesh['colon carcinoma'] = 'Colonic Neoplasms'
        alias2mesh['oral carcinoma'] = 'Mouth Neoplasms'
        alias2mesh['gallbladder carcinoma'] = 'Gallbladder Neoplasms'
        alias2mesh['prostate carcinoma'] = 'Prostatic Neoplasms'
        alias2mesh['lung squamous cell carcinoma'] = 'Carcinoma, Non-Small-Cell Lung'
        alias2mesh['esophageal carcinoma'] = 'Esophageal Neoplasms'
        alias2mesh['laryngeal squamous cell carcinoma'] = 'Laryngeal Neoplasms'
        alias2mesh['colonic adenocarcinoma'] = 'Colonic Neoplasms'
        alias2mesh['primary gallbladder carcinoma'] = 'Gallbladder Neoplasms'
        alias2mesh['head and neck squamous cell carcinoma'] = 'Squamous Cell Carcinoma, Head And Neck'
        alias2mesh['oral squamous cell carcinoma'] = 'Mouth Neoplasms'
        alias2mesh['laryngeal carcinoma'] = 'Laryngeal Neoplasms'
        alias2mesh['cervical carcinoma'] = 'Uterine Cervical Neoplasms'
        alias2mesh['esophageal adenocarcinoma'] = 'Esophageal Neoplasms'
        alias2mesh['primary thyroid lymphoma'] = 'Thyroid Neoplasms'
        alias2mesh['endometrial serous adenocarcinoma'] = 'Endometrial Neoplasms'
        alias2mesh['cervical squamous cell carcinoma'] = 'Uterine Cervical Neoplasms'
        alias2mesh['colorectal adenocarcinoma'] = 'Colorectal Neoplasms'
        alias2mesh['hypopharyngeal squamous cell carcinoma'] = 'Hypopharyngeal Neoplasms'
        alias2mesh['pancreatic ductal adenocarcinoma'] = 'Carcinoma, Pancreatic Ductal'
        alias2mesh['pancreatic adenocarcinoma'] = 'Pancreatic Neoplasms'
        alias2mesh['gastric adenocarcinoma'] = 'Stomach Neoplasms'
        alias2mesh['ovarian carcinoma'] = 'Ovarian Neoplasms'
        return alias2mesh

    @staticmethod
    def add_hmdd_aliases():
        """add data for HMDD database"""
        alias2mesh = {}
        alias2mesh['Ovary Syndrome'] = 'Polycystic Ovary Syndrome'
        alias2mesh['HBV Infection'] = 'Hepatitis B virus'
        alias2mesh['Carcinoma, Oral'] = 'Mouth Neoplasms'
        alias2mesh['HEV'] = 'Hepatitis E virus'
        alias2mesh['HPV Infection'] = 'Papillomavirus Infections'
        alias2mesh['Papilary thyroid carcinoma'] = 'Thyroid Neoplasms|Carcinoma, Papillary'
        alias2mesh['PRRSV Infection'] = 'Porcine respiratory and reproductive syndrome virus'
        alias2mesh['Fatty Liver, Non-Alcoholic'] = 'Non-alcoholic Fatty Liver Disease'
        alias2mesh['SIV Infection'] = 'Simian immunodeficiency virus'
        alias2mesh['Leukemia, Acute'] = 'Leukemia, Acute'
        alias2mesh['HCMV Infection'] = 'Cytomegalovirus Infections'
        alias2mesh['HCV'] = 'Hepacivirus'
        return alias2mesh

    @staticmethod
    def add_mir2disease_aliases():
        """add data for mir2disease database"""
        alias2mesh = {}
        alias2mesh['limb-girdle muscular dystrophies types 2A (LGMD2A)'] = 'Muscular Dystrophies, Limb-Girdle'
        alias2mesh['Head and neck cancer'] = 'Head and Neck Neoplasms'
        alias2mesh['head and neck squamous cell carcinoma (HNSCC)'] = 'Head and Neck Neoplasms'
        alias2mesh['pancreatic ductal adenocarcinoma (PDAC)'] = 'Carcinoma, Pancreatic Ductal'
        alias2mesh['type 2 diabetes'] = 'Diabetes Mellitus, Type 2'
        alias2mesh['HCV infection'] = 'Hepatitis C'
        alias2mesh['epithelial ovarian cancer (EOC)'] = 'Ovarian epithelial cancer'
        alias2mesh['non-alcoholic fatty liver disease (NAFLD)'] = 'Non-alcoholic Fatty Liver Disease'
        alias2mesh['neutrophilia'] = 'Leukocyte Disorders'
        alias2mesh['Intrahepatic cholangiocarcinoma (ICC)'] = 'Intrahepatic cholangiocarcinoma'
        alias2mesh["Huntington's disease (HD)"] = 'Huntington Disease'
        alias2mesh['HBV-related cirrhosis'] = 'Liver Cirrhosis'
        alias2mesh["Hodgkin's lymphoma"] = 'Hodgkin Disease'
        alias2mesh['acute myelogeneous leukemia (AML)'] = 'Leukemia, Myeloid'
        alias2mesh['serous ovarian cancer'] = 'Ovarian Neoplasms'
        alias2mesh['papillary thyroid carcinoma (PTC)'] = 'Thyroid cancer, papillary'
        alias2mesh['prostate cance'] = 'Prostatic Neoplasms'
        alias2mesh["Alzheimer's disease"] = 'Alzheimer Disease'
        alias2mesh['nasopharyngeal carcinoma (NPC)'] = 'Nasopharyngeal carcinoma'
        alias2mesh['homozygous sickle cell disease (HbSS)'] = 'Anemia, Sickle Cell'
        alias2mesh['PFV-1 infection'] = 'Spumavirus'
        alias2mesh['Oral Squamous Cell Carcinoma (OSCC)'] = 'Mouth Neoplasms'
        alias2mesh['Malignant mesothelioma (MM)'] = 'Mesothelioma, Malignant'
        alias2mesh['neurodegeneration'] = 'Nerve Degeneration'
        alias2mesh['tongue squamous cell carcinoma'] = 'Tongue Neoplasms'
        alias2mesh['recurrent ovarian cancer'] = 'Ovarian Neoplasms'
        alias2mesh['uterine leiomyoma (ULM)'] = 'Leiomyoma'
        alias2mesh["Parkinson's disease"] = 'Parkinson Disease'
        alias2mesh['miyoshi myopathy (MM)'] = 'Miyoshi myopathy'
        alias2mesh['diffuse large B-cell lymphoma (DLBCL)'] = 'Lymphoma, Large B-Cell, Diffuse'
        alias2mesh['Cerebellar neurodegeneration'] = 'Paraneoplastic Cerebellar Degeneration'
        alias2mesh['myocardial injury'] = 'Cardiomyopathies'
        alias2mesh['Glomerulosclerosis'] = 'Glomerulonephritis'
        alias2mesh['diarrhea predominant irritable bowel syndrome (IBS-D)'] = 'Irritable Bowel Syndrome'
        alias2mesh['Oral Carcinoma'] = 'Mouth Neoplasms'
        alias2mesh['glomerular disease'] = 'Glomerulonephritis'
        return alias2mesh

if __name__ == '__main__':
    pass

