### MiRAI
Prediction of miRNA-disease associations with a vector space model

11 January 2019

Copyright 2015-2019 UNS-CNRS

Author: Claude Pasquier (I3S Laboratory, UCA)

Contact: claude.pasquier@unice.fr 

### Preamble

MiRAI is a method that uses available data produced by biologists to perform predictions.
We first applied it to the micro-RNA topic.

Our basic approach is to represent distributional information on miRNAs and diseases in
a high-dimensional vector space and to define associations between miRNAs and diseases
in terms of their vector similarity. Parameters of MiRAI were tuned using an evolutionary
algorithm. Cross validations performed on a dataset of known miRNA-disease associations
demonstrate the excellent performance of our method.

### Prerequisite
MiRAI needs python3 and uses the following libraries:
* **boto**   MIT   https://pypi.python.org/pypi/boto 
* **bz2file**  Apache 2.0 https://pypi.python.org/pypi/bz2file 
* **certifi** MPL https://github.com/certifi/python-certifi/blob/master/LICENSE 
* **cycler** BSD https://pypi.python.org/pypi/Cycler 
* **chardet** GNU LGPLv2.1 https://pypi.python.org/pypi/chardet 
* **gensim** GNU   LGPLv2.1 https://radimrehurek.com/gensim/about.html
* **idna** BSD like https://pypi.python.org/pypi/idna 
* **matplotlib** Python Software Foundation (PSF) http://matplotlib.org/
* **numpy**  BSD-3 https://docs.scipy.org/doc/numpy-1.10.0/license.html 
* **pandas** BSD-3 http://pandas.pydata.org/
* **pyparsing** MIT https://pypi.python.org/pypi/pyparsing/2.2.0 
* **python-dateutil** BSD-3 https://pypi.python.org/pypi/python-dateutil/2.6.0 
* **pytz** MIT https://pypi.python.org/pypi/pytz 
* **PyYAML** MIT http://pyyaml.org/
* **requests** Apache 2.0 http://docs.python-requests.org/en/master/user/intro/ 
* **scikit-learn** BSD http://scikit-learn.org
* **smart-open** MIT https://pypi.python.org/pypi/smart_open 
* **scipy**  BSD-3 https://www.scipy.org/scipylib/license.html 
* **six**    MIT   https://pypi.python.org/pypi/six/ 
* **urllib3** MIT https://urllib3.readthedocs.io/en/latest/ 

It also uses the following databases:
* **mirtarbase** Academic use free of charge http://mirtarbase.mbc.nctu.edu.tw/

  BECAUSE MIRTARBASE DATA ARE LICENSED FOR ACADEMIC USE FREE OF CHARGE, IN 
  NO EVENT SHALL ANY MEMBER OF THE MIRTARBASE DEVELOPMENT TEAM BE LIABLE TO
  ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL
  DAMAGES ARISING OUT OF THE USE OF THE MIRTARBASE DATABASE, EVEN IF THE
  MIRTARBASE TEAM HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES [...].
  (see http://mirtarbase.mbc.nctu.edu.tw/cache/download/LICENSE for details)

* **MIRbase** http://www.mirbase.org/

  miRBase is in the public domain. It is not copyrighted.  You may
  freely modify, redistribute, or use it for any purpose.  See
  ftp://mirbase.org/pub/mirbase/CURRENT/LICENSE for details, and
  http://www.mirbase.org/ for information on how to cite miRBse in
  publications.

* **MeSH** https://www.nlm.nih.gov/mesh
  Downloading data from the National Library of Medicine FTP servers
  indicates your acceptance of the following Terms and Conditions:
  No charges, usage fees or royalties are paid to NLM for this data.

  NLM freely provides PubMed/MEDLINE data. Please note some PubMed/MEDLINE
  abstracts may be protected by copyright
  (see https://www.nlm.nih.gov/databases/download/terms_and_conditions.html for details)
  
* **Pubmed** https://www.ncbi.nlm.nih.gov/pubmed
  Downloading data from the National Library of Medicine FTP servers
  indicates your acceptance of the following Terms and Conditions:
  No charges, usage fees or royalties are paid to NLM for this data.

  NLM freely provides PubMed/MEDLINE data. Please note some PubMed/MEDLINE
  abstracts may be protected by copyright
  (see https://www.nlm.nih.gov/databases/download/terms_and_conditions.html for details)
    
* **HMDD v2.0** http://www.cuilab.cn/hmdd
  all data in the database can be freely downloaded.
  (See https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3964961/ for details) 

### Running the program

1. Install the following libraries:
    * NumPy >= 1.3 http://www.numpy.org/
    * SciPy >= 0.7 https://www.scipy.org/
    * scikit-learn >= v0.18.1  http://scikit-learn.org
    * pandas >= v0.20.1  http://pandas.pydata.org/
    * PyYAML >= v3.12    http://pyyaml.org/
    * matplotlib >= v2.0.2 http://matplotlib.org/
    * gensim >= v2.1.0   http://radimrehurek.com/gensim/
    
2. download or clone the software in a working directory
3. in this directory, download the required data from public databases by executing 'get_data.sh'
4. execute "python mirai.py" from the working directory

### Academic citing¶
If you use this software, please cite:

    Pasquier C, Gardès J Prediction of miRNA-disease associations with a vector space model.
    Scientific Reports 2016, 6(June):27036

**BibTeX entry**

    @article{Pasquier2016,
    author = {Pasquier, Claude and Gard{\`{e}}s, Julien},
    title = {{Prediction of miRNA-disease associations with a vector space model}},
    journal = {Scientific Reports},
    number = {June},
    pages = {27036},
    volume = {6},
    year = {2016}
    publisher = {Nature Publishing Group},
    doi = {10.1038/srep27036},
    issn = {2045-2322},
    url = {http://www.nature.com/articles/srep27036},
    }
