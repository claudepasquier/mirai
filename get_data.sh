DATADIR=./data

curl ftp://nlmpubs.nlm.nih.gov/online/mesh/2015/supp2015.xml -o $DATADIR/www.nlm.nih.gov/mesh/supp2015.xml --create-dirs
curl ftp://nlmpubs.nlm.nih.gov/online/mesh/2015/desc2015.xml -o $DATADIR/www.nlm.nih.gov/mesh/desc2015.xml --create-dirs

curl ftp://mirbase.org/pub/mirbase/21/aliases.txt.gz -o $DATADIR/www.mirbase.org/aliases.txt.gz --create-dirs
curl ftp://mirbase.org/pub/mirbase/21/mature.fa.gz -o $DATADIR/www.mirbase.org/mature.fa.gz --create-dirs
curl ftp://mirbase.org/pub/mirbase/21/miRNA.dead.gz -o $DATADIR/www.mirbase.org/miRNA.dead.gz --create-dirs
curl ftp://mirbase.org/pub/mirbase/21/hairpin.fa.gz -o $DATADIR/www.mirbase.org/hairpin.fa.gz --create-dirs
curl ftp://mirbase.org/pub/mirbase/21/miFam.dat.gz -o $DATADIR/www.mirbase.org/miFam.dat.gz --create-dirs
curl ftp://mirbase.org/pub/mirbase/21/miRNA.dat.gz -o $DATADIR/www.mirbase.org/miRNA.dat.gz --create-dirs
curl ftp://mirbase.org/pub/mirbase/21/genomes/hsa.gff3 -o $DATADIR/www.mirbase.org/hsa.gff3 --create-dirs

gunzip $DATADIR/www.mirbase.org/*.gz
gunzip $DATADIR/generated/pubmed/abstracts_of_mirbase.xml.gz
gunzip $DATADIR/generated/disease-similarities.json.gz
gunzip $DATADIR/generated/meshData.json.gz
gunzip $DATADIR/miRNA-assoc/hmdd_v2.txt.gz
gunzip $DATADIR/mirtarbase.mbc.nctu.edu.tw/hsa_MTI.csv.gz
