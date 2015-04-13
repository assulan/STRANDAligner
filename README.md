# STRAND Aligner

An implementation of the HTML alignment algorithm of STRAND:

P. Resnik and N. A Smith. The web as a parallel corpus. Computational
Linguistics, 29(3):349-380, 2003.

This also contains an implementation of the Gale and Church sentence alignment
algorithm:

W. A. Gale and K. W. Church. A program for aligning sentences in
bilingual corpora. Computational Linguistics, 19:75-102, March 1993.

NOTE: You might need sudo permissions for some of the below steps.

# Installation prerequisites
* pip install numpy
* pip install cython
* pip install beautifulsoup4
* sudo apt-get install libxml2-dev libxslt-dev python-dev
* pip install lxml
* apt-get install python-scipy
* pip install nltk
* checkout maxent from https://github.com/lzhang10/maxent.git
* cd into maxent folder
* follow the installation instructions here https://github.com/lzhang10/maxent/blob/master/INSTALL.
   This will install c++ maxent lib.
* cd into python dir within maxent, and follow installation instructions in README file.

# Installation
1. Checkout the project from git
2. cd into dir called 'strand'
3. python setup.py install

