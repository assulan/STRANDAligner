# STRAND Aligner

An implementation of the HTML alignment algorithm of STRAND:

P. Resnik and N. A Smith. The web as a parallel corpus. Computational
Linguistics, 29(3):349-380, 2003.

This also contains an implementation of the Gale and Church sentence alignment
algorithm:

W. A. Gale and K. W. Church. A program for aligning sentences in
bilingual corpora. Computational Linguistics, 19:75-102, March 1993.


# Installation prerequisites
NOTE: You might need sudo permissions for some of the steps below.

* pip install numpy
* pip install cython
* pip install beautifulsoup4
* sudo apt-get install libxml2-dev libxslt-dev python-dev
* pip install lxml
* apt-get install python-scipy
* pip install nltk
* open python shell
* import nltk
* nltk.download()
* type in 'd' for download, and enter 'punkt'
* checkout maxent from https://github.com/lzhang10/maxent.git
* cd into maxent folder
* follow the installation instructions here https://github.com/lzhang10/maxent/blob/master/INSTALL.
   This will install c++ maxent lib.
* cd into python dir within maxent, and follow installation instructions in README file.

# Installation
1. Checkout the project from git
2. cd into dir called 'strand'
3. python setup.py install

# Running the program
1. cd into dir called 'data' under STRANDAligner. There is a folder named 'articles' that contains 360 parallel pages named 1.{en,kz}, etc.
2. sudo python ../strand/run_strand_batch -i articles
3. when the above commands completes the following will be created in 'data' directory:
	* chunks_output directory - contains files 1.chunks, 2.chunks, etc. Each file contains lengths of aligned chunks from parallel texts.
	* chunks_tagged_output directory - contains files 1.tagged_chunks, etc. Each file contains aligned chunks together with tags.
	* df_percentage file - contains data needed for learning algorithm used in our paper.
	* strand_output directory - contains files 1.{en, kz}, etc. with chunks.

# TODO
* Remove unused code from run_strand_batch.py
* Do refactoring to make code more flexible, remove hard coded stuff.

	
