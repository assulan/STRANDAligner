#!/usr/bin/python

# Runs STRAND in batch mode give directory of input pages.
import codecs
import gzip
import itertools
import optparse
import os
import re
import sys

from random import shuffle
from StringIO import StringIO

# Used for parsing HTML
import bs4
import lxml.html
from lxml import etree

import parsers
import strand
from collections import Counter
from segmenter import Segmenter
from py_aligner import PyGaleChurchAligner


def main():
    parser = optparse.OptionParser()
    parser.add_option("-i", "--input-dir", dest="input_dir", default="",
                        type="string", help="Location of the uncompressed mined webpages, e.g. articles folder")
    parser.add_option("-o", "--output-dir", dest="output_dir", default="strand_output",
                        type="string", help="Location where to store parallel sentences")
    parser.add_option("-c", "--chunks-output-dir", dest="chunks_dir", default="chunks_output",
                        type="string", help="Location where to store chunk lengths for parallel pages")
    parser.add_option("-t", "--chunks-tagged-dir", dest="chunks_tagged_dir", default="chunks_tagged_output",
                        type="string", help="Location where to store tagged and aligned parallel pages")
    parser.add_option("-p", "--out-df", dest="out_df", default="df_percentage",
                        help="File name where difference percentage is stored (located in current dir)")

    (opts, args) = parser.parse_args()
    if opts.input_dir == "":
        print "No input dir given"
        return

    # Create output folder if does not exist
    if not os.path.exists(opts.output_dir):
        os.mkdir(opts.output_dir)

    # Create chunks folder if does not exist
    if not os.path.exists(opts.chunks_dir):
        os.mkdir(opts.chunks_dir)

    # Create chunks tagged folder if does not exist
    if not os.path.exists(opts.chunks_tagged_dir):
        os.mkdir(opts.chunks_tagged_dir)

    # List files in output dir
    files = os.listdir(opts.input_dir)
    en_files = [f for f in files if f.endswith("en")]
    kz_files = [f for f in files if f.endswith("kz")]

    df_percentage_file = open(opts.out_df, "w")
    df_percentage_file.write("filename diff_strand diff_nu n M N alignment_cost\n")
    for i, en_file in enumerate(en_files):
        fname = en_file[:-2]
        kz_file = [f for f in kz_files if f.startswith(fname)][0]
        (correlation, p_value, difference_percentage, n, N, M, Y, X, alignment_cost) = run_strand(en_file, kz_file, opts.input_dir, opts.output_dir, opts.chunks_dir, opts.chunks_tagged_dir)
        alignment_cost = -float(alignment_cost)
        df_percentage = (alignment_cost) / (M + N)
        out_str = "%s %s %s %s %s %s %s\n" % (fname, str(difference_percentage), str(df_percentage), str(n), str(M), str(N), str(alignment_cost))
        df_percentage_file.write(out_str)
        if i % 1000 == 0:
            print "Processing file %i" % i
    df_percentage_file.close()
    print "Done."

# Return set from list of tag chunks.
def get_set(tag_chunks):
    myset = []
    for line in tag_chunks:
        if line.startswith("[") and line.endswith("]"):
            # tag
            myset.append(line.strip("\n").encode())
    return Counter(myset)


# Return the length of string s without whitespace characters.
def mylen(s):
    from string import whitespace
    return len([c for c in s if c not in whitespace])

# Run strand algo on English and Kazakh webpage pairs.
def run_strand(file_en, file_kz, input_dir, out_dir, chunks_dir, chunks_tagged_dir):
  # Initialize the HTML parser and aligners
  strand_parser = etree.HTMLParser(encoding="utf-8", target=parsers.StrandTarget())
  gc_aligner = PyGaleChurchAligner()
  strand_aligner = strand.StrandAligner()

  # We will always be working with English
  segmenter = Segmenter("English")

  webpages = [open(f, 'r').read() for f in [os.path.join(input_dir, file_en), os.path.join(input_dir, file_kz)]]
  tag_chunks = {"en": None, "kz": None}
  for i, webpage in enumerate(webpages):
    try:
      tagchunks = apply_parser(webpage, strand_parser)
      if i == 0:
        tag_chunks['en'] = tagchunks
      else:
        tag_chunks['kz'] = tagchunks
    except:
      print "Error parsing input file"
      return

  source_strand = tag_chunks["en"].split('\n')
  target_strand = tag_chunks["kz"].split('\n')
  source_set = get_set(source_strand)
  source_set_copy = source_set.copy()
  target_set = get_set(target_strand)

  source_set.subtract(target_set)
  Y = float(sum(source_set.values())) if sum(source_set.values()) > 0 else 0.0

  target_set.subtract(source_set_copy)
  X = float(sum(target_set.values())) if sum(target_set.values()) > 0 else 0.0
  # Extract parallel sentences and similarity metrics using STRAND
  (source_sents, target_sents, correlation, p_value, diff_percentage, n, N, M, alignment_cost) = \
      strand_extract_and_clean(strand_aligner, gc_aligner,
                               source_strand,  target_strand,
                               segmenter,      segmenter, os.path.join(chunks_dir, '%schunks' % file_en[:-2]),
                               os.path.join(chunks_tagged_dir, "%stagged_chunks") % file_en[:-2])

  # Write sentences to output files
  for out_file in [os.path.join(out_dir, file_en), os.path.join(out_dir, file_kz)]:
    out = codecs.open(out_file, mode="w", encoding="utf-8")
    sents = []
    if out_file.endswith("en"):
        sents = source_sents
    elif out_file.endswith("kz"):
        sents = target_sents
    for s in sents:
      out.write(s + "\n")
    out.close()
  return (correlation, p_value, diff_percentage, n, N, M, Y, X,alignment_cost)


# Run STRAND on the source and target HTML (parsed to tagchunks), and return
# a tuple (sentences_1, sentences_2, correlation, p_value, difference_percentage, n)
def strand_extract_and_clean(strand_aligner, gc_aligner, source, target, source_seg, target_seg, tagged_output_file, tagged_aligned_out_file):
  source_tagchunks = strand_aligner.create_tag_chunk_stream(source)
  target_tagchunks = strand_aligner.create_tag_chunk_stream(target)
  grid_size = len(source_tagchunks) * len(target_tagchunks)
  if grid_size > 1000000000:
    return ([], [])
  alignment = strand_aligner.align(source_tagchunks, target_tagchunks)
  correlation, p_value, diff_percentage, n, s_size, t_size, alignment_cost = strand_aligner.get_similarity_metrics(source_tagchunks, target_tagchunks)
  source_out = []
  target_out = []

  # open file for writing tagged aligned output
  f_tagged_aligned = codecs.open(tagged_aligned_out_file, mode='w', encoding='utf-8')
  f_out = codecs.open(tagged_output_file, mode="w", encoding="utf-8")
  f_out.write('length_kz,length_en\n')
  for (s, t) in alignment:
    # write tagged aligned data
    f_tagged_aligned.write("%s,%s\n" % (str(t).decode("utf-8"), str(s).decode("utf-8")))
    if (not s or not t):
        # one of the chunks is None
        if not s and t and t.tc_type == strand.TCType.CHUNK:
            # source chunk is None
            f_out.write("%s,%s\n" % (len(t.chunk_data), str(0)))
            # f_out.write("%s,%s\n" % (mylen(t.chunk_data), str(0)))
        elif not t and s and s.tc_type == strand.TCType.CHUNK:
            # target chunk is None
            f_out.write("%s,%s\n" % (str(0), len(s.chunk_data)))
            # f_out.write("%s,%s\n" % (str(0), mylen(s.chunk_data)))
        elif not s and not t:
            # both are None
            f_out.write("%s,%s\n" % (str(0), str(0)))
    if (s and s.tc_type == strand.TCType.CHUNK
        and t and t.tc_type == strand.TCType.CHUNK):
      f_out.write("%s,%s\n" % (len(t.chunk_data), len(s.chunk_data)))
      # f_out.write("%s,%s\n" % (mylen(t.chunk_data), mylen(s.chunk_data)))
      source_sents = source_seg.process(unicode(s.chunk_data))
      target_sents = target_seg.process(unicode(t.chunk_data))
      grid_size = len(source_sents) * len(target_sents)
      if grid_size > 1000000000:
        continue
      (cost, aligned_source, aligned_target) = gc_aligner.align(
          source_sents, target_sents)
      for i in xrange(0, len(aligned_source)):
        s_sent = aligned_source[i]
        t_sent = aligned_target[i]
        #if s_sent != t_sent and alpha_min_length(s_sent, t_sent) >= 5 and end_punc(s_sent, t_sent) == 1:
        if s_sent != t_sent:
          source_out.append(s_sent)
          target_out.append(t_sent)
  f_tagged_aligned.close()
  f_out.close()
  return (source_out, target_out, correlation, p_value, diff_percentage, n, s_size, t_size, alignment_cost)

# Usese BeautifulSoup to handle encodings (taken from lxml tutorial)
def decode_html(html_string):
  converted = bs4.UnicodeDammit(html_string, isHTML=True)
  if not converted.unicode_markup:
    raise UnicodeDecodeError(
      "Failed to detect encoding, tried [%s]",
      ', '.join(converted.tried_encodings))
  return converted.unicode_markup

# Passes the HTML through the given parser. Uses the BeautifulSoup parser as a
# failsafe.
def apply_parser(html, parser):
  result = ""
  try:
    result = etree.parse(StringIO(html), parser)
  except: # TODO: find the specific error
    try:
      result = etree.parse(StringIO(decode_html(html)), parser)
    except:
      soup = bs4.BeautifulSoup(html, "lxml")
      result = etree.parse(StringIO(str(soup)), parser)

  return result


if __name__ == "__main__":
    main()


