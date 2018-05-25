# To change just drag filenames between the triple quotes """
# Make sure each entry is on a separate line

files_str = """
example_data/interactBSC_FSAP.prot.xml  example_data/interactbsc_fsap.pep.xml

"""

import os
import glob
from protxml2csv import convert_to_csv

for line in files_str.splitlines():
  words = line.split()
  if not words:
    continue
  protxml = words[0]
  pepxmls = words[1:]

  print "\n\n" + "="*60
  if not os.path.isfile(protxml):
    print "Couldn't find protxml", protxml
    continue
  for pepxml in pepxmls:
    if not os.path.isfile(pepxml):
      print "Couldn't find protxml", pepxml
      continue

  params = {
    'protxml': protxml,
    'pepxmls': pepxmls,
    'csv': protxml + '.csv',
    'protein_probability_cutoff': None, # or None
    'protein_error_cutoff': 0.01, # or None
    'peptide_probability_cutoff': None, # or None
    'peptide_error_cutoff': 0.01, # or None
    'is_skip_no_unique': False,
    'is_only_one_sibling': False,
  }
  convert_to_csv(params)

