 # -*- coding: utf-8 -*-
from __future__ import print_function
from pprint import pprint
import xml.etree.ElementTree as etree
import json
#import logging
import re
import os
import ntpath, posixpath, macpath

#import parse
#import proteins as parse_proteins

"""
Parsers to read TPP output files: PEPXML and PROTXML.
"""

#logger = logging.getLogger('tpp')


def parse_scan(scan_elem, nsmap):
  scan = parse_attrib(scan_elem)
  scan['matches'] = []
  tag = lambda tag_id: fixtag('', tag_id, nsmap)
  for search_elem in scan_elem.findall(fixtag('', "search_result", nsmap)):
    search_hit_elem = search_elem[0] 
    match = parse_attrib(search_hit_elem)
    match['modified_sequence'] = match['peptide']

    match['other_seqids'] = []
    for alt_protein in search_hit_elem.findall(fixtag('', 'alternative_protein', nsmap)):
      match['other_seqids'].append(alt_protein.attrib['protein'])

    match['modifications'] = []
    for modified_elem in search_hit_elem.findall(fixtag('', 'modification_info', nsmap)):
      attr = parse_attrib(modified_elem)
      match['modified_sequence'] = attr['modified_peptide']
      for modification_elem in modified_elem.findall(fixtag('', 'mod_aminoacid_mass', nsmap)):
        attr = parse_attrib(modification_elem)
        #attr['i'] = attr['position'] - 1
        attr['i'] = attr['position']
        del attr['position']
        match['modifications'].append(attr)
    
    n_term = ''
    n_term_mod = None
    mp_list = filter(None, re.split("[\[\]]", match['modified_sequence']))
    
    if mp_list[0] == 'n':
      n_term_mod = {'i':0, 'mass':mp_list[1]}
      n_term = mp_list[0] + '[' + mp_list[1] + ']'
      mp_list = mp_list[2:]
    
    idx = 0
    idx_list = [x['i'] for x in match['modifications']]
    for i in range(0, (len(mp_list) - 1), 2):
        idx += len(mp_list[i])
        if idx not in idx_list:
            match['modifications'].append({'i':idx, 'mass':mp_list[i+1]})
    
    mp = match['peptide']
    for mod in sorted(match['modifications'],
                key=lambda m: m['i'],
                reverse=True):

        #print (mod['i'] in range(0, 1+len(match['peptide'])))
        if mod['i'] in range(0, len(match['peptide'])):
          #p = mod['i'] + 1
          p = mod['i']
            
          mp = mp[:p] + '[{}]'.format(int(mod['mass'])) + mp[p:]
    match['modified_peptide'] = n_term + mp
    if n_term_mod != None:
      match['modifications'].append(n_term_mod)

    for score_elem in search_hit_elem.findall(tag('search_score')):
      match.update(parse_name_value(score_elem))
    
    #for analysis_elem in search_hit_elem.find(fixtag('', 'analysis_result', nsmap)):
    #  if analysis_elem.tag == fixtag('', 'peptideprophet_result', nsmap): #can change to interprophet here
    #    match.update(parse_attrib(analysis_elem))
    #    for param_elem in analysis_elem[0]:
    #      match.update(parse_name_value(param_elem))
    
    for analysis_elem in search_hit_elem.findall(fixtag('', 'analysis_result', nsmap)):
      elem = analysis_elem.find(fixtag('', 'interprophet_result', nsmap)) 
      if elem != None:
        if elem.tag == fixtag('', 'interprophet_result', nsmap):
          match.update(parse_attrib(elem))
      elem_peptideprophet = analysis_elem.find(fixtag('', 'peptideprophet_result', nsmap))
      if elem_peptideprophet != None:
        match['peptideprophet_probability'] = elem_peptideprophet.attrib['probability']
    scan['matches'].append(match)
  return scan


def parse_peptide_probabilities(elem, nsmap):
  # try with error_point
  error_points = elem.findall(fixtag('', 'error_point', nsmap))
  if len(error_points) == 0:
    charge = 0
    for charge_elem in elem.findall(fixtag('', 'roc_error_data', nsmap)):
      if charge_elem.attrib['charge'] == 'all':
        error_points = charge_elem.findall(fixtag('', 'error_point', nsmap))
        break
  probs = []
  for elem in error_points:
      attrib = parse_attrib(elem)
      probs.append({
        'error': attrib['error'],
        'prob': attrib['min_prob'],
      })
  probs.sort(key=lambda d:d['error'])
  return probs


class PepxmlReader(object):
  def __init__(self, pepxml, fdr_method, unique, prob_cutoff=None, error_cutoff=None):
    self.pepxml = pepxml
    self.probs = None
    self.source_names = []
    self.i_source = None
    self.unique = unique
    self.filter_method = fdr_method
    self.prob_cutoff = float(prob_cutoff)
    self.error_cutoff = float(error_cutoff)
    self.nsmap = {}
    
  def iter(self):
    
    for event, elem in etree.iterparse(self.pepxml, events=('start', 'end', 'start-ns')):
      if event == 'start-ns':
        self.nsmap.update({elem})
      elif event == 'start':
        if elem.tag == fixtag('', 'msms_run_summary', self.nsmap):
          fname = elem.attrib['base_name']
          self.source_names.append(fname)
          self.i_source = len(self.source_names) - 1
      elif event == 'end':
        if elem.tag == fixtag('', 'spectrum_query', self.nsmap):
          scan = parse_scan(elem, self.nsmap)
          if self.i_source is not None:
            scan['source'] = self.source_names[self.i_source]
          if self.prob_cutoff is None:
            yield scan
          else:
            if len(scan['matches']) == 1:
              #print ('probability: ' + str(scan['matches'][0]['probability']))
              #print (scan['matches'][0]['probability'] >= self.prob_cutoff)

              if scan['matches'][0]['probability'] >= self.prob_cutoff:
                yield scan
            else:
              remove = []
              for i in range(len(scan['matches'])):
                if scan['matches'][i]['probability'] <= self.prob_cutoff:
                  remove.append(i)
                if len(remove) == 0:
                  yield scan
                else:
                  for i in remove:
                    del(scan['matches'][i])
                  if len(scan['matches']) > 0:
                    yield scan
          elem.clear()
        
        #elif elem.tag == fixtag('', 'peptideprophet_summary', self.nsmap):
        #  print ('found peptideprophet tag')
        elif elem.tag == fixtag('', 'interprophet_summary', self.nsmap):
          self.probs = parse_peptide_probabilities(elem, self.nsmap)
          if (self.filter_method == 'iFDR'):
            print("[iFDR] Determining Spectrum Probability Cutoff without using Mayu...")
            self.prob_cutoff = error_to_probability(self.probs, self.error_cutoff)
          elif (self.filter_method == 'noFDR'):
            print('[noFDR] Spectrum Probability Cutoff is used directly...')
          else: 
            if self.filter_method == 'mayu':
              print("[mayu] Using Mayu to determine Peptide Probability cutoff...")
            else:
              print("[mayu_iFDR] Using Mayu to determine Peptide Probability cutoff...")
          print ("Spectrum Probability Cutoff: " + str(self.prob_cutoff))
          elem.clear()


def parse_protein_probabilities(elem, nsmap):
  probs = []
  for data_point in elem.findall(fixtag('', 'protein_summary_data_filter', nsmap)):
    attrib = parse_attrib(data_point)
    probs.append({
      'error': attrib['false_positive_error_rate'],
      'prob': attrib['min_probability'],
    })
  probs.sort(key=lambda d:d['error'])
  return probs


def parse_protein_group(elem, nsmap):
  group = parse_attrib(elem)
  group['proteins'] = []
  for protein_elem in elem.findall(fixtag('', 'protein', nsmap)):
    protein = parse_attrib(protein_elem)
    protein['group_number'] = group['group_number']
    
    length_elem = protein_elem.find(fixtag('', 'parameter', nsmap))
    if length_elem is not None:
      protein['prot_length'] = length_elem.attrib['value']

    annotation_elem = protein_elem.find(fixtag('', 'annotation', nsmap))
    if annotation_elem is not None:
      protein['description'] = annotation_elem.attrib['protein_description']

    protein['other_seqids'] = []
    for alt_protein in protein_elem.findall(fixtag('', 'indistinguishable_protein', nsmap)):
      protein['other_seqids'].append(alt_protein.attrib['protein_name'])

    protein['other_seqids'] = protein['other_seqids']
    protein['protein_name'] = protein['protein_name']

    protein['peptides'] = []
    n_unique_peptide = 0
    for peptide_elem in protein_elem.findall(fixtag('', 'peptide', nsmap)):
      peptide = parse_attrib(peptide_elem)
      protein['peptides'].append(peptide)
      peptide['modifications'] = []
      peptide['modified_sequence'] = peptide['peptide_sequence']
      for modified_elem in peptide_elem.findall(fixtag('', 'modification_info', nsmap)):
        attr = parse_attrib(modified_elem)
        peptide['modified_sequence'] = attr['modified_peptide']
        for modification_elem in modified_elem.findall(fixtag('', 'mod_aminoacid_mass', nsmap)):
          attr = parse_attrib(modification_elem)
          peptide['modifications'].append(attr)

    group['proteins'].append(protein)
  return group


def read_protxml(protxml):
  nsmap = {}
  distribution = {}
  protein_groups = {}
  for event, elem in etree.iterparse(protxml, events=('end', 'start-ns')):
    if event == 'start-ns':
      nsmap.update({elem})
    if event == 'end':
      if elem.tag == fixtag('', 'protein_group', nsmap):
        group = parse_protein_group(elem, nsmap)
        protein_groups[group['group_number']] = group
        elem.clear()
      elif elem.tag == fixtag('', 'proteinprophet_details', nsmap):
        distribution = parse_protein_probabilities(elem, nsmap)
        elem.clear()
  return protein_groups, distribution


def error_to_probability(distribution, error):
  """
  Given a False-Positive-Error vs. Probability distribution,
  Cacluates the probability for an acceptable FPE.
  """
  fractionate = lambda a0, a1, a: (a-a0)/(a1-a0)
  interpolate = lambda a0, a1, f: a0 + f*(a1-a0)
  n = len(distribution)
  for i in range(1, n):
    error0 = distribution[i-1]['error']
    error1 = distribution[i]['error']
    #print ('error0: ' + str(error0))
    #print ('error: ' + str(error))
    #print ('error1: ' + str(error1))
    if error0 <= error < error1:
      prob0 = distribution[i-1]['prob']
      prob1 = distribution[i]['prob']
      f = fractionate(error0, error1, error)
      return interpolate(prob0, prob1, f)
  return None


"""
Utility parsing functions that is useful for
different proteomics data parsers.
"""

# string parsers

float_regex_pattern = r"""
^ 
[-+]? # optional sign
(?:
    (?: \d* \. \d+ ) # .1 .12 .123 etc 9.1 etc 98.1 etc
    |
    (?: \d+ \.? ) # 1. 12. 123. etc 1 12 123 etc
)
# followed by optional exponent part if desired
(?: [Ee] [+-]? \d+ ) ?
$
"""
float_regex = re.compile(float_regex_pattern, re.VERBOSE)
    
    
def basename(fname):
  for m in ntpath, posixpath, macpath:
    if m.sep in fname:
      return m.basename(fname)
  return fname


def parse_string(s):
  if re.search(r'^[-+]?\d+$', s):
    return int(s)
  elif float_regex.match(s):
    return float(s)
  else:
    return s


# tsv helper function


# def split_tab(line):
#   if line[-1] == '\n':
#     line = line[:-1]
#   return line.split('\t')


# XML helper functions


def fixtag(ns, tag, nsmap):
  return '{' + nsmap[ns] + '}' + tag


def parse_attrib(elem):
  result = {}
  for key, value in elem.attrib.items():
    result[key] = parse_string(value)
  return result


def parse_name_value(elem):
  attrib = parse_attrib(elem)
  return { attrib['name']: attrib['value'] }


def save_data_dict(data, fname):
  with open(fname, 'w') as f:
    pprint(data, stream=f, indent=2)
