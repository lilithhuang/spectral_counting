# -*- coding: utf-8 -*-
# (c) 2013 Bosco Ho
# from __future__ import print_function
import os
import sys
import getopt
import csv
import json
from pprint import pprint
import re
import urllib2
from collections import Counter

#import uniprot
import peptagram.tpp
#import peptagram.parse

import time

start = time.time()

usage = """
Reads a TPP Protein XML file and outputs it into a convenient CSV format.
For better protein name description, interrogates the Uniprot website to
extract more useful meta-data.

Usage: %s protein_xml [output_csv]
""" % os.path.basename(__file__)


def convert_to_csv(params):
    protxml_file, pepxml_file, mayu_file, outfile, peptide_quant_csv, peptide_csv, protein_quant_csv, protein_report_csv, spectrum_csv, prefix, protein_quant_count, protein_probability_cutoff, protein_error_cutoff, peptide_probability_cutoff, peptide_error_cutoff, spectrum_error_cutoff, fdr_method, cluster, unique = get_opt(sys.argv[1:])
    
    pepxmls = [pepxml_file]
    proteins_csv = outfile #remove later
    if not proteins_csv:
        proteins_csv = protxml_file + '.csv'
    is_skip_no_unique = params['is_skip_no_unique']  
   
    if (fdr_method == 'mayu' or fdr_method == 'mayu_iFDR') and mayu_file:
        with open(mayu_file, 'rb') as csvfile:
            reader = csv.DictReader(csvfile)
            
            if protein_error_cutoff == None:
                protein_error_cutoff = 1 # set to 1
            if spectrum_error_cutoff == None:
                spectrum_error_cutoff = 1 # set to 1
            for row in reader:

                if float(row['protFDR']) > protein_error_cutoff or float(row['pepFDR']) > peptide_error_cutoff or float(row['mFDR']) > spectrum_error_cutoff:
                    break
                else:
                    peptide_probability_cutoff = float(row['IP/PPs'])

    print "LOADING:", protxml_file
    protxml_groups, protein_probs = peptagram.tpp.read_protxml(protxml_file)

    # create an auxillary mapping to access proteins more easily when loading pepxml files
    protein_by_seqid = {}
    for group_id, protxml_group in protxml_groups.items():
        for protein in protxml_group['proteins']:
            seqids = [protein['protein_name']] + protein['other_seqids']
            for seqid in seqids:
                protein_by_seqid[seqid] = protein
    
    ####################
    benchmark1 = time.time()
    
    print ('finished parsing protxml')
    print (benchmark1 - start)
    ####################
    
    # load pepxml, run through scans, and match by seqids, then loads scans
    # into the peptide['scans'] structure
    for group_id, protxml_group in protxml_groups.items():
        for protein in protxml_group['proteins']:
            for peptide in protein['peptides']:
                peptide['scans'] = []
    for pepxml in pepxmls:
        print "MATCHING:", pepxml
        pepxml_reader = peptagram.tpp.PepxmlReader(
            pepxml, fdr_method, unique, error_cutoff=spectrum_error_cutoff, prob_cutoff=peptide_probability_cutoff)
        
        #print ('from reader')
        #print (pepxml_reader.prob_cutoff)
        #pepxml_reader = peptagram.tpp.PepxmlReader(
            #pepxml, error_cutoff=peptide_error_cutoff, prob_cutoff=None)
        for scan in pepxml_reader.iter():

            for match in scan['matches']:

                seqids = [match['protein']] + match['other_seqids']
                for seqid in seqids:
                    if seqid not in protein_by_seqid:
                        print('Warning {} not found in prot.XML'.format(seqid))
                        continue
                    protein = protein_by_seqid[seqid]
                    test_charge = scan['assumed_charge']
                    #test_sequence = match['modified_sequence']
                    test_sequence = match['peptide']
                    
                    for peptide in protein['peptides']:
                        charge = peptide['charge']
                        sequence = peptide['peptide_sequence']
                        
                        #if test_charge == charge and test_sequence == sequence:
                        if test_sequence == sequence:

                            if scan not in peptide['scans']:
                                peptide['scans'].append(scan)
                                break
    
    ####################
    benchmark2 = time.time()
    
    print ('finished parsing pepxml')
    print (benchmark2 - benchmark1)
    ####################
    
    n_unique_scan_total = 0
    reps = []
    for group_id, protxml_group in protxml_groups.items():
        protxml_group['cluster'] = {'group_number': group_id, 'protein_list': [], 'protein_group_id': '', 'description': '', 'counts': {}}
    
        for protein in protxml_group['proteins']:
            if protein['protein_name'] not in protxml_group['cluster']['protein_list']:
                protxml_group['cluster']['protein_list'].append(protein['protein_name'])
            if protein['group_sibling_id'] == 'a':
                protxml_group['cluster']['protein_group_id'] = protein['protein_name']
                protxml_group['cluster']['description'] = 'Cluster of ' + protein['description']

            protein['counts'] = {}
            #protein['n_peptide'] = len(protein['peptides']) # this can be wrong -> need to check length of scans (might not need this)
            protein['n_peptide'] = 0
            
            for i in reversed(range(len(protein['peptides']))):
                peptide = protein['peptides'][i]
      
      
                #for peptide in protein['peptides']: 
                #if float(peptide['weight']) > 0.50:
                #  protein['n_unique_peptide'] += 1 # should check shared here
                #  if 'scans' in peptide:
                #    protein['n_unique_scan'] += len(peptide['scans']) # should check charge and modification here

                
                if unique == True:
                    if peptide['is_nondegenerate_evidence'] == 'N':
                        del protein['peptides'][i]
                        continue
                if len(peptide['scans']) == 0: #handles case where peptide has no scans
                    continue
                else:
                    protein['n_peptide']+=1
                    for scan in peptide['scans']:
                        rep = scan['spectrum'].split('.',1)[0]
                        if rep not in reps: #list of reps used to generate output files
                            reps.append(rep)
                        
                        if rep not in protein['counts']:
                            protein['counts'][rep] = {'unique_peptide': [], 'unique_spectrum': [], 'total_spectrum': []}
                        #protein total spectrum count
                        if scan['start_scan'] not in protein['counts'][rep]['total_spectrum']:
                            protein['counts'][rep]['total_spectrum'].append(scan['start_scan'])
                        if peptide['is_nondegenerate_evidence'] == 'Y':
                            #protein exclusive unique peptide count
                            if peptide['peptide_sequence'] not in protein['counts'][rep]['unique_peptide']:
                                protein['counts'][rep]['unique_peptide'].append(peptide['peptide_sequence'])
                            #protein exclusive unique spectrum count
                            if (peptide['peptide_sequence'], peptide['modifications'], peptide['charge']) not in protein['counts'][rep]['unique_spectrum']:
                                protein['counts'][rep]['unique_spectrum'].append((peptide['peptide_sequence'], peptide['modifications'], peptide['charge']))
                        
                        if rep not in protxml_group['cluster']['counts']:
                            protxml_group['cluster']['counts'][rep] = {'unique_peptide': [], 'unique_spectrum': [], 'total_spectrum': []}
                        #cluster total spectrum count
                        if scan['start_scan'] not in protxml_group['cluster']['counts'][rep]['total_spectrum']:
                            protxml_group['cluster']['counts'][rep]['total_spectrum'].append(scan['start_scan'])
                        #cluster exclusive unique peptide count
                        if peptide['peptide_sequence'] not in protxml_group['cluster']['counts'][rep]['unique_peptide']:
                            protxml_group['cluster']['counts'][rep]['unique_peptide'].append(peptide['peptide_sequence'])
                        #cluster exclusive unique spectrum count
                        if (peptide['peptide_sequence'], peptide['modifications'], peptide['charge']) not in protxml_group['cluster']['counts'][rep]['unique_spectrum']:
                            protxml_group['cluster']['counts'][rep]['unique_spectrum'].append((peptide['peptide_sequence'], peptide['modifications'], peptide['charge']))
                        
                        for rep in protein['counts']:
                            n_unique_scan_total += len(protein['counts'][rep]['total_spectrum'])

    ####################
    benchmark3 = time.time()
    
    print ('finished counting spectrums')
    print (benchmark3 - benchmark2)
    ####################

    #for group_id, protxml_group in protxml_groups.items():
    #    for protein in protxml_group['proteins']:
    #        percent = 100.0*len(protein['total_spectrum'])/n_unique_scan_total
    #        percent = float('%.2f' % percent)
    #        protein['percent_spectra'] = percent
    
    if fdr_method == 'noFDR':
        print('[noFDR] Protein Probability Cutoff is used directly...')
        prob_cutoff = float(protein_probability_cutoff)
    elif fdr_method == 'mayu':
        print('[mayu] No Protein Probability Cutoff used...')
        prob_cutoff = 0.0
    else: #fdr_method is iFDR or mayu_iFDR
        if fdr_method == 'iFDR':
            print("[iFDR] Determining Protein Probability Cutoff from Protein FDR...")
        else:
            print("[mayu_iFDR] Determining Protein Probability Cutoff from Protein FDR...")
        prob_cutoff = peptagram.tpp.error_to_probability(protein_probs, protein_error_cutoff)
    if fdr_method != 'mayu':
        print('Protein Probability Cutoff: ' + str(prob_cutoff))
        
#    for group_id in protxml_groups.keys():
#        #if group_id == 20:
#        #    print ('checking')
#        #    print (group_id)
#
#        protxml_group = protxml_groups[group_id]
#        proteins = protxml_group['proteins']
#        for i in reversed(range(len(proteins))):
#            protein = proteins[i]
#            #print ('protein probability: ' + str(protein['probability']))
#            #print (protein['probability'] >= prob_cutoff)
#            if protein['probability'] < prob_cutoff or protein['protein_name'].startswith('DECOY'):
#                del proteins[i]
#                protxml_group['cluster']['protein_list'].remove(protein['protein_name'])
#            elif 'n_peptide' not in protein or protein['n_peptide'] == 0:
#                del proteins[i]
#                protxml_group['cluster']['protein_list'].remove(protein['protein_name'])
#            else:
#                remove = True
#                for rep in protein['counts']:
#                    if len(protein['counts'][rep]['unique_peptide']) != 0:
#                        remove = False
#                if remove == True:
#                    del proteins[i]
#                    protxml_group['cluster']['protein_list'].remove(protein['protein_name'])
#
#            protxml_group['cluster']['protein_name'] = ','.join(protxml_group['cluster']['protein_list'])
#            
#            #TODO
#            #checking here
#            if (len(protxml_group['cluster']['protein_list']) != len(proteins)):
#                print ('----')
#                print ('error')
#                print protxml_group['cluster']['protein_list']
#                print [protein['protein_name'] for protein in proteins]
#                print ('----')
#        if len(proteins) == 0:
#            del protxml_groups[group_id]
     
    peptide_dict = {}
    summary_dict = {}
    
    for rep in reps:
        summary_dict[rep] = {'groups': [], 'proteins': [], 'peptides': [], 'modified': [], 'spectra': []}
    
    n_unique_scan_total = 0
    for group_id, protxml_group in protxml_groups.items():



        #filtering proteins
        proteins = protxml_group['proteins']
        for i in reversed(range(len(proteins))):
            protein = proteins[i]
            if protein['probability'] < prob_cutoff or protein['protein_name'].startswith('DECOY'):
                del proteins[i]
                protxml_group['cluster']['protein_list'].remove(protein['protein_name'])
            elif 'n_peptide' not in protein or protein['n_peptide'] == 0:
                del proteins[i]
                protxml_group['cluster']['protein_list'].remove(protein['protein_name'])
            else:
                remove = True
                for rep in protein['counts']:
                    if len(protein['counts'][rep]['unique_peptide']) != 0:
                        remove = False
                if remove == True:
                    del proteins[i]
                    protxml_group['cluster']['protein_list'].remove(protein['protein_name'])
                
                #polulating peptide_dict
                else: 
                    for peptide in proteins[i]['peptides']: 
                        if len(peptide['scans']) == 0: #handles case where peptide has no scans
                            continue
                        else:
                            #print ('adding scans here')
                            if peptide['peptide_sequence'] not in peptide_dict:
                                peptide_dict[peptide['peptide_sequence']] = {'protein_list': [], 'counts': {}}
                                #peptide_dict[peptide['peptide_sequence']] = {'counts': {}}
                                peptide_dict[peptide['peptide_sequence']]['group_id'] = protxml_group['cluster']['protein_group_id']
                                #peptide_dict[peptide['peptide_sequence']]['protein_name'] = protxml_group['cluster']['protein_name']
                                #peptide_dict[peptide['peptide_sequence']]['protein_list'] = []
                                peptide_dict[peptide['peptide_sequence']]['group_name'] = protxml_group['cluster']['description']
                            peptide_entry = peptide_dict[peptide['peptide_sequence']]
                            if protein['protein_name'] not in peptide_entry['protein_list']:
                                peptide_entry['protein_list'].append(protein['protein_name'])
                            for scan in peptide['scans']:
                                
                                rep = scan['spectrum'].split('.',1)[0]
                                
                                #populating summary dict
                                if protein['group_number'] not in summary_dict[rep]['groups']:
                                    summary_dict[rep]['groups'].append(protein['group_number'])
                                if protein['protein_name'] not in summary_dict[rep]['proteins']:
                                    summary_dict[rep]['proteins'].append(protein['protein_name'])
                                if peptide['peptide_sequence'] not in summary_dict[rep]['peptides']:
                                    summary_dict[rep]['peptides'].append(peptide['peptide_sequence'])
                                if (peptide['peptide_sequence'], peptide['modifications']) not in summary_dict[rep]['modified']:
                                    summary_dict[rep]['modified'].append((peptide['peptide_sequence'], peptide['modifications']))
                                if scan['start_scan'] not in summary_dict[rep]['spectra']:
                                    summary_dict[rep]['spectra'].append(scan['start_scan'])
        
                                if rep not in peptide_entry['counts']:
        
                                    peptide_entry['counts'][rep] = {'scans': [], 'modified_peptides': [], 'prev_aa': scan['matches'][0]['peptide_prev_aa'], 'next_aa': scan['matches'][0]['peptide_next_aa'], 'num_tot_proteins': scan['matches'][0]['num_tot_proteins'], 'num_missed_cleavages': scan['matches'][0]['num_missed_cleavages'], 'best_prob': scan['matches'][0]['probability'], 'is_nondegenerate_evidence': peptide['is_nondegenerate_evidence'], 'n_enzymatic_termini': peptide['n_enzymatic_termini']}
                                    
                                    #scan['assumed_charge'],
                                    peptide_entry['counts'][rep]['scans'].append((scan['start_scan'], scan['matches'][0]['modified_peptide'], scan['assumed_charge'], scan['matches'][0]['peptide_prev_aa'], scan['matches'][0]['peptide_next_aa'], scan['precursor_neutral_mass'], scan['retention_time_sec'], scan['matches'][0]['num_tot_proteins'], scan['matches'][0]['num_matched_ions'], scan['matches'][0]['tot_num_ions'], scan['matches'][0]['num_tol_term'], scan['matches'][0]['calc_neutral_pep_mass'], scan['matches'][0]['num_missed_cleavages'], scan['matches'][0]['massdiff'], scan['matches'][0]['expect'], scan['matches'][0]['deltacn'] if 'deltacn' in scan['matches'][0] else "-", scan['matches'][0]['spscore'] if 'spscore' in scan['matches'][0] else "-", scan['matches'][0]['xcorr'] if 'xcorr' in scan['matches'][0] else "-", scan['matches'][0]['num_matched_peptides'] if 'num_matched_peptides' in scan['matches'][0] else "-", scan['matches'][0]['probability'], scan['matches'][0]['peptideprophet_probability']))
                                    peptide_entry['mods'] = [scan['matches'][0]['modified_peptide']] #list of unique modifications
                                    peptide_entry['counts'][rep]['modified_peptides'].append(scan['matches'][0]['modified_peptide']) #list of modified sequences
        
                                else:
                                    if scan['matches'][0]['modified_peptide'] not in  peptide_entry['mods']:
                                        peptide_entry['mods'].append(scan['matches'][0]['modified_peptide'])
                                        
                                    if (scan['start_scan'], scan['matches'][0]['modified_peptide'], scan['assumed_charge'], scan['matches'][0]['peptide_prev_aa'], scan['matches'][0]['peptide_next_aa'], scan['precursor_neutral_mass'], scan['retention_time_sec'], scan['matches'][0]['num_tot_proteins'], scan['matches'][0]['num_matched_ions'], scan['matches'][0]['tot_num_ions'], scan['matches'][0]['num_tol_term'], scan['matches'][0]['calc_neutral_pep_mass'], scan['matches'][0]['num_missed_cleavages'], scan['matches'][0]['massdiff'], scan['matches'][0]['expect'], scan['matches'][0]['deltacn'] if 'deltacn' in scan['matches'][0] else "-", scan['matches'][0]['spscore'] if 'spscore' in scan['matches'][0] else "-", scan['matches'][0]['xcorr'] if 'xcorr' in scan['matches'][0] else "-", scan['matches'][0]['num_matched_peptides'] if 'num_matched_peptides' in scan['matches'][0] else "-", scan['matches'][0]['probability'], scan['matches'][0]['peptideprophet_probability']) not in peptide_entry['counts'][rep]['scans']:
                                        peptide_entry['counts'][rep]['scans'].append((scan['start_scan'], scan['matches'][0]['modified_peptide'], scan['assumed_charge'], scan['matches'][0]['peptide_prev_aa'], scan['matches'][0]['peptide_next_aa'], scan['precursor_neutral_mass'], scan['retention_time_sec'], scan['matches'][0]['num_tot_proteins'], scan['matches'][0]['num_matched_ions'], scan['matches'][0]['tot_num_ions'], scan['matches'][0]['num_tol_term'], scan['matches'][0]['calc_neutral_pep_mass'], scan['matches'][0]['num_missed_cleavages'], scan['matches'][0]['massdiff'], scan['matches'][0]['expect'], scan['matches'][0]['deltacn'] if 'deltacn' in scan['matches'][0] else "-", scan['matches'][0]['spscore'] if 'spscore' in scan['matches'][0] else "-", scan['matches'][0]['xcorr'] if 'xcorr' in scan['matches'][0] else "-", scan['matches'][0]['num_matched_peptides'] if 'num_matched_peptides' in scan['matches'][0] else "-", scan['matches'][0]['probability'], scan['matches'][0]['peptideprophet_probability']))
                                        
                                        peptide_entry['counts'][rep]['modified_peptides'].append(scan['matches'][0]['modified_peptide'])
        
                                        if peptide_entry['counts'][rep]['best_prob'] < scan['matches'][0]['probability']:
                                            peptide_entry['counts'][rep]['best_prob'] = scan['matches'][0]['probability']
                            
                    protxml_group['cluster']['protein_name'] = ','.join(protxml_group['cluster']['protein_list'])
            
            #TODO
            #checking here
            if (len(protxml_group['cluster']['protein_list']) != len(proteins)):
                print ('----')
                print ('error')
                print protxml_group['cluster']['protein_list']
                print [protein['protein_name'] for protein in proteins]
                print ('----')
        if len(proteins) == 0:
            del protxml_groups[group_id]

    ####################
    benchmark4 = time.time()
    
    print ('finished peptide dict')
    print (benchmark4 - benchmark3)
    ####################
#------        

#            for peptide in protein['peptides']: 
#                if len(peptide['scans']) == 0: #handles case where peptide has no scans
#                    continue
#                else:
#                    print ('adding scans here')
#                    if peptide['peptide_sequence'] not in peptide_dict:
#                        peptide_dict[peptide['peptide_sequence']] = {'protein_list': [], 'counts': {}}
#                        #peptide_dict[peptide['peptide_sequence']] = {'counts': {}}
#                        peptide_dict[peptide['peptide_sequence']]['group_id'] = protxml_group['cluster']['protein_group_id']
#                        #peptide_dict[peptide['peptide_sequence']]['protein_name'] = protxml_group['cluster']['protein_name']
#                        #peptide_dict[peptide['peptide_sequence']]['protein_list'] = []
#                        peptide_dict[peptide['peptide_sequence']]['group_name'] = protxml_group['cluster']['description']
#                    peptide_entry = peptide_dict[peptide['peptide_sequence']]
#                    if protein['protein_name'] not in peptide_entry['protein_list']:
#                        peptide_entry['protein_list'].append(protein['protein_name'])
#                    for scan in peptide['scans']:
#                        
#                        rep = scan['spectrum'].split('.',1)[0]
#                        
#                        #summary dict
#                        if protein['group_number'] not in summary_dict[rep]['groups']:
#                            summary_dict[rep]['groups'].append(protein['group_number'])
#                        if protein['protein_name'] not in summary_dict[rep]['proteins']:
#                            summary_dict[rep]['proteins'].append(protein['protein_name'])
#                        if peptide['peptide_sequence'] not in summary_dict[rep]['peptides']:
#                            summary_dict[rep]['peptides'].append(peptide['peptide_sequence'])
#                        if (peptide['peptide_sequence'], peptide['modifications']) not in summary_dict[rep]['modified']:
#                            summary_dict[rep]['modified'].append((peptide['peptide_sequence'], peptide['modifications']))
#                        if scan['start_scan'] not in summary_dict[rep]['spectra']:
#                            summary_dict[rep]['spectra'].append(scan['start_scan'])
#
#                        if rep not in peptide_entry['counts']:
#
#                            peptide_entry['counts'][rep] = {'scans': [], 'modified_peptides': [], 'prev_aa': scan['matches'][0]['peptide_prev_aa'], 'next_aa': scan['matches'][0]['peptide_next_aa'], 'num_tot_proteins': scan['matches'][0]['num_tot_proteins'], 'num_missed_cleavages': scan['matches'][0]['num_missed_cleavages'], 'best_prob': scan['matches'][0]['probability'], 'is_nondegenerate_evidence': peptide['is_nondegenerate_evidence'], 'n_enzymatic_termini': peptide['n_enzymatic_termini']}
#                            
#                            #scan['assumed_charge'],
#                            peptide_entry['counts'][rep]['scans'].append((scan['start_scan'], scan['matches'][0]['modified_peptide'], scan['assumed_charge'], scan['matches'][0]['peptide_prev_aa'], scan['matches'][0]['peptide_next_aa'], scan['precursor_neutral_mass'], scan['retention_time_sec'], scan['matches'][0]['num_tot_proteins'], scan['matches'][0]['num_matched_ions'], scan['matches'][0]['tot_num_ions'], scan['matches'][0]['num_tol_term'], scan['matches'][0]['calc_neutral_pep_mass'], scan['matches'][0]['num_missed_cleavages'], scan['matches'][0]['massdiff'], scan['matches'][0]['expect'], scan['matches'][0]['deltacn'] if 'deltacn' in scan['matches'][0] else "-", scan['matches'][0]['spscore'] if 'spscore' in scan['matches'][0] else "-", scan['matches'][0]['xcorr'] if 'xcorr' in scan['matches'][0] else "-", scan['matches'][0]['num_matched_peptides'] if 'num_matched_peptides' in scan['matches'][0] else "-", scan['matches'][0]['probability'], scan['matches'][0]['peptideprophet_probability']))
#                            peptide_entry['mods'] = [scan['matches'][0]['modified_peptide']] #list of unique modifications
#                            peptide_entry['counts'][rep]['modified_peptides'].append(scan['matches'][0]['modified_peptide']) #list of modified sequences
#
#                        else:
#                            if scan['matches'][0]['modified_peptide'] not in  peptide_entry['mods']:
#                                peptide_entry['mods'].append(scan['matches'][0]['modified_peptide'])
#                                
#                            if (scan['start_scan'], scan['matches'][0]['modified_peptide'], scan['assumed_charge'], scan['matches'][0]['peptide_prev_aa'], scan['matches'][0]['peptide_next_aa'], scan['precursor_neutral_mass'], scan['retention_time_sec'], scan['matches'][0]['num_tot_proteins'], scan['matches'][0]['num_matched_ions'], scan['matches'][0]['tot_num_ions'], scan['matches'][0]['num_tol_term'], scan['matches'][0]['calc_neutral_pep_mass'], scan['matches'][0]['num_missed_cleavages'], scan['matches'][0]['massdiff'], scan['matches'][0]['expect'], scan['matches'][0]['deltacn'] if 'deltacn' in scan['matches'][0] else "-", scan['matches'][0]['spscore'] if 'spscore' in scan['matches'][0] else "-", scan['matches'][0]['xcorr'] if 'xcorr' in scan['matches'][0] else "-", scan['matches'][0]['num_matched_peptides'] if 'num_matched_peptides' in scan['matches'][0] else "-", scan['matches'][0]['probability'], scan['matches'][0]['peptideprophet_probability']) not in peptide_entry['counts'][rep]['scans']:
#                                peptide_entry['counts'][rep]['scans'].append((scan['start_scan'], scan['matches'][0]['modified_peptide'], scan['assumed_charge'], scan['matches'][0]['peptide_prev_aa'], scan['matches'][0]['peptide_next_aa'], scan['precursor_neutral_mass'], scan['retention_time_sec'], scan['matches'][0]['num_tot_proteins'], scan['matches'][0]['num_matched_ions'], scan['matches'][0]['tot_num_ions'], scan['matches'][0]['num_tol_term'], scan['matches'][0]['calc_neutral_pep_mass'], scan['matches'][0]['num_missed_cleavages'], scan['matches'][0]['massdiff'], scan['matches'][0]['expect'], scan['matches'][0]['deltacn'] if 'deltacn' in scan['matches'][0] else "-", scan['matches'][0]['spscore'] if 'spscore' in scan['matches'][0] else "-", scan['matches'][0]['xcorr'] if 'xcorr' in scan['matches'][0] else "-", scan['matches'][0]['num_matched_peptides'] if 'num_matched_peptides' in scan['matches'][0] else "-", scan['matches'][0]['probability'], scan['matches'][0]['peptideprophet_probability']))
#                                
#                                peptide_entry['counts'][rep]['modified_peptides'].append(scan['matches'][0]['modified_peptide'])
#
#                                if peptide_entry['counts'][rep]['best_prob'] < scan['matches'][0]['probability']:
#                                    peptide_entry['counts'][rep]['best_prob'] = scan['matches'][0]['probability']

    title_key_pairs = [
        ('group', 'group_number'),
        ('sibling', 'group_sibling_id'),
        # ('seqid', 'acc'),
        ('protein', 'protein_name'),
        # ('uniprot_id', 'id'),
        # ('url', 'link'),
        ('description', 'description'),
        ('gene', 'gene'),
        ('length', 'prot_length'),
        ('percent_coverage', 'percent_coverage'),
        ('probability', 'probability'),
        ('percent_spectra', 'percent_spectra'),
        # ('organism', 'organism'),
        ('other_seqids', 'other_seqids'),
        ('peptides', 'unique_stripped_peptides')
        # ('peptide_list', 'peptides')
        ]
    for rep in reps:
        title_key_pairs.append((rep + 'exclusive_unique_peptide_count', (rep, 'unique_peptide'))) 
        title_key_pairs.append((rep + 'exclusive_unique_spectrum_count', (rep, 'unique_spectrum')))
        title_key_pairs.append((rep + 'total_spectrum_count', (rep, 'total_spectrum')))
   
    cluster_key_pairs = [
        ('group', 'group_number'),
        ('sibling', 'protein_group_id'),
        ('protein', 'protein_list'),
        ('description', ''),
        ('gene', ''),
        ('length', ''),
        ('percent_coverage', ''),
        ('probability', ''),
        #('exclusive_unique_peptide_count', 'unique_peptide'),
        #('exclusive_unique_spectrum_count', 'unique_spectrum'),
        #('total_spectrum_count', 'total_spectrum'),
        ('percent_spectra', ''),
        ('other_seqids', ''),
        ('peptides', '')
        ]
    for rep in reps:
        cluster_key_pairs.append((rep + 'exclusive_unique_peptide_count', (rep, 'unique_peptide')))
        cluster_key_pairs.append((rep + 'exclusive_unique_spectrum_count', (rep, 'unique_spectrum')))
        cluster_key_pairs.append((rep + 'total_spectrum_count', (rep, 'total_spectrum')))
    
    headings = [title for title, key in title_key_pairs]
    rows = []
    
    #Protein Quantitative Report
    protein_quant_pairs = [
            ('Group ID', 'group_number'),
            ('Protein Group ID', 'protein_group_id'),
            ('Protein Name', 'description'),
            ('Protein Accession Number', 'protein_name')
    ]
    for rep in reps:
        protein_quant_pairs.append((rep + '_' + protein_quant_count + '_count', (rep, protein_quant_count)))

    protein_quant_header = [title for title, key in protein_quant_pairs]
    protein_quant_rows = []
    
    #Protein Report
    protein_pairs = [
            ('Group ID', 'group_number'),
            ('Protein Group ID', 'protein_group_id'),
            ('Protein Name', 'description'),
            ('Protein Accession Number', 'protein_name'),
            ('Length', 'prot_length'),
            ('Percent Coverage', 'percent_coverage'),
            #('Total Number of Peptides', 'total_number_peptides'),
            #('Total Number of Distinct Peptides', 'total_number_distinct_peptides'),
            ('Confidence', 'confidence'),
            ('Percent of Spectrum IDs', 'pct_spectrum_ids'),
            ('Number of Indistinguishable Proteins', 'n_indistinguishable_proteins'),
            ('Protein Identification Probability', 'probability'), #protxml
            ('Exclusive Unique Peptide Count', 'unique_peptide'),
            ('Exclusive Unique Spectrum Count', 'unique_spectrum'),
            ('Total Spectrum Count', 'total_spectrum')
    ]
    
    protein_report_header = [title for title, key in protein_pairs]
    protein_report_header.insert(0, 'MS/MS Sample Name')
    protein_report_rows = []
    
    for group_id, group in protxml_groups.items():
        row = []
        prot_quant_row = []
        
        skip_protein = False
        num_proteins = len(group['cluster']['protein_list']) # if cluster contains one protein, only output cluster
        if num_proteins == 1 or cluster == True:
            skip_protein = True
            
        # CLUSTER LEVEL
        for title, key in cluster_key_pairs:
            if isinstance(key, tuple):
                row.append(len(group['cluster']['counts'][key[0]][key[1]]) if key[0] in group['cluster']['counts'] else "0")
            else:
                row.append(group['cluster'][key] if key in group['cluster'] else "")
        rows.append(row)
        
        #protein quant
        for title, key in protein_quant_pairs:
            if isinstance(key, tuple):
                prot_quant_row.append(len(group['cluster']['counts'][key[0]][key[1]]) if key[0] in group['cluster']['counts'] else "0")
            else:
                prot_quant_row.append(group['cluster'][key] if key in group['cluster'] else "-")
        protein_quant_rows.append(prot_quant_row)
        
        #protein report
#        for rep in reps:
#            prot_report_row = [rep]
#            for title, key in protein_pairs:
#                if (key == 'unique_peptide' or key == 'unique_spectrum' or key == 'total_spectrum'):
#                    prot_report_row.append(len(group['cluster']['counts'][rep][key]) if rep in group['cluster']['counts'] else "0")
#                elif key == 'probability':
#                    prot_report_row.append(group[key] if key in group else "")
#                else:
#                    prot_report_row.append(group['cluster'][key] if key in group['cluster'] else "-")
#            protein_report_rows.append(prot_report_row)

        # PROTEIN LEVEL
        for protein in group['proteins']:
            row = []
            prot_quant_row = []
            prot_report_row = []
            for title, key in title_key_pairs:
                if isinstance(key, tuple):
                    row.append(len(protein['counts'][key[0]][key[1]]) if key[0] in protein['counts'] else "0")
                #if (key == 'unique_peptide' or key == 'unique_spectrum' or key == 'total_spectrum'):
                #  row.append(len(protein[key]) if key in protein else "")
                else:
                    row.append(protein[key] if key in protein else "-")
            rows.append(row)

            #protein quant
            if skip_protein == False:
                for title, key in protein_quant_pairs:
                    if isinstance(key, tuple):
                        prot_quant_row.append(len(protein['counts'][key[0]][key[1]]) if key[0] in protein['counts'] else "0")
                    else:
                        prot_quant_row.append(protein[key] if key in protein else "-")
                protein_quant_rows.append(prot_quant_row)
    
            #protein report
            for rep in reps:
                prot_report_row = [rep]
                skip = False # only output proteins/clusters with total_spectrum > 0
                
                #protein report cluster entry
                if protein['group_sibling_id'] == 'a':
                    for title, key in protein_pairs:
                        if key == 'total_spectrum':
                            if rep in group['cluster']['counts']:
                                if len(group['cluster']['counts'][rep][key]) == 0:
                                    skip = True
                            else:
                                skip = True
                        if (key == 'unique_peptide' or key == 'unique_spectrum' or key == 'total_spectrum'):
                            prot_report_row.append(len(group['cluster']['counts'][rep][key]) if rep in group['cluster']['counts'] else "0")
                        elif key == 'probability':
                            prot_report_row.append(group[key] if key in group else "")
                        elif key in ['prot_length', 'percent_coverage', 'total_number_peptides', 'total_number_distinct_peptides', 'confidence', 'pct_spectrum_ids', 'n_indistinguishable_proteins']:
                             prot_report_row.append(protein[key] if key in protein else "-")
                        else:
                            prot_report_row.append(group['cluster'][key] if key in group['cluster'] else "-")
                    if skip == False:
                        protein_report_rows.append(prot_report_row)
                
                #protein report protein entry
                if skip_protein == False:
                    prot_report_row = [rep]
                    for title, key in protein_pairs:
                        if key == 'total_spectrum':
                            if rep in protein['counts']:
                                if len(protein['counts'][rep][key]) == 0:
                                    skip = True
                            else:
                                skip = True
                        if (key == 'unique_peptide' or key == 'unique_spectrum' or key == 'total_spectrum'):
                            prot_report_row.append(len(protein['counts'][rep][key]) if rep in protein['counts'] else "0")
                        elif (key == 'protein_group_id'):
                            prot_report_row.append(group['cluster'][key] if key in group['cluster'] else "-")
                        else:
                            prot_report_row.append(protein[key] if key in protein else "-")
                    if skip == False:
                        protein_report_rows.append(prot_report_row)
    
    #Peptide Quantitative Report
    peptide_quant_pairs = [
            ('Protein Group ID', 'group_id'),
            ('Protein Group Name', 'group_name'),
            ('Protein Accession Numbers', 'protein_list')
    ]

    peptide_quant_header = [title for title, key in peptide_quant_pairs]
    peptide_quant_header.append('Peptide Sequence')
    peptide_quant_header.append('Modified Sequence')
    for rep in reps:
        peptide_quant_header.append(rep)
    peptide_quant_rows = []
    
    #Peptide Report
    peptide_report_pairs = [
            ('Protein Group ID', 'group_id'),
            ('Protein Group Name', 'group_name'),
            ('Protein Accession Number', 'protein_list')
            ]
    peptide_report_scan_pairs = [
            #('Previous Amino Acid', 'prev_aa'), 
            #('Next Amino Acid', 'next_aa'),
            ('Number of Enzymatic Termini', 'n_enzymatic_termini'),
            ('Unique', 'is_nondegenerate_evidence'),
            #('Number of Total Proteins', 'num_tot_proteins'),
            #('Number of Missed Cleavages', 'num_missed_cleavages'),
            ('Best Peptide Identification Probability', 'best_prob')
    ]

    peptide_report_header = [title for title, key in peptide_report_pairs]
    peptide_report_header.insert(0, 'MS/MS Sample Name')
    #peptide_report_header.append('Peptide Sequence')
    peptide_report_header += [title for title, key in peptide_report_scan_pairs]
    peptide_report_rows = []
    
    #Spectrum Report
    spectrum_report_scan_header = [
            'Scan',
            'Peptide Sequence',
            'Modified Sequence',
            'Charge State',
            #'Modifications',
            'Previous Amino Acid',
            'Next Amino Acid',
            'Precursor Neutral Mass',
            'Retention Time Sec',
            'Number of Total Proteins',
            'Number of Matched Ions',
            'Total Number of Ions',
            'Number of Total Terms',
            'Calculated Neutral Peptide Mass',
            'Number of Missed Cleavages',
            'Mass Difference',
            'Expectation',
            'Deltacn', 
            'Spscore',
            'Xcorr',
            'Number of Matched Peptides',
            'Interprophet Identification Probability',
            'Peptideprophet Identification Probability'
    ]
    
    spectrum_report_header = [title for title, key in peptide_report_pairs]
    spectrum_report_header.insert(0, 'MS/MS Sample Name')
    spectrum_report_header += spectrum_report_scan_header
    peptide_report_header = peptide_report_header[:-3] + spectrum_report_scan_header[:-2] + peptide_report_header[-3:] + spectrum_report_scan_header[-2:]
    spectrum_report_rows = []

    for peptide_name, peptide in peptide_dict.items():
        
        pep_quant_row = []
        #pep_report_row = []
        
        for title, key in peptide_quant_pairs:
            if key == 'protein_list':
                pep_quant_row.append(','.join(peptide[key]) if key in peptide else "-")
            else:
                pep_quant_row.append(peptide[key] if key in peptide else "-")
                #pep_report_row.append(peptide[key] if key in peptide else "-")
        pep_quant_row.append(peptide_name)
        #pep_report_row.append(peptide_name)
        
        counter_list = []
        
        for rep in reps:
            counter_list.append(Counter(peptide['counts'][rep]['modified_peptides']) if rep in peptide['counts'] else 0)
            
            if rep in peptide['counts']:
                
                peptide['counts'][rep]['scans'].sort(key=lambda r:r[-2], reverse=True)
                pep_mod_list = []
                
                #pep_row = []
                #for title, key in peptide_report_scan_pairs:
                #    pep_row.append(peptide['counts'][rep][key] if key in peptide['counts'][rep] else "0")
                #peptide_report_rows.append([rep] + pep_quant_row[0:4] + pep_row)
                
                for scan in peptide['counts'][rep]['scans']:
                    
                    scan_row = [x for x in scan]
                    scan_row.insert(1, peptide_name)
                    
                    if (scan_row[1], scan_row[2]) not in pep_mod_list:
                        if scan_row[1] == "FTGSEIR":
                            print (scan_row[1] + ', '+ scan_row[2])
                        
                        pep_mod_list.append((scan_row[1], scan_row[2]))
                        pep_row = []
                        for title, key in peptide_report_scan_pairs:
                            pep_row.append(peptide['counts'][rep][key] if key in peptide['counts'][rep] else "0")
                        peptide_report_rows.append([rep] + pep_quant_row[0:3] + pep_row[:-3] + scan_row[:-2] + pep_row[-3:] + scan_row[-2:])
                        
                    #scan_row = [x for x in scan]
                    #scan_row.insert(2, peptide_name)
                    spectrum_report_rows.append([rep] + pep_quant_row[0:3] + scan_row) 

        
        for mod in peptide['mods']:
            pep_mod_row = [mod]
            
            for counter_dict in counter_list:
                if counter_dict != 0:
                    pep_mod_row.append(counter_dict[mod] if mod in counter_dict else 0)
                else:
                    pep_mod_row.append(0)
            peptide_quant_rows.append(pep_quant_row + pep_mod_row)
        
#        for rep in reps:
            
#            if rep in peptide['counts']:
#                for key, value in Counter(peptide['counts'][rep]['scans']).items():
#                    mods.append(key)
#                    
#                    pep_quant_row.append(len(peptide['counts'][rep]['scans']) if rep in peptide['counts'] else 0)
            
#            if rep in peptide['counts']:
#                
#                peptide['counts'][rep]['scans'].sort(key=lambda r:r[-2], reverse=True)
#                pep_mod_list = []
#                
#                #pep_row = []
#                #for title, key in peptide_report_scan_pairs:
#                #    pep_row.append(peptide['counts'][rep][key] if key in peptide['counts'][rep] else "0")
#                #peptide_report_rows.append([rep] + pep_quant_row[0:4] + pep_row)
#                
#                for scan in peptide['counts'][rep]['scans']:
#                    
#                    scan_row = [x for x in scan]
#                    scan_row.insert(1, peptide_name)
#                    
#                    if (scan_row[1], scan_row[2]) not in pep_mod_list:
#                        if scan_row[1] == "FTGSEIR":
#                            print (scan_row[1] + ', '+ scan_row[2])
#                        
#                        pep_mod_list.append((scan_row[1], scan_row[2]))
#                        pep_row = []
#                        for title, key in peptide_report_scan_pairs:
#                            pep_row.append(peptide['counts'][rep][key] if key in peptide['counts'][rep] else "0")
#                        peptide_report_rows.append([rep] + pep_quant_row[0:3] + pep_row[:-3] + scan_row[:-2] + pep_row[-3:] + scan_row[-2:])
#                        
#                    #scan_row = [x for x in scan]
#                    #scan_row.insert(2, peptide_name)
#                    spectrum_report_rows.append([rep] + pep_quant_row[0:3] + scan_row) 
#        #peptide_quant_rows.append(pep_quant_row)
    
    #SUMMARY
    summary_report_pairs = [
            ('Total Protein Groups', 'groups'),
            ('Total Proteins', 'proteins'),
            ('Total Stripped Peptides', 'peptides'),
            ('Total Modified Peptides', 'modified'),
            ('Total Spectra', 'spectra')
            ]

    summary_report_header = [title for title, key in summary_report_pairs]
    summary_report_header.insert(0, 'MS/MS Sample Name')
    summary_report_rows = []

    for rep, entry in summary_dict.items():
        summary_report_row = [rep]
        summary_report_row += [len(entry[value]) if value in entry else "-" for key,value in summary_report_pairs]
        summary_report_rows.append(summary_report_row)

    rows.sort(key=lambda r:r[0])
    rows.insert(0, headings)
    protein_quant_rows.insert(0, protein_quant_header)
    protein_report_rows.sort(key=lambda r:r[0])
    protein_report_rows.sort(key=lambda r:r[1])
    protein_report_rows.insert(0, protein_report_header)
    peptide_quant_rows.insert(0, peptide_quant_header)
    peptide_report_rows.insert(0, peptide_report_header)
    spectrum_report_rows.insert(0, spectrum_report_header)
    summary_report_rows.sort(key=lambda r:r[0])
    summary_report_rows.insert(0, summary_report_header)

    ####################
    benchmark5 = time.time()
    
    print ('finished output files')
    print (benchmark5 - benchmark4)
    ####################

    print "-"*60
    print "WRITING:", os.path.abspath(proteins_csv)
    csv_writer = csv.writer(open(proteins_csv, 'wb'))
    for row in rows:
        csv_writer.writerow(row)

    print "WRITING:", os.path.abspath(protein_quant_csv)
    csv_writer = csv.writer(open(protein_quant_csv, 'wb'))
    for row in protein_quant_rows:
        csv_writer.writerow(row)

    print "WRITING:", os.path.abspath(protein_report_csv)
    csv_writer = csv.writer(open(protein_report_csv, 'wb'))
    for row in protein_report_rows:
        csv_writer.writerow(row)

    print "WRITING:", os.path.abspath(peptide_quant_csv)
    csv_writer = csv.writer(open(peptide_quant_csv, 'wb'))
    for row in peptide_quant_rows:
        csv_writer.writerow(row)
    
    print "WRITING:", os.path.abspath(peptide_csv)
    csv_writer = csv.writer(open(peptide_csv, 'wb'))
    for row in peptide_report_rows:
        csv_writer.writerow(row)
    
    print "WRITING:", os.path.abspath(spectrum_csv)
    csv_writer = csv.writer(open(spectrum_csv, 'wb'))
    for row in spectrum_report_rows:
        csv_writer.writerow(row)
    
    summary_csv = prefix + '_' + str(pepxml_reader.prob_cutoff) + '_' + 'summary_report.csv'
    
    print "WRITING:", os.path.abspath(summary_csv)
    csv_writer = csv.writer(open(summary_csv, 'wb'))
    for row in summary_report_rows:
        csv_writer.writerow(row)
    
    
    ####################
    benchmark6 = time.time()
    
    print ('total time')
    print (benchmark6 - start)
    ####################

def get_opt(argv):

    if len(argv) == 0:
        print("Error: No Protein XML file specified\n")
        usage()
        
    protxml_file = ''
    root = ''
    pepxml_file = ''
    mayu_file = ''
    
    protein_quant_count = 'total_spectrum'   
    fdr_method = 'mayu'
    
    protein_probability_cutoff = 0.95 # or None
    peptide_probability_cutoff = 0.95 # or None
    protein_error_cutoff = 0.01 # or None
    peptide_error_cutoff = 0.01 # or None
    spectrum_error_cutoff = 0.01
    
    cluster = False
    unique = False
    help = False
    error = False
    
    # parse the command line inputs
    try:
        opts, args = getopt.getopt(argv, "i:j:m:p:q:f:a:b:r:e:s:cuh",
                                   ["protXML=","pepXML=","mayuFile=","prefix=","quant=","filter=","protProb=","pepProb=","protFDR=","pepFDR=","psmFDR=","cluster","unique","help"])
        print ('debug')                            
        print (opts)
    
    except getopt.GetoptError:
        help = True
        opts = (("", ""),)

    if help:
        print_help()

    for opt, arg in opts:
        if opt in ("-i", "--protXML"):
            protxml_file = str(arg).strip()
            root = os.path.splitext(os.path.splitext(protxml_file)[0])[0]
            if pepxml_file == '':
                pepxml_file = root + '.pep.xml'
        elif opt in ("-j", "--pepXML"):
            pepxml_file = str(arg).strip()
        elif opt in ("-m", "--mayuFile"):
            mayu_file = str(arg).strip()
        elif opt in ("-p", "--prefix"):
            root = str(arg).strip()
        elif opt in ("-q", "--quant"):
            if int(arg) == 1:
                protein_quant_count = 'unique_peptide'
            elif int(arg) == 2:
                protein_quant_count = 'unique_spectrum'
            elif int(arg) == 3:
                protein_quant_count = 'total_spectrum'
            else:
                print("Error: The quant options for the Protein Quantitative Report are:")
                print("    1      Exclusive Unique Peptide Count")
                print("    2      Exclusive Unique Spectrum Count")
                print("    3      Total Spectrum Count")
                error = True
        elif opt in ("-f", "--filter"):
            if str(arg).strip() == 'mayu':
                fdr_method = 'mayu'
            elif str(arg).strip() == 'iFDR':
                fdr_method = 'iFDR'
            elif str(arg).strip() == 'noFDR':
                fdr_method = 'noFDR'
            elif str(arg).strip() == 'mayu_iFDR':
                fdr_method = 'mayu_iFDR'    
            else:
                print("Error: The calculation options for the Spectrum and Protein Probability Cutoffs are:")
                print("    mayu   Mayu is used to compute Peptide Probability Cutoff from protein, peptide, and spectrum FDR values")
                print("    iFDR   Spectrum and Protein Probability Cutoffs are determined from spectrum and protein FDR values without using Mayu")
                print("    noFDR  spectrum and Protein Probability Cutoffs are used directly")
                print("    mayu_iFDR  Mayu used to compute Peptide Probability Cutoff, protFDR used to compute Protein Probability Cutoff\n")
                error = True
        elif opt in ("-a", "--protProb"):
            protein_probability_cutoff = arg
        elif opt in ("-b", "--pepProb"):
            peptide_probability_cutoff = arg
        elif opt in ("-r", "--protFDR"):
            protein_error_cutoff = arg
        elif opt in ("-e", "--pepFDR"):
            peptide_error_cutoff = arg
        elif opt in ("-s", "--psmFDR"):
            spectrum_error_cutoff = arg
        elif opt in ("-c", "--cluster"):
            cluster = True
        elif opt in ("-u", "--unique"):
            unique = True
        elif opt in ("-h", "--help"):
            print_help()
            
    # check if input file exists
    if protxml_file == '':
        print("Error: No Protein XML file specified\n")
        usage()
    if not os.path.exists(protxml_file):
        print("Error: %s does not exist.\n") % os.path.abspath(protxml_file)
        usage()
    if not os.path.exists(pepxml_file):
        print("Error: Cound not find Peptide XML file %s\n") % os.path.abspath(pepxml_file)
        usage()
    if fdr_method == 'mayu':
        if mayu_file == '':
            print("Error: Please specify Mayu file\n")
            usage()
        elif not os.path.exists(mayu_file):
            print("Error: Cound not find Mayu file %s\n") % os.path.abspath(mayu_file)
            usage()
    if error == True:
        usage()
    
    #OUTPUT FILES    
    outfile = root + 'test.prot.csv'
    peptide_quant_csv = root + '.peptide_quantitative_report.csv'
    peptide_csv = root + '.peptide_report.csv'
    protein_quant_csv = root + '.protein_quantitative_report.csv'
    protein_report_csv = root + '.protein_report.csv'
    spectrum_csv = root + '.spectrum_report.csv'
    #summary_csv = root + '.summary_report.csv'

    return protxml_file, pepxml_file, mayu_file, outfile, peptide_quant_csv, peptide_csv, protein_quant_csv, protein_report_csv, spectrum_csv, root, protein_quant_count, protein_probability_cutoff, protein_error_cutoff, peptide_probability_cutoff, peptide_error_cutoff, spectrum_error_cutoff, fdr_method, cluster, unique

def usage():
    print("Usage:      protxml2csv.py [-i protXML] [-j pepXML] [-m mayuFile] [-p prefix] [-q quant] [-f filter] [-a protProb] [-b pepProb] [-r protFDR] [-e pepFDR] [-s psmFDR] [-c] [-u] \n")
    sys.exit()

def print_help():
    print("protxml2csv.py")
    print("---------------------------------------------------------------------------------------------")
    print("Usage:      protxml2csv.py [-i protXML] [-j pepXML] [-m mayuFile] [-p prefix] [-q quant] [-f filter] [-a protProb] [-b pepProb] [-r protFDR] [-e pepFDR] [-s psmFDR] [-c] [-u] \n")
    print("            Performs spectral counting on TPP Protein XML file and corresponding Peptide XML file and outputs results in CSV format\n")
    print("                       -i [--protXML]   STR    TPP Protein XML filename to perform spectral counting calculations on")
    print("            (optional) -j [--pepXML]    STR    Corresponding TPP Peptide XML filename [Default: searches for pepXML file with matching protXML filename]")
    print("            (optional) -m [--mayuFile]  STR    Mayu file used to calculate the Peptide and Protein Probability Cutoffs. Only used when mayu option specified for [-f] flag")
    print("            (optional) -p [--prefix]    STR    Prefix of the output files [Default: same as protXML filename]")
    print("            (optional) -q [--quant]     INT    Spectral counting method used in Protein Quantitative Report [Default: Total Spectrum Count]")
    print("                                               Available options are:")
    print("                                                 1      Exclusive Unique Peptide Count")
    print("                                                 2      Exclusive Unique Spectrum Count")
    print("                                                 3      Total Spectrum Count\n")
    print("            (optional) -f [--filter]    STR    Filtering method used to calculate Peptide and Protein Probability Cutoffs from FDR values [Default: mayu]")
    print("                                               Available options are:")
    print("                                                 mayu   Mayu is used to compute Peptide Probability Cutoff from protein, peptide, and spectrum FDR values")
    print("                                                 iFDR   Protein and Peptide Probability Cutoffs are determined from FDR values without using Mayu")
    print("                                                 noFDR  Protein and Peptide Probability Cutoffs are used directly\n")
    print("                                                 mayu_iFDR  Mayu used to compute Peptide Probability Cutoff, protFDR used to compute Protein Probability Cutoff\n")
    print("            (optional) -a [--protProb]  FLOAT  Protein Probability Cutoff [Default: 0.95]")
    print("            (optional) -b [--pepProb]   FLOAT  Peptide Probability Cutoff [Default: 0.95]")
    print("            (optional) -r [--protFDR]   FLOAT  Protein FDR [Default: 0.01]")
    print("            (optional) -e [--pepFDR]    FLOAT  Peptide FDR [Default: 0.01]")
    print("            (optional) -s [--psmFDR]    FLOAT  Spectrum FDR [Default: 0.01]\n")
    print("            (optional) -c [--cluster]          If flag is specified, only outputs protein clusters")
    print("            (optional) -u [--unique]           If flag is specified, only inclues unique peptides in spectral counting calculation")

    sys.exit()

if __name__ == "__main__":


    params = {
        'is_skip_no_unique': False,
    }
    convert_to_csv(params)




