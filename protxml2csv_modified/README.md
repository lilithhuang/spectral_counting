

====================================================
PEPTO - proteomic analysis and visualization scripts
====================================================

This is a library of scripts, ostensibly designed to generate visualizations of PEPTOGRAPHS, a protocol for peptide-analysis of proteomic mass-spec experiments, designed by Benjamin Cravatt. There are also many convenience scripts for manipulating mass-spec data formts

## Dependencies

MAKO
UNIPROT

## Peptographs

Peptographs is an extremely powerful way of analysis proteomic data. PEPTO generates interactive Peptographs with integrated protein-lists/top-spectra displays in HTML format. This is designed to be extremely portable - you can send it as an attachment, or dispaly on any web server. Different scripts are provided for the different types of input proteomics search results.

## Python parsers for major proteomic platforms

PEPTO is, unashamedly, a Python library that provides parsers for MASCOT, X!Tandem, MAXQUANT and TPP search results. These search engines report likely peptides, taken from a reference databases of known protein sequences, that match the mass-spectra of an experiment. 

PEPTO reads these results into a unified Python data structure organized in terms of SeqIDs.

## Label-free protein identification

protxml_to_csv.py is a convenience script that converts PROT.XML into EXCEL-friendly CSV files. In label-free protein identification analyses, MASCOT and X!Tandem results are post-processed using the TPP to identify the most likely proteins from the set of identified proteins. In such analyses, a list of proteins with their associated probabilities are reported in the TPP's PROTXML format. 

Where the relevant SeqIDs are reconginzably Uniprot identifiers, the script will generate a useful HTTP link in the CSV file, and extract relevant metadata, such as GENE and ORGANISM.

## TPP's ProteinProphet wrapper

TPP is a complex package, with a fully-featured interactive web-interface. The most useful module w.r.t. to PEPTO is the ProteinProphet analysis tool for inducing the most like proteins from a peptide analysis of mass-spec data. 

protein_prophet.py is a python script that takes MASCOT data, and a reference database, and generates the PEPXML and PROTXML file for ProteinProphet analysis. This allows the TPP to be embedded in a Unix-style pipeline.

## Sequence overview

coverage_html.py generates a simple view of the peptide coverage of protein predictions in an attractive HTML file. 




