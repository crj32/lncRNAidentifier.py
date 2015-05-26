# lncRNAidentifier.py

Program to find non coding RNAs

Developed in Python by molecular physiology at plant science cambridge 

Authors: Chris and Ivan

This is a pipeline to identify long non coding RNAs from diverse monocot and dicot plant transcriptomes. It uses PLEK as its core algorithm to define non coding RNAs, which is trained using human data. It has additional filters to remove tRNAs, rRNAs, and precursor miRNAs, using a series of usearch databases explained in the script. It uses a blast search against a high quality uniprot protein database to remove potential coding proteins as an additional step to PLEK. It removes short transcripts (< 200) and those with maxiumum ORF > 100.

Update 18.05.14 - Program is now operational (v1)

Update 26.05.14 - Program now properly parses command line options and has optional blast against uniprot

