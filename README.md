# lncRNAidentifier.py

Program to find non coding RNAs
Developed in Python by molecular physiology at plant science cambridge 
Authors: Chris and Ivan

This is a pipeline to identify long non coding RNAs from diverse monocot and dicot plant transcriptomes. It uses PLEK as its core algorithm to define non coding RNAs. It has additional filters to remove tRNAs, rRNAs, and precursor miRNAs, using a series of arabidopsis & ensembl databases. It uses a blast search against a high quality uniprot protein database to remove potential coding proteins as an additional step to PLEK. This triple lock serves to give a high confidence series of transcripts for being lncRNAs.

Update 18.05.14 - Program is now operational (v1)

