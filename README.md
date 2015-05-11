# noncodingRNAfinder.py

Developed by molecular physiology, lab101, plant science, cambridge, UK

Program to find non coding RNAs (Chris and Ivans project)

This is a pipeline to identify long non coding RNAs from plant transcriptomes. It uses PLEK are it's core algorithm to
define non coding RNAs. It has additional filters to remove tRNAs, rRNAs, and precursor miRNAs, using a series of
databases. It uses a blast search against a high quality uniprot protein database to remove potential coding sequences
as an additional step to PLEK.

It is not finished.

