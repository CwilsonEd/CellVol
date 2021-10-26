# CellVol
Repository for the CellVol code

For a full introduction to the software see DOI:xxxx

This code has been designed to run with the standard distrubution of the Cambridge Crystallographic Data Centre. 
When installing the program no additional libraries are required outside of the standard CCDC install.

This code features functionality to run both through Mercury and directly from the command line on single or multiple cif files.

Command line inputs should be formulated >>python CellVol cifname pointsperrun no.runs precision
For example:
   >>python GLYCIN05 1 16 0.1

will search the directory for GLYCIN05 before taking this CIF from the CSD and running it with 16 runs of a million points to a precision of 0.1%. 

When running from Mercury the program runs on the current entry in the visualisation window.

This code is reguarly updated and maintained by the Parsons group.
If you have any problems running the software or encounter any bugs, please contact professor Parsons: simon.parsons@ed.ac.uk
