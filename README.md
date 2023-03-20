# CellVol
Repository for the CellVol code

For a full introduction to the software see Wilson, C.J.G., Cervenka, T., Wood, P.A. and Parsons, S., 2022. Behavior of Occupied and Void Space in Molecular Crystal Structures at High Pressure. Crystal Growth & Design, 22(4), pp.2328-2341. 

This should be cited in all instances of use.

This code has been designed to run with the standard distrubution of the Cambridge Crystallographic Data Centre. 
When installing the program no additional libraries are required outside of the standard CCDC install. Compatibility last checked Mercury (2022.3.0) and CSD Python API (3.7.9).

This code features functionality to run both through Mercury and directly from the command line on single or multiple cif files.

Command line inputs should be formulated >>python CellVol cifname pointsperrun no.runs precision.

For example:
   >>python GLYCIN05 1 16 0.1

will search the directory for GLYCIN05 before taking this CIF from the CSD and running it with 16 runs of a million points to a precision of 0.1%. 

When running from Mercury the program runs on the current entry in the visualisation window. Mercury currently features a timeout error when scripts take longer tha 180 seconds, the code will still run and output to C:\CellVol_Output

This code is reguarly updated and maintained by the Parsons group.
If you have any problems running the software or encounter any bugs, please contact Professor Parsons: simon.parsons@ed.ac.uk
