This program can be used to generate the coordinates for a supercell needed for an

INTERFACES calculation. The program required python and ASE to be installed.



start the program with:



py supercell\_maker.py



It will then ask for the replications of the a, b, and c directions requested.



The output will print out the cell matrix that can be pasted directly into your

INTERFACES input file in addition to an xyz file. This will need to be converted

into a mol2 format using software such as openbabel, Avogadro, Mercury, etc to add

and save the bonds to the structure.

