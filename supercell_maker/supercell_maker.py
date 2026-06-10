#!/usr/bin/env python3

from ase.io import read, write
import numpy as np
import sys
import os

def cif_to_supercell_xyz(infile, outfile_xyz=None, outfile_cell=None, rep=(1, 1, 1)):
    # Read CIF
    atoms = read(infile)

    # Build supercell
    supercell = atoms * rep

    # Extract 3×3 Cartesian cell matrix
    cell = supercell.get_cell()

    # Default filenames
    base, _ = os.path.splitext(infile)
    if outfile_xyz is None:
        outfile_xyz = f"{base}_{rep[0]}x{rep[1]}x{rep[2]}.xyz"
    if outfile_cell is None:
        outfile_cell = f"{base}_{rep[0]}x{rep[1]}x{rep[2]}_cell.txt"

    # Write XYZ
    write(outfile_xyz, supercell, format="xyz")

    # Write cell matrix
    with open(outfile_cell, "w") as f:
        f.write("cell_cart\n")
        for row in cell:
            f.write(" ".join(f"{float(x):.6f}" for x in row) + "\n")


    print(f"Wrote XYZ file: {outfile_xyz}")
    print(f"Wrote cell matrix file: {outfile_cell}")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} input.cif [output.xyz] [cell.txt]")
        sys.exit(1)

    infile = sys.argv[1]
    outfile_xyz = sys.argv[2] if len(sys.argv) > 2 else None
    outfile_cell = sys.argv[3] if len(sys.argv) > 3 else None

    # Ask user for supercell dimensions
    print("Enter supercell replication as three integers (e.g., 5 5 3):")
    try:
        nx, ny, nz = map(int, input().split())
    except Exception:
        print("Invalid input. Expected three integers like: 5 5 3")
        sys.exit(1)

    rep = (nx, ny, nz)

    cif_to_supercell_xyz(infile, outfile_xyz, outfile_cell, rep=rep)
