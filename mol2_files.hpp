#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
using namespace std;

int write_mol2(char *filename, int N_atoms, int N_bonds, int *atom_id, char (*element)[3], double (*xyz)[3],char (*atom_type)[8], int *bond_id, int *ori_atom_id, int *tar_atom_id, char (*bond_type)[20]) {
    //This function writes out the structure given by xyz as a mol2 file.
    //All other variables in the mol2 file are taken from the initial reading of the starting structure
    int i;
    FILE *out, *error_file;
    char error_filename[128];

    sprintf(error_filename,"Errors.txt");

    remove(filename);
    out=fopen(filename,"w");

        if(out==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: Mol2 file '%s' could not be written\n", filename);
            fclose(error_file);
            exit(1);
        }
        else{
            fprintf(out,"@<TRIPOS>MOLECULE\n*****\n %d %d 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n", N_atoms,N_bonds);
            for(i=0;i<N_atoms;i++){
            fprintf(out,"%d %s %.4f %.4f %.4f %s     1  UNL1        0.0000\n",atom_id[i], element[i], xyz[i][0], xyz[i][1], xyz[i][2], atom_type[i]);
            }
            fprintf(out,"\n@<TRIPOS>BOND\n");
            for(i=0;i<N_bonds;i++){
                fprintf(out,"%d %d %d %s\n",bond_id[i], ori_atom_id[i], tar_atom_id[i], bond_type[i]);
            }
        fclose(out);
        return 0;
        }
}

void add_surface(char *mol2_filename){
    //This function adds 5 bonded silicon atoms at the material surface
    //These are useful for generating Figures and visualizing the orientation of the complex on the surface
    FILE *fp, *out;
    int i, j, k, line_Atoms, line_Bonds, N_atoms, N_bonds;
    char error_filename[128], buffer[120];
    sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    fp=fopen(mol2_filename,"r");
    if(fp==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: File '%s' could not be found\n", mol2_filename);
        fclose(error_file);
        exit(1);
    }

    i=j=1;
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
            if(strcmp(buffer, "@<TRIPOS>MOLECULE\n")==0){
                fgets(buffer, sizeof(buffer), fp);
                fgets(buffer, sizeof(buffer), fp);
                sscanf(buffer,"%d %d",&N_atoms, &N_bonds);
                i++; j++;
                i++; j++;
            }
            if(strcmp(buffer, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
                }
            if(strcmp(buffer, "@<TRIPOS>BOND\n")==0){
                line_Bonds = j;
            }
            i++; j++;
    }
    fclose(fp);


    char add_surface_filename[100];
    int structure_filename_len=strlen(mol2_filename)-5;
    char base_filename[structure_filename_len+1];

    sprintf(base_filename,"%.*s",structure_filename_len, mol2_filename);
    sprintf(add_surface_filename, "%s_with_surface.mol2", base_filename);
    remove(add_surface_filename);

    fp=fopen(mol2_filename, "r");
    out=fopen(add_surface_filename,"a");

    if(out==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Could not write file '%s'\n", add_surface_filename);
        fclose(error_file);
        exit(1);
    }

    i=1;
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
        if(i==1){
            fprintf(out,"@<TRIPOS>MOLECULE\n*****\n %d %d 0 0 0\n", N_atoms+5, N_bonds+4);
            fgets(buffer, sizeof(buffer), fp);
            fgets(buffer, sizeof(buffer), fp);
            i=i+2;
        }
        else if(i==line_Atoms+N_atoms){
            fprintf(out,"%d Si 0.0000 0.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+1);
            fprintf(out,"%d Si 50.0000 50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+2);
            fprintf(out,"%d Si -50.0000 -50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+3);
            fprintf(out,"%d Si -50.0000 50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+4);
            fprintf(out,"%d Si 50.0000 -50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+5);
            fprintf(out, "%s", buffer);
            i++;
        }
        else{
            fprintf(out, "%s", buffer);
            i++;
        }
    }
    fclose(fp);

    fprintf(out,"%d %d %d 1\n", N_bonds+1, N_atoms+1, N_atoms+2);
    fprintf(out,"%d %d %d 1\n", N_bonds+2, N_atoms+1, N_atoms+3);
    fprintf(out,"%d %d %d 1\n", N_bonds+3, N_atoms+1, N_atoms+4);
    fprintf(out,"%d %d %d 1\n", N_bonds+4, N_atoms+1, N_atoms+5);
    fclose(out);

    remove(mol2_filename);
    rename(add_surface_filename,mol2_filename);
}

void compile_mol2_files(char *base_filename, int num_files){
    //This function combines a sequence of mol2 files into a single mol2 file showing their overlay

    char error_filename[128], buffer[120];
    sprintf(error_filename,"Errors.txt");
    FILE *error_file, *fp, *in, *out;
    int atom_index = 1, i, j, k = 1, line_Molecule, line_Atoms, line_Bonds, N_atoms, N_bonds;
    int structure_filename_len=strlen(base_filename)+20;
    char structure_filename[structure_filename_len];
    sprintf(structure_filename, "%s_struct1.mol2", base_filename);

    for(i=1; i<num_files+1; i++){
        //add surface molecules to each structure before reading
        char structure_filename_temp[100];
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, i);
        add_surface(structure_filename_temp);
    }

    //1) Open structure mol2 file and find Number of Atoms/Bonds & line  numbers for atom and bond sections
    fp=fopen(structure_filename,"r");
    if(fp==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename);
        fclose(error_file);
        exit(1);
    }

    i=1;
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
            if(strcmp(buffer, "@<TRIPOS>MOLECULE\n")==0){
                fgets(buffer, sizeof(buffer), fp);
                fgets(buffer, sizeof(buffer), fp);
                sscanf(buffer,"%d %d",&N_atoms, &N_bonds);
                i++;
                i++;
            }
            else if(strcmp(buffer, "@<TRIPOS>BOND\n")==0){
                line_Bonds = i;
            }
            else if(strcmp(buffer, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
            }
            i++;
    }
    fclose(fp);

    //2) Scan the bond section and store bond information
    int bond_id[N_bonds], ori_atom_id[N_bonds], tar_atom_id[N_bonds];
    char bond_type[N_bonds][20];

    j=1;
    fp=fopen(structure_filename,"r");
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
        if(j==line_Bonds+1){
            for(i=0; i<N_bonds; i++){
                sscanf(buffer,"%d %d %d %s", &bond_id[i], &ori_atom_id[i], &tar_atom_id[i], &bond_type[i]);
                fgets(buffer, sizeof(buffer), fp);
                j++;
            }}
        else{
        j++;
        }
    }
    fclose(fp);

    //3) create the structure overlay file.
    char overlay_filename[128];
    sprintf(overlay_filename, "%s_overlay.mol2", base_filename);
    remove(overlay_filename);
    out=fopen(overlay_filename,"w");

    //4) Add the molecule section
    if(out==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file,"\nstructure overlay file could not be written.\n");
        fclose(error_file);
    }

    else
        fprintf(out,"@<TRIPOS>MOLECULE\n*****\n %d %d 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n", (num_files)*N_atoms, (num_files)*N_bonds);

    fclose(out);

    //5) Add the atom section from each individual file(Loop until m reaches the num_files)
    for(k=1; k<num_files+1; k++){
        char structure_filename_temp[100];
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, k);

        in=fopen(structure_filename_temp, "r");
        out=fopen(overlay_filename, "a");

        if(in==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename_temp);
            fclose(error_file);
        }

        j=1;
        char arguments[8][16];

        while(fgets(buffer, sizeof(buffer), in) != NULL){
            if(j==line_Atoms+1){
                for(i=1; i<N_atoms+1; i++){
                    sscanf(buffer,"%s %s %s %s %s %s %s %s",&arguments[0],&arguments[1],&arguments[2],&arguments[3],&arguments[4],&arguments[5],&arguments[6],&arguments[7]);
                    fprintf(out, "%d %s %s %s %s %s %s %s\n", atom_index,arguments[1],arguments[2],arguments[3],arguments[4],arguments[5],arguments[6],arguments[7]);
                    atom_index++;
                    fgets(buffer, sizeof(buffer), in);
                }
            j++;
            }
            else
                j++;
        }
        fclose(in);
        fclose(out);
    }
    atom_index=1;

    //Add the bond section using the bond information stored from earlier
    out=fopen(overlay_filename, "a");
    fprintf(out, "\n@<TRIPOS>BOND\n");

    for(j=0; j<num_files; j++){
        for(i=0; i<N_bonds; i++){
            fprintf(out, "%d %d %d %s\n", bond_id[i]+(j*N_bonds), ori_atom_id[i]+(j*N_atoms), tar_atom_id[i]+(j*N_atoms), bond_type[i]);
        }
    }
    fclose(out);
}

void compile_all_mol2_files(char *base_filename, int num_primary_files, int num_other_files){
    //This function performs the same task as compile_mol2_files, however it will look for files also starting with other_
    //It also removes the said files after they have been read.
    char error_filename[128], other_base_filename[128],overlay_filename[128],structure_filename_temp[128], buffer[128], arguments[8][16];;
    sprintf(error_filename,"Errors.txt");
    sprintf(other_base_filename,"other_%s_struct",base_filename);
    FILE *error_file,*fp, *in, *out;
    int atom_index = 1, i, j, k = 1, line_Molecule, line_Atoms, line_Bonds;
    int N_atoms, N_bonds;
    int structure_filename_len=strlen(base_filename)+20;
    char structure_filename[structure_filename_len];
    sprintf(structure_filename, "%s_struct1.mol2", base_filename);
    int total_num_files = num_primary_files + num_other_files;

    //1) Open structure mol2 file and find Number of Atoms/Bonds & line for bond section
    fp=fopen(structure_filename,"r");
    if(fp==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename);
        fclose(error_file);
        exit(1);
    }

    i=1;
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
            if(strcmp(buffer, "@<TRIPOS>MOLECULE\n")==0){
                fgets(buffer, sizeof(buffer), fp);
                fgets(buffer, sizeof(buffer), fp);
                sscanf(buffer,"%d %d",&N_atoms, &N_bonds);
                i++;
                i++;
            }
            else if(strcmp(buffer, "@<TRIPOS>BOND\n")==0){
                line_Bonds = i;
            }
            else if(strcmp(buffer, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
            }
            i++;
    }
    fclose(fp);

    //2) Scan the bond section and store bond information
    int bond_id[N_bonds], ori_atom_id[N_bonds], tar_atom_id[N_bonds];
    char bond_type[N_bonds][20];

    fp=fopen(structure_filename,"r");
    j=1;
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
        if(j==line_Bonds+1){
            for(i=0; i<N_bonds; i++){
                sscanf(buffer,"%d %d %d %s", &bond_id[i], &ori_atom_id[i], &tar_atom_id[i], &bond_type[i]);
                fgets(buffer, sizeof(buffer), fp);
                j++;
            }}
        else{
        j++;
        }
    }
    fclose(fp);

    //3) make structure overlay
    sprintf(overlay_filename, "%s_overlay_all.mol2", base_filename);
    remove(overlay_filename);
    out=fopen(overlay_filename,"w");

    //4) Add the molecule section from primary structures
    if(out==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file,"structure overlay file could not be written.\n");
        fclose(error_file);
    }

    else
        fprintf(out,"@<TRIPOS>MOLECULE\n*****\n %d %d 0 0 0\nSMALL\nGASTEIGER\n\n@<TRIPOS>ATOM\n", (total_num_files)*N_atoms, (total_num_files)*N_bonds);

    fclose(out);

    //5) Add the atom section from each individual primary file(Loop until m reaches the num_primary_files)
    for(k=1; k<num_primary_files+1; k++){
        //add surface molecules to each structure before reading...used for .cif file
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, k);
        in=fopen(structure_filename_temp, "r");
        out=fopen(overlay_filename, "a");

        if(in==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: File '%s' not found could not create structure overlay\n", structure_filename_temp);
            fclose(error_file);
        }

        j=1;
        while(fgets(buffer, sizeof(buffer), in) != NULL){
            if(j==line_Atoms+1){
                for(i=1; i<N_atoms+1; i++){
                    sscanf(buffer,"%s %s %s %s %s %s %s %s",&arguments[0],&arguments[1],&arguments[2],&arguments[3],&arguments[4],&arguments[5],&arguments[6],&arguments[7]);
                    fprintf(out, "%d %s %s %s %s %s %s %s\n", atom_index,arguments[1],arguments[2],arguments[3],arguments[4],arguments[5],arguments[6],arguments[7]);
                    atom_index++;
                    fgets(buffer, sizeof(buffer), in);
                }
            j++;
            }
            else
                j++;
        }
        fclose(in);
        fclose(out);
        remove(structure_filename_temp);
    }

    //6) Add the atom section from each individual other file
    for(k=1; k<num_other_files+1; k++){
        sprintf(structure_filename_temp, "%s%d.mol2", other_base_filename, k);
        add_surface(structure_filename_temp);
        in=fopen(structure_filename_temp, "r");
        out=fopen(overlay_filename, "a");

        if(in==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename_temp);
            fclose(error_file);
        }

        j=1;
        while(fgets(buffer, sizeof(buffer), in) != NULL){
            if(j==line_Atoms+1){
                for(i=1; i<N_atoms+1; i++){
                    sscanf(buffer,"%s %s %s %s %s %s %s %s",&arguments[0],&arguments[1],&arguments[2],&arguments[3],&arguments[4],&arguments[5],&arguments[6],&arguments[7]);
                    fprintf(out, "%d %s %s %s %s %s %s %s\n", atom_index,arguments[1],arguments[2],arguments[3],arguments[4],arguments[5],arguments[6],arguments[7]);
                    atom_index++;
                    fgets(buffer, sizeof(buffer), in);
                }
            j++;
            }
            else
                j++;
        }
        fclose(in);
        fclose(out);
        remove(structure_filename_temp);
    }

    //7) Add the bond section using the bond information stored from earlier
    out=fopen(overlay_filename, "a");
    fprintf(out, "\n@<TRIPOS>BOND\n");

    for(j=0; j<total_num_files; j++){
        for(i=0; i<N_bonds; i++){
            fprintf(out, "%d %d %d %s\n", bond_id[i]+(j*N_bonds), ori_atom_id[i]+(j*N_atoms), tar_atom_id[i]+(j*N_atoms), bond_type[i]);
        }
    }
    fclose(out);
}
