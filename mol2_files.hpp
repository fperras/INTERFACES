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
    int atom_index = 1, i = 1, j = 1, k = 1, line_Atoms, line_Bonds;
    int N_atoms, N_bonds;
    char error_filename[128];
    sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    fp=fopen(mol2_filename,"r");
    if(fp==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: File '%s' could not be found\n", mol2_filename);
        fclose(error_file);
        exit(1);
    }
    char line_content[120];

    while(fgets(line_content, sizeof(line_content), fp) != NULL){
            if(strcmp(line_content, "@<TRIPOS>MOLECULE\n")==0){
                fgets(line_content, sizeof(line_content), fp);
                fgets(line_content, sizeof(line_content), fp);
                sscanf(line_content,"%d %d",&N_atoms, &N_bonds);
                i++; j++;
                i++; j++;
            }
            if(strcmp(line_content, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
                }
            if(strcmp(line_content, "@<TRIPOS>BOND\n")==0){
                line_Bonds = j;
            }
            i++; j++;
    }
    k = 1;
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

    while(fgets(line_content, sizeof(line_content), fp) != NULL){
        if(k==1){
            fprintf(out,"@<TRIPOS>MOLECULE\n*****\n %d %d 0 0 0\n", N_atoms+5, N_bonds+4);
            fgets(line_content, sizeof(line_content), fp);
            fgets(line_content, sizeof(line_content), fp);
            k=k+2;
        }
        else if(k==line_Atoms+N_atoms){
            fprintf(out,"%d Si 0.0000 0.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+1);
            fprintf(out,"%d Si 50.0000 50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+2);
            fprintf(out,"%d Si -50.0000 -50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+3);
            fprintf(out,"%d Si -50.0000 50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+4);
            fprintf(out,"%d Si 50.0000 -50.0000 0.0000 Si 1 UNL1 0.0000\n", N_atoms+5);
            fprintf(out, "%s", line_content);
            k++;
        }
        else{
            fprintf(out, "%s", line_content);
            k++;
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

    char error_filename[128];
    sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    FILE *fp, *in, *out;
    int atom_index = 1, i = 1, j = 1, k = 1, m, line_Molecule, line_Atoms, line_Bonds;
    int N_atoms, N_bonds;
    int structure_filename_len=strlen(base_filename)+20;
    char structure_filename[structure_filename_len];
    sprintf(structure_filename, "%s_struct1.mol2", base_filename);

    for(m=1; m<num_files+1; m++){
        //add surface molecules to each structure before reading
        char structure_filename_temp[100];
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, m);
        add_surface(structure_filename_temp);
    }

    //1) Open structure mol2 file and find Number of Atoms/Bonds & line for bond section
    fp=fopen(structure_filename,"r");
    if(fp==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename);
        fclose(error_file);
        exit(1);
    }
    char line_content[120];

    while(fgets(line_content, sizeof(line_content), fp) != NULL){
            if(strcmp(line_content, "@<TRIPOS>MOLECULE\n")==0){
                fgets(line_content, sizeof(line_content), fp);
                fgets(line_content, sizeof(line_content), fp);
                sscanf(line_content,"%d %d",&N_atoms, &N_bonds);
                i++;
                i++;
            }
            else if(strcmp(line_content, "@<TRIPOS>BOND\n")==0){
                line_Bonds = i;
            }
            else if(strcmp(line_content, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
            }
            i++;
    }

    //2) Scan the bond section and store bond information
    int bond_id[N_bonds], ori_atom_id[N_bonds], tar_atom_id[N_bonds];
    char bond_type[N_bonds][20];

    fclose(fp);

    fp=fopen(structure_filename,"r");
    while(fgets(line_content, sizeof(line_content), fp) != NULL){
        if(k==line_Bonds+1){
            for(i=0; i<N_bonds; i++){
                sscanf(line_content,"%d %d %d %s", &bond_id[i], &ori_atom_id[i], &tar_atom_id[i], &bond_type[i]);
                fgets(line_content, sizeof(line_content), fp);
                k++;
            }}
        else{
        k++;
        }
    }
    fclose(fp);

    //3) make structure overlay
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
    for(m=1; m<num_files+1; m++){
        char structure_filename_temp[100];
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, m);

        in=fopen(structure_filename_temp, "r");
        out=fopen(overlay_filename, "a");

        if(in==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename_temp);
            fclose(error_file);
        }

        k=1;
        char arguments[8][16];

        while(fgets(line_content, sizeof(line_content), in) != NULL){
            if(k==line_Atoms+1){
                for(i=1; i<N_atoms+1; i++){
                    sscanf(line_content,"%s %s %s %s %s %s %s %s",&arguments[0],&arguments[1],&arguments[2],&arguments[3],&arguments[4],&arguments[5],&arguments[6],&arguments[7]);
                    fprintf(out, "%d %s %s %s %s %s %s %s\n", atom_index,arguments[1],arguments[2],arguments[3],arguments[4],arguments[5],arguments[6],arguments[7]);
                    atom_index++;
                    fgets(line_content, sizeof(line_content), in);
                }
            k++;
            }
            else{
            k++;}
        }
        fclose(in);
        fclose(out);
    }
    atom_index=1;

    //Add the bond section using the bond information stored from earlier
    out=fopen(overlay_filename, "a");
    fprintf(out, "\n@<TRIPOS>BOND\n");

    for(m=0; m<num_files; m++){
        for(i=0; i<N_bonds; i++){
            fprintf(out, "%d %d %d %s\n", bond_id[i]+(m*N_bonds), ori_atom_id[i]+(m*N_atoms), tar_atom_id[i]+(m*N_atoms), bond_type[i]);
          //  printf("%d %d %d %s\n",bond_id[i], ori_atom_id[i], tar_atom_id[i], bond_type[i]);
        }
    }
    fclose(out);

}

void compile_all_mol2_files(char *base_filename, int num_primary_files, int num_other_files){
    //This function performs the same task as compile_mol2_files, however it will look for files also starting with other_
    //It also removes the said files after they have been read.
    char error_filename[128], other_base_filename[128];
    sprintf(error_filename,"Errors.txt");
    sprintf(other_base_filename,"other_%s_struct",base_filename);
    FILE *error_file;

    FILE *fp, *in, *out;
    int atom_index = 1, i = 1, j = 1, k = 1, line_Molecule, line_Atoms, line_Bonds;
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
    char line_content[120];

    while(fgets(line_content, sizeof(line_content), fp) != NULL){
            if(strcmp(line_content, "@<TRIPOS>MOLECULE\n")==0){
                fgets(line_content, sizeof(line_content), fp);
                fgets(line_content, sizeof(line_content), fp);
                sscanf(line_content,"%d %d",&N_atoms, &N_bonds);
                i++;
                i++;
            }
            else if(strcmp(line_content, "@<TRIPOS>BOND\n")==0){
                line_Bonds = i;
            }
            else if(strcmp(line_content, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
            }
            i++;
    }

    //2) Scan the bond section and store bond information
    int bond_id[N_bonds], ori_atom_id[N_bonds], tar_atom_id[N_bonds];
    char bond_type[N_bonds][20];

    fclose(fp);

    fp=fopen(structure_filename,"r");
    while(fgets(line_content, sizeof(line_content), fp) != NULL){
        if(k==line_Bonds+1){
            for(i=0; i<N_bonds; i++){
                sscanf(line_content,"%d %d %d %s", &bond_id[i], &ori_atom_id[i], &tar_atom_id[i], &bond_type[i]);
                fgets(line_content, sizeof(line_content), fp);
                k++;
            }}
        else{
        k++;
        }
    }
    fclose(fp);

    //3) make structure overlay
    char overlay_filename[128];

    sprintf(overlay_filename, "%s_overlay_all.mol2", base_filename);

    remove(overlay_filename);

    int m;
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
    for(m=1; m<num_primary_files+1; m++){
        char structure_filename_temp[100];
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, m);

        //add surface molecules to each structure before reading...used for .cif file
        sprintf(structure_filename_temp, "%s_struct%d.mol2", base_filename, m);

        in=fopen(structure_filename_temp, "r");
        out=fopen(overlay_filename, "a");

        if(in==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: File '%s' not found could not create structure overlay\n", structure_filename_temp);
            fclose(error_file);
        }

        k=1;
        char arguments[8][16];

        while(fgets(line_content, sizeof(line_content), in) != NULL){
            if(k==line_Atoms+1){
                for(i=1; i<N_atoms+1; i++){
                    sscanf(line_content,"%s %s %s %s %s %s %s %s",&arguments[0],&arguments[1],&arguments[2],&arguments[3],&arguments[4],&arguments[5],&arguments[6],&arguments[7]);
                    fprintf(out, "%d %s %s %s %s %s %s %s\n", atom_index,arguments[1],arguments[2],arguments[3],arguments[4],arguments[5],arguments[6],arguments[7]);
                    atom_index++;
                    fgets(line_content, sizeof(line_content), in);
                }
            k++;
            }
            else{
            k++;}
        }
        fclose(in);
        fclose(out);
        remove(structure_filename_temp);
    }

    //6) Add the atom section from each individual other file(Loop until m reaches the num_other_files)
    for(m=1; m<num_other_files+1; m++){
        char structure_filename_temp[100];
        sprintf(structure_filename_temp, "%s%d.mol2", other_base_filename, m);

        //add surface molecules to each structure before reading...used for .cif file
        add_surface(structure_filename_temp);
        sprintf(structure_filename_temp, "%s%d.mol2", other_base_filename, m);

        in=fopen(structure_filename_temp, "r");
        out=fopen(overlay_filename, "a");

        if(in==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: File '%s' not found, could not create structure overlay\n", structure_filename_temp);
            fclose(error_file);
        }

        k=1;
        char arguments[8][16];

        while(fgets(line_content, sizeof(line_content), in) != NULL){
            if(k==line_Atoms+1){
                for(i=1; i<N_atoms+1; i++){
                    sscanf(line_content,"%s %s %s %s %s %s %s %s",&arguments[0],&arguments[1],&arguments[2],&arguments[3],&arguments[4],&arguments[5],&arguments[6],&arguments[7]);
                    fprintf(out, "%d %s %s %s %s %s %s %s\n", atom_index,arguments[1],arguments[2],arguments[3],arguments[4],arguments[5],arguments[6],arguments[7]);
                    atom_index++;
                    fgets(line_content, sizeof(line_content), in);
                }
            k++;
            }
            else{
            k++;}
        }
        fclose(in);
        fclose(out);
        remove(structure_filename_temp);
    }
    atom_index=1;

    //7) Add the bond section using the bond information stored from earlier
    out=fopen(overlay_filename, "a");
    fprintf(out, "\n@<TRIPOS>BOND\n");

    for(m=0; m<total_num_files; m++){
        for(i=0; i<N_bonds; i++){
            fprintf(out, "%d %d %d %s\n", bond_id[i]+(m*N_bonds), ori_atom_id[i]+(m*N_atoms), tar_atom_id[i]+(m*N_atoms), bond_type[i]);
        }
    }
    fclose(out);
}
