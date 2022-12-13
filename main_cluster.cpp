#include "affected_atoms.hpp"
#include "alterations.hpp"
#include "overlay_structures.hpp"
#include "probability_ellipsoids.hpp"
#include <ctype.h>
/*Use this version of the main.cpp file when using a cluster. It uses an argument for the input file rather than asking for it in the command line*/
int main(int argc, char *argv[]){
    char input_filename[120], mol2_filename[120], error_filename[128], buffer[256], keyword[64], support[32];
    int  i, j, k, l, line_Atoms, line_Bonds, N_atoms=0, N_bonds=0, N_curves=0, N_constraints=0, meticulous=0, found_structures=0;
    long long int  N_steps_Z=1, N_steps_X=1, N_steps_Y=1, N_rotatable_bonds=0,max_acceptable_struct = 1000;
    double threshold_accuracy=90., z_min=0., z_max=0., cutoff_RMSD=2.5, minor_structures_CL;
    double surface_collision_distance = 1.5, interatomic_collision_distance = 1.5;
    vector<vector<int> > REDOR_det_index, REDOR_rec_index;
    vector< REDOR_dataset > REDOR;
    FILE *input, *mol2_file, *error_file, *log_file;
    log_file=fopen("log.txt","w");

    fprintf(log_file,"\n8888888 888b    888 88888888888 8888888888 8888888b.  8888888888     d8888  .d8888b.  8888888888  .d8888b.  \n");
    fprintf(log_file,"  888   8888b   888     888     888        888   Y88b 888           d88888 d88P  Y88b 888        d88P  Y88b \n");
    fprintf(log_file,"  888   88888b  888     888     888        888    888 888          d88P888 888    888 888        Y88b.      \n");
    fprintf(log_file,"  888   888Y88b 888     888     8888888    888   d88P 8888888     d88P 888 888        8888888     ^Y888b.   \n");
    fprintf(log_file,"  888   888 Y88b888     888     888        8888888P^  888        d88P  888 888        888            ^Y88b. \n");
    fprintf(log_file,"  888   888  Y88888     888     888        888 T88b   888       d88P   888 888    888 888              ^888 \n");
    fprintf(log_file,"  888   888   Y8888     888     888        888  T88b  888      d8888888888 Y88b  d88P 888        Y88b  d88P \n");
    fprintf(log_file,"8888888 888    Y888     888     8888888888 888   T88b 888     d88P     888  ^Y8888P^  8888888888  ^Y8888P^  \n");

    fprintf(log_file,"\n(Interpret NMR To Elucidate or Reconstruct the Full Atomistic Configurations of External Surfaces)\n");
    fprintf(log_file,"_____________________________________________________________________________________________________\n");
    fprintf(log_file,"\nA program for the automated structure elucidation of surface sites using RE(SP)DOR NMR, or other data\n");
    fprintf(log_file,"\nWritten by James Cunningham and Frederic A. Perras\n");
    fprintf(log_file,"US DOE, Ames National Laboratory, 2022\n");
    fprintf(log_file,"\nIf used for a publication, please cite: J. Magn. Reson. Open 2022, 12-13, 100066.\n");
    fprintf(log_file,"_____________________________________________________________________________________________________\n");

    printf("\n8888888 888b    888 88888888888 8888888888 8888888b.  8888888888     d8888  .d8888b.  8888888888  .d8888b.  \n");
    printf("  888   8888b   888     888     888        888   Y88b 888           d88888 d88P  Y88b 888        d88P  Y88b \n");
    printf("  888   88888b  888     888     888        888    888 888          d88P888 888    888 888        Y88b.      \n");
    printf("  888   888Y88b 888     888     8888888    888   d88P 8888888     d88P 888 888        8888888     ^Y888b.   \n");
    printf("  888   888 Y88b888     888     888        8888888P^  888        d88P  888 888        888            ^Y88b. \n");
    printf("  888   888  Y88888     888     888        888 T88b   888       d88P   888 888    888 888              ^888 \n");
    printf("  888   888   Y8888     888     888        888  T88b  888      d8888888888 Y88b  d88P 888        Y88b  d88P \n");
    printf("8888888 888    Y888     888     8888888888 888   T88b 888     d88P     888  ^Y8888P^  8888888888  ^Y8888P^  \n");

    printf("\n(Interpret NMR To Elucidate or Reconstruct the Full Atomistic Configurations of External Surfaces)\n");
    printf("_____________________________________________________________________________________________________\n");
    printf("\nA program for the automated structure elucidation of surface sites using RE(SP)DOR NMR, or other data\n");
    printf("\nWritten by James Cunningham and Frederic A. Perras\n");
    printf("US DOE, Ames National Laboratory, 2022\n");


    //Opening the error file to print out any issues that come up
    sprintf(error_filename,"Errors.txt");
    remove(error_filename);
    error_file=fopen(error_filename,"w");
    fclose(error_file);

  if (argc < 2) {
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Missing input file declaration\n");
        fclose(error_file);
        exit(1);
   }

   else{
   	printf("\nreading %s\n\n",argv[1]);
   	sprintf(input_filename,"%s",argv[1]);
   }

    //Here the provided input file is read
    input=fopen(input_filename,"r");
        if(input==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: Input file '%s' not found\n", input_filename);
            fclose(error_file);
            exit(1);
        }
    while(fgets(buffer, sizeof(buffer), input) != NULL){
        sscanf(buffer,"%s",keyword);
        if(strcmp(keyword, "meticulous")==0){
            meticulous=1;
            sprintf(keyword,"void");
        }

        else if(strcmp(keyword, "minor_structures_CL")==0){
            sscanf(buffer,"%s %lf",keyword, &minor_structures_CL);
            found_structures=1;
            meticulous=1; //the minor structures code uses the meticulous version of the REDOR calculations
            sprintf(keyword,"void");
        }

        else if(strcmp(keyword, "structure")==0){
            sscanf(buffer,"%s %s",keyword,mol2_filename);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "support")==0){
            sscanf(buffer,"%s %s",keyword,support);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "surface_collision_distance")==0){
            sscanf(buffer,"%s %lf",keyword,&surface_collision_distance);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "interatomic_collision_distance")==0){
            sscanf(buffer,"%s %lf",keyword,&interatomic_collision_distance);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "z_distance")==0){
            sscanf(buffer,"%s %lf %lf %d",keyword,&z_min, &z_max, &N_steps_Z);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "rotate_x")==0){
            sscanf(buffer,"%s %d",keyword, &N_steps_X);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "rotate_y")==0){
            sscanf(buffer,"%s %d",keyword, &N_steps_Y);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "confidence_level")==0){
            sscanf(buffer,"%s %lf",keyword,&threshold_accuracy);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "max_structures")==0){
            sscanf(buffer,"%s %d",keyword,&max_acceptable_struct);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "cutoff_rmsd")==0){
            sscanf(buffer,"%s %lf",keyword,&cutoff_RMSD);
            sprintf(keyword,"void");
        }
        //These calls are counted to determine the number of curves, bonds, and constraints, but they are not read yet.
        else if(strcmp(keyword, "surface-REDOR")==0){
            N_curves++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "intramolecular-REDOR")==0){
            N_curves++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "revolve")==0){
            N_rotatable_bonds++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "stretch")==0){
            N_rotatable_bonds++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "stretch_symmetric")==0){
            N_rotatable_bonds++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "bend")==0){
            N_rotatable_bonds++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "bend_symmetric")==0){
            N_rotatable_bonds++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "distance_constraint")==0 || strcmp(keyword, "angle_constraint")==0 || strcmp(keyword, "dihedral_constraint")==0 || strcmp(keyword, "surface_distance_constraint")==0){
            N_constraints++;
            sprintf(keyword,"void");
        }
    }
    fclose(input);

    //Knowing the quantity of these variables the arrays to contain their values are created and the file is read a second time to extract the values.
    Bond bond[N_rotatable_bonds];
    Constraint constraint[N_constraints];
    int bond_index=0, Nconst=0;
    input=fopen(input_filename,"r");
    REDOR.resize(N_curves);

    for(i=0;i<N_curves;i++){
        REDOR[i].order_parameter=1.;
    }

    int counter=0;
    i=j=0;
    while(fgets(buffer, sizeof(buffer), input) != NULL){
        sscanf(buffer," %s",keyword);
            if(strcmp(keyword, "revolve")==0){
                sscanf(buffer, "%s %d %d %d", keyword, &bond[bond_index].atom1, &bond[bond_index].atom2, &bond[bond_index].N_steps);
                bond[bond_index].atom1=bond[bond_index].atom1-1;
                bond[bond_index].atom2=bond[bond_index].atom2-1;
                bond[bond_index].type=0;
                bond_index++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "stretch")==0){
                sscanf(buffer, "%s %d %d %lf %lf %d", keyword, &bond[bond_index].atom1, &bond[bond_index].atom2, &bond[bond_index].dmin, &bond[bond_index].dmax, &bond[bond_index].N_steps);
                bond[bond_index].atom1=bond[bond_index].atom1-1;
                bond[bond_index].atom2=bond[bond_index].atom2-1;
                bond[bond_index].type=1;
                bond_index++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "stretch_symmetric")==0){
                sscanf(buffer, "%s %d %d %lf %lf %d", keyword, &bond[bond_index].atom1, &bond[bond_index].atom2, &bond[bond_index].dmin, &bond[bond_index].dmax, &bond[bond_index].N_steps);
                bond[bond_index].atom1=bond[bond_index].atom1-1;
                bond[bond_index].atom2=bond[bond_index].atom2-1;
                bond[bond_index].type=4;
                bond_index++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "bend")==0){
                sscanf(buffer, "%s %d %d %d %lf %lf %d", keyword, &bond[bond_index].atom0, &bond[bond_index].atom1, &bond[bond_index].atom2, &bond[bond_index].dmin, &bond[bond_index].dmax, &bond[bond_index].N_steps);
                bond[bond_index].atom0=bond[bond_index].atom0-1;
                bond[bond_index].atom1=bond[bond_index].atom1-1;
                bond[bond_index].atom2=bond[bond_index].atom2-1;
                bond[bond_index].type=2;
                bond_index++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "bend_symmetric")==0){
                sscanf(buffer, "%s %d %d %d %lf %lf %d", keyword, &bond[bond_index].atom0, &bond[bond_index].atom1, &bond[bond_index].atom2, &bond[bond_index].dmin, &bond[bond_index].dmax, &bond[bond_index].N_steps);
                bond[bond_index].atom0=bond[bond_index].atom0-1;
                bond[bond_index].atom1=bond[bond_index].atom1-1;
                bond[bond_index].atom2=bond[bond_index].atom2-1;
                bond[bond_index].type=3;
                bond_index++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "surface-REDOR")==0){
                sscanf(buffer,"%s %s %lf",keyword,REDOR[counter].filename, &REDOR[counter].scaling_factor);
                REDOR[counter].type=0;

                if((REDOR[counter].scaling_factor>1.)||(REDOR[counter].scaling_factor<=0.)){
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: Illegal REDOR curve scaling factor of %lf in %s, set it to default (1.0) \n", REDOR[counter].scaling_factor,REDOR[counter].filename);
                    fclose(error_file);
                    REDOR[counter].scaling_factor=1.0;
                }

                counter++;
                j=0;

                fgets(buffer, sizeof(buffer), input);
                sscanf(buffer,"%s",keyword);
                if(strcmp(keyword, "detected_spins")==0){
                    char *ptr = buffer, word[32];
                    sscanf(ptr,"%s",word);
                    ptr = strstr(ptr, word);
                    ptr += strlen(word); //skip the keyword
                    j=0;
                    while((sscanf(ptr,"%s",word)) == 1) {
                        REDOR[i].detected.push_back(0);
                        REDOR[i].recoupled.push_back(0);
                        REDOR[i].detected[j]=atoi(word)-1;
                        REDOR[i].recoupled[j]=0;
                        ptr = strstr(ptr, word);  // Find where the current word starts.
                        ptr += strlen(word); // Skip past the current word.
                        j++;
                    }

                }
                else{
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", REDOR[counter].filename);
                    fclose(error_file);
                    exit(1);
                }
                i++;
                sprintf(keyword,"void");
            }//surface_curve

            else if(strcmp(keyword, "intramolecular-REDOR")==0){
                sscanf(buffer,"%s %s %lf",keyword,&REDOR[counter].filename, &REDOR[counter].scaling_factor);
                REDOR[counter].type=1;

                if((REDOR[counter].scaling_factor>1.)||(REDOR[counter].scaling_factor<0.)){
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: Illegal REDOR curve scaling factor of %lf in %s, set it to default (1.0) \n", REDOR[counter].scaling_factor,REDOR[counter].filename);
                    fclose(error_file);
                    REDOR[counter].scaling_factor=1.0;
                }
                j=0;

                    fgets(buffer, sizeof(buffer), input);
                    sscanf(buffer,"%s",keyword);

                    if(strcmp(keyword, "detected_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR[i].detected.push_back(0);
                            REDOR[i].detected[j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }

                        fgets(buffer, sizeof(buffer), input);
                        sscanf(buffer,"%s",keyword);

                       if(strcmp(keyword, "recoupled_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR[i].recoupled.push_back(0);
                            REDOR[i].recoupled[j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }
                        if(REDOR[i].recoupled.size()>1){
                           error_file=fopen(error_filename,"a");
                           fprintf(error_file, "\nWARNING: calculated curve for %s might be inacurate at longer recoupling times, using root-sum-square dipole over recoupled spins.\n", REDOR[counter].filename);
                           fclose(error_file);
                        }

                    }
                       else{
                           error_file=fopen(error_filename,"a");
                           fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", REDOR[counter].filename);
                           fclose(error_file);
                           exit(1);
                        }
                    }

                    else if(strcmp(keyword, "recoupled_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR[i].recoupled.push_back(0);
                            REDOR[i].recoupled[j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }

                        fgets(buffer, sizeof(buffer), input);
                        sscanf(buffer,"%s",keyword);

                       if(strcmp(keyword, "detected_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR[i].detected.push_back(0);
                            REDOR[i].detected[j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }

                    }
                       else{
                           error_file=fopen(error_filename,"a");
                           fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", REDOR[counter].filename);
                           fclose(error_file);
                           exit(1);
                        }
                    }

                    else{
                        error_file=fopen(error_filename,"a");
                        fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", REDOR[counter].filename);
                        fclose(error_file);
                        exit(1);
                    }
                    i++;
                    counter++;
            }//intramolecular curve

            else if(strcmp(keyword, "distance_constraint")==0){
                sscanf(buffer,"%s %d %d %lf %lf",keyword,&constraint[Nconst].atom1, &constraint[Nconst].atom2, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].atom2--;
                constraint[Nconst].type=0;

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].maximum=value;
                }

                Nconst++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "angle_constraint")==0){
                sscanf(buffer,"%s %d %d %d %lf %lf",keyword,&constraint[Nconst].atom1, &constraint[Nconst].atom2, &constraint[Nconst].atom3, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].atom2--;
                constraint[Nconst].atom3--;
                constraint[Nconst].type=1;

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].maximum=value;
                }

                //Resetting the values so that they fall in the expected range of angles
                while(constraint[Nconst].minimum<0.){
                    constraint[Nconst].minimum = constraint[Nconst].minimum+180.;
                }
                while(constraint[Nconst].minimum>180.0){
                    constraint[Nconst].minimum = constraint[Nconst].minimum-180.;
                }
                while(constraint[Nconst].maximum<0.0){
                    constraint[Nconst].maximum = constraint[Nconst].maximum+180.;
                }
                while(constraint[Nconst].maximum>180.0){
                    constraint[Nconst].maximum = constraint[Nconst].maximum-180.;
                }

                Nconst++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "dihedral_constraint")==0){
                sscanf(buffer,"%s %d %d %d %d %lf %lf",keyword,&constraint[Nconst].atom1, &constraint[Nconst].atom2,&constraint[Nconst].atom3,&constraint[Nconst].atom4, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].atom2--;
                constraint[Nconst].atom3--;
                constraint[Nconst].atom4--;
                constraint[Nconst].type=2;

                //Resetting the values so that they fall in the expected range of dihedrals
                while(constraint[Nconst].minimum<-180.0){
                    constraint[Nconst].minimum = constraint[Nconst].minimum+360.;
                }
                while(constraint[Nconst].minimum>180.0){
                    constraint[Nconst].minimum = constraint[Nconst].minimum-360.;
                }
                while(constraint[Nconst].maximum<-180.0){
                    constraint[Nconst].maximum = constraint[Nconst].maximum+360.;
                }
                while(constraint[Nconst].maximum>180.0){
                    constraint[Nconst].maximum = constraint[Nconst].maximum-360.;
                }

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].maximum=value;
                }

                Nconst++;
                sprintf(keyword,"void");
            }

            else if(strcmp(keyword, "surface_distance_constraint")==0){
                sscanf(buffer,"%s %d %lf %lf",keyword,&constraint[Nconst].atom1, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].type=3;

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].maximum=value;
                }

                Nconst++;
                sprintf(keyword,"void");
            }

            else if(strcmp(keyword, "order_parameter")==0){
                int index;
                double S;
                sscanf(buffer,"%s %d %lf",keyword,&index, &S);

                if(index>N_curves){
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: Order parameter given for an inexistent REDOR curve\n");
                    fclose(error_file);
                }
                else if(S>1.0 || S<=0.0){
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: Illegal order parameter for curve %d; must be between 0 and 1.\n",index);
                    fclose(error_file);
                    exit(1);
                }
                else{
                    REDOR[index-1].order_parameter=S;
                }
                sprintf(keyword,"void");
            }
        }//end while
    fclose(input);

    //In the event that there are only 1 or 2 spins the fast approach is equivalent to the meticulous one.
    int metic[N_curves];
    for(i=0;i<N_curves;i++){
        metic[i]=0;
    }
    if(meticulous==1){
        for(i=0;i<N_curves;i++){
            if(REDOR[i].detected.size()>2)
                metic[i]=1;
            else
                metic[i]=0;
        }
    }

    //Next we read the provided starting mol2 file to extract the atomic coordinates and bonding connectivities.
    //The bonding connectivities are used when determining which atoms are affected by a given rotation.
    int mol2_len=strlen(mol2_filename);
    const char *filetype = &mol2_filename[mol2_len-5];
    if(strcmp(filetype,".mol2")!=0){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Structure must be provided as a *.mol2 file\n");
        fclose(error_file);
        exit(1);
    }

    mol2_file=fopen(mol2_filename, "r");

    if(mol2_file==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: mol2 file '%s' not found\n", mol2_filename);
        fclose(error_file);
        exit(1);
    }

    //First the program extracts the line numbers where the atom and bond tables begin, in addition to the number of atoms and bonds.
    i=j=0;
    while(fgets(buffer, sizeof(buffer), mol2_file) != NULL){
        sscanf(buffer,"%s",keyword);
            if(strcmp(keyword, "@<TRIPOS>MOLECULE")==0){
                fgets(buffer, sizeof(buffer), mol2_file);
                fgets(buffer, sizeof(buffer), mol2_file);
                sscanf(buffer,"%d %d",&N_atoms, &N_bonds);
                i++; j++;
                i++; j++;
            }
            if(strcmp(keyword, "@<TRIPOS>ATOM")==0){
                line_Atoms = i;
                }
            if(strcmp(keyword, "@<TRIPOS>BOND")==0){
                line_Bonds = j;
            }
            i++; j++;
            }
    fclose(mol2_file);

    //Now knowing how many bonds and atoms there are the following variables are created.
    //Atom variables
    vector< vector<double> > xyz;
    xyz.resize(N_atoms, vector<double>(3,0.));
    char element[N_atoms][3], atom_type[N_atoms][8];
    int atom_id[N_atoms];
    //bond variables
    int bond_id[N_bonds], ori_atom_id[N_bonds], tar_atom_id[N_bonds];
    char bond_type[N_bonds][20];
    vector< vector<int> > neighbors;
    neighbors.resize(N_atoms, vector<int> (1,0));
    for(i=0;i<N_atoms;i++){
        neighbors[i][0]=i;
    }

    //The mol2 file is opened for a second time to extract the coordinates and bond connections.
    k=1;
    mol2_file=fopen(mol2_filename,"r");
    while(fgets(buffer, sizeof(buffer), mol2_file) != NULL){
        if(k==line_Atoms+1){
            for(i=0; i<N_atoms; i++){
                fgets(buffer, sizeof(buffer), mol2_file);
                sscanf(buffer,"%d %s %lf %lf %lf %s", &atom_id[i], element[i],&xyz[i][0],&xyz[i][1],&xyz[i][2], atom_type[i]);
                //MOL2 files tend to mix the cases so this will normalize them to something like Si, not SI.
                element[i][0]=toupper(element[i][0]);
                element[i][1]=tolower(element[i][1]);
                k++;
            }}
        else if(k==line_Bonds){
            for(i=0; i<N_bonds; i++){
                fgets(buffer, sizeof(buffer), mol2_file);
                sscanf(buffer,"%d %d %d %s", &bond_id[i], &ori_atom_id[i], &tar_atom_id[i], bond_type[i]);
                neighbors[ori_atom_id[i]-1].push_back(tar_atom_id[i]-1);
                neighbors[tar_atom_id[i]-1].push_back(ori_atom_id[i]-1);
                k++;
            }}
        else{
            k++;
        }
    }
    fclose(mol2_file);

    center_structure(N_atoms,xyz);
    k = 0;

    //Here the program created tables of Chi^2 values as a function of distance and standard deviation of distance
    //between the atom and the surface plane.  There is one table per atom.
    //In the special case where there are two atoms contributing to a given REDOR curve, the program will calculate the
    //average REDOR curve for the pair, otherwise a gaussian distribution is used.
    vector< vector<double> > REDORs;
    for(i=0;i<REDOR.size();i++){
        REDOR[i].X2.resize(200,vector<double>(101,0.));
    }

    REDORs.resize(10000, vector<double>(9,0.));
    generate_REDORs(REDORs);
    printf("\n");
    fprintf(log_file,"\n");

    for(i=0; i<N_curves; i++){
        int Npoints=find_Npoints(REDOR[i].filename);
        REDOR[i].DSS0.resize(Npoints, 0.);
        REDOR[i].tmix.resize(Npoints, 0.);
        load_exp_curve(REDOR[i]);//populates DSS0 and tmix, the experimental data
        REDOR[i].DSS0_lib.resize(250, vector<double>(200,0.));
        for(j=0;j<N_atoms;j++){
            sprintf(REDOR[i].det_element,"%s",element[REDOR[i].detected[0]]);
            sprintf(REDOR[i].rec_element,"%s",element[REDOR[i].recoupled[0]]);
        }
        if(REDOR[i].type==0)
            load_simulations(support, REDOR[i]);
        else{
            REDOR[i].RDD1A=RDD_1A(REDOR[i].det_element, REDOR[i].rec_element);
            REDOR[i].spin=spin(REDOR[i].rec_element);
        }

        if(metic[i]==0){
            create_X2_table(REDOR[i],support);
            fprintf(log_file,"Creating Chi2 table for %s\n",REDOR[i].filename);
        }
        else{
            printf("Curve %s will be calculated on-the-fly\n",REDOR[i].filename);
            fprintf(log_file,"Curve %s will be calculated on-the-fly\n",REDOR[i].filename);
        }
    }

    //This function uses the bond list from the mol2 file to determine what atoms will be affected by
    //the rotation or elongation of a given bond.
    get_affected_atoms(N_rotatable_bonds,bond,neighbors);

    //After all the affected atoms were calculated we add the pairs with distance constraints to the lists
    //of neighbors so that collisions aren't triggered when working with things such as flexible rings
    //with short 1-bond distance constraints.
    for(i=0;i<N_constraints;i++){
        switch (constraint[i].type){
            case 0 : //distance
                neighbors[constraint[i].atom1].push_back(constraint[i].atom2);
                neighbors[constraint[i].atom2].push_back(constraint[i].atom1);
                break;
        }//end switch
    }

    k=1;
    int mol2_filename_len=strlen(mol2_filename)-5;
    char filename_base[mol2_filename_len+1];
    sprintf(filename_base,"%.*s",mol2_filename_len,mol2_filename);

    //We calculate the total number of structural iterations that will be done.
    //In order to make the program parallel it was necessary to have a single loop
    //over the various structural variables. To achieve this the program loops over
    //a generic iterator: it and uses modulo operations to determine this iteration's
    //bond angles, distances, etc. So each variable has an increment and a modulo (mod) variable
    long long int iterations=(long long int)1;
    for(i=0;i<N_rotatable_bonds;i++){
        bond[i].mod=iterations;
        iterations=iterations*(long long int)bond[i].N_steps;
    }

    long long int x_mod=iterations;
    iterations = iterations * (long long int)N_steps_X;
    double x_angle = 2*Pi/N_steps_X;
    if(N_steps_X==1)
        x_angle=0.;

    long long int y_mod=iterations;
    iterations = iterations * (long long int)N_steps_Y;
    double y_angle = 2*Pi/N_steps_Y;
    if(N_steps_Y==1)
        y_angle=0.;

    long long int z_mod=iterations;
    iterations = iterations * (long long int)N_steps_Z;
    double z_step = (z_max-z_min)/(N_steps_Z-1);
    if(N_steps_Z==1)
        z_step=0.;

    vector< vector<double> > xyz_ref;
    vector< vector<double> > xyz_best;
    xyz_ref.resize(N_atoms, vector<double>(3,0.));
    xyz_best.resize(N_atoms, vector<double>(3,0.));
    copy_structure(N_atoms,xyz,xyz_ref);

    //set chi2 values to a very high starting values to ensure a fit is found
    double chi2_min = 5000000;

    for(i=0; i<N_curves; i++){
        REDOR[i].chi2_min=5000000;
        REDOR[i].chi2_best=5000000;
    }

    vector< vector< vector< vector<double> > > > xyz_it;
    xyz_it.resize(N_curves, vector< vector< vector<double> > > (3, vector< vector<double> > (N_atoms, vector<double>(3,0.))));

    for(i=0;i<N_curves;i++){
        REDOR[i].d_index[1]=1000;
        REDOR[i].d_index[2]=0;
        REDOR[i].std_index[1]=1000;
        REDOR[i].std_index[2]=0;
    }

    //A loop begins below where all variables are sampled to find the best-fit structure
    //We also store the distance and std indices for this structure in order to supply its
    //simulated curves.
    printf("_____________________________________________________________________________________________________\n");
    printf("\nWill perform a search over a total of %d conformations\n",iterations);
    printf("Searching for the best-fit structure\n");
    printf("_____________________________________________________________________________________________________\n\n");

    fprintf(log_file,"_____________________________________________________________________________________________________\n");
    fprintf(log_file,"\nWill perform a search over a total of %d conformations\n",iterations);
    fprintf(log_file,"Searching for the best-fit structure\n");
    fprintf(log_file,"_____________________________________________________________________________________________________\n\n");
    int top_thread;

    #pragma omp parallel for
    for(long long int it=0;it<iterations;it++){
        //processor-specific variables
        vector< vector<double> > xyz_priv;
        xyz_priv.resize(N_atoms, vector<double>(3,0.));
        int ii, jj, kk, bond_position[N_rotatable_bonds];
        int d_indices[N_curves], std_indices[N_curves];
        double nominator, angle, R[4][3], bondvector[3], step;
        char min_chi2_output_filename[128];
        double curve_chi2[N_curves];
        copy_structure(N_atoms,xyz_ref,xyz_priv);

        //The calculations return the index for each of the structural variations:
        //z distance, y rotation, x rotation, bond rotations
        int z_position = int(floor(it/z_mod));
        int y_position = int(floor((it-z_position*z_mod)/y_mod));
        int x_position = int(floor((it-z_position*z_mod-y_position*y_mod)/x_mod));

        for(jj=N_rotatable_bonds-1;jj>=0;jj--){
            nominator=it-z_position*z_mod-y_position*y_mod-x_position*x_mod;
            for(kk=jj+1; kk<N_rotatable_bonds;kk++){
                nominator=nominator - bond_position[kk] * bond[kk].mod;
            }
            bond_position[jj] = int(floor(nominator/bond[jj].mod));
        }

        //Having found the indices for the different variations, the structural modifications are now done.
        translate_molecule_Z(N_atoms,xyz_priv,z_step*z_position+z_min);
        rotate_molecule_around_Y(N_atoms,xyz_priv,y_angle*y_position);
        rotate_molecule_around_X(N_atoms,xyz_priv,x_angle*x_position);

        for(jj=N_rotatable_bonds-1;jj>=0;jj--){
            if(bond[jj].type==0){//bond rotation
                angle=2.*Pi/bond[jj].N_steps*bond_position[jj];
                generate_bond_rot_matrix(R, xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2], angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                }
            }
            else if(bond[jj].type==1){//bond elongation
                get_internuclear_vector(bondvector,xyz_priv[bond[jj].atom1],xyz_priv[bond[jj].atom2]);
                step=bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin;
                bondvector[0]=bondvector[0]*step;
                bondvector[1]=bondvector[1]*step;
                bondvector[2]=bondvector[2]*step;
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    translate_atom_X(xyz_priv[bond[jj].affected_atom[ii]],bondvector[0]);
                    translate_atom_Y(xyz_priv[bond[jj].affected_atom[ii]],bondvector[1]);
                    translate_atom_Z(xyz_priv[bond[jj].affected_atom[ii]],bondvector[2]);
                }
            }
            else if(bond[jj].type==4){//bond elongation
                get_internuclear_vector(bondvector,xyz_priv[bond[jj].atom1],xyz_priv[bond[jj].atom2]);
                step=bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin;
                step=step/2.;
                bondvector[0]=bondvector[0]*step;
                bondvector[1]=bondvector[1]*step;
                bondvector[2]=bondvector[2]*step;
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    translate_atom_X(xyz_priv[bond[jj].affected_atom[ii]],bondvector[0]);
                    translate_atom_Y(xyz_priv[bond[jj].affected_atom[ii]],bondvector[1]);
                    translate_atom_Z(xyz_priv[bond[jj].affected_atom[ii]],bondvector[2]);
                }
                bondvector[0]=bondvector[0]*-1.;
                bondvector[1]=bondvector[1]*-1.;
                bondvector[2]=bondvector[2]*-1.;
                for(ii=0;ii<bond[jj].N_aff_atoms2; ii++){
                    translate_atom_X(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[0]);
                    translate_atom_Y(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[1]);
                    translate_atom_Z(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[2]);
                }
            }
            else if(bond[jj].type==2){//bond angle
                angle=(Pi/180.)*(bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin);
                generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                }
            }
            else{//bend_symmetric
                angle=(Pi/180.)*(bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin);
                angle=angle/2.;
                generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                }
                angle=-1.*angle;
                generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms2; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom2[ii]], R);
                }
            }
        }

        //Next we verify whether this structure satisfies all the specified constrains.  If not the processor moves on to the next structure.
        int check_constraints=0;
        double value;
        ii=0;
        for(ii=0;ii<N_constraints;ii++){
            switch (constraint[ii].type){
                case 0 : //distance
                    value = distance_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2]);
                    break;
                case 1: //angle
                    value = angle_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2],xyz_priv[constraint[ii].atom3]);
                    break;
                case 2: //dihedral
                    value = dihedral_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2],xyz_priv[constraint[ii].atom3],xyz_priv[constraint[ii].atom4]);
                    break;
                case 3: //surface distance
                    value=xyz_priv[constraint[ii].atom1][2];
                    break;

            }//end switch

            if(value > constraint[ii].maximum){
                check_constraints=1;
                break;
            }
            if(value < constraint[ii].minimum){
                check_constraints=1;
                break;
            }
        }

        //The collisions function returns true if 2 atoms are within 1 Angstrom from each other, or if
        //one atom falls at a distance below collision_distance from the surface.  A distance of exactly
        //zero is ignored as this is assumed to be a surface atom.
        if(check_constraints == 0){
        if(!collisions(xyz_priv, neighbors, surface_collision_distance, interatomic_collision_distance)){
            //Here we loop over the different curves to calculate the distance and STD index of that structure
            //in the Chi^2 table and take a sum to calculate the total chi 2
            double chi2=0.;

                for(kk=0; kk<N_curves; kk++){
                    d_indices[kk] = get_distance_index(REDOR[kk],xyz_priv);
                    std_indices[kk] = get_STDEV_index(REDOR[kk],xyz_priv);

                    if(metic[kk]==0){
                        curve_chi2[kk] = REDOR[kk].X2[d_indices[kk]][std_indices[kk]];
                     }
                    else
                        curve_chi2[kk] = calculate_curve_Chi2(REDOR[kk],xyz_priv,REDORs);

                    REDOR[kk].chi2_min = (curve_chi2[kk]<REDOR[kk].chi2_min)*curve_chi2[kk] + (curve_chi2[kk]>=REDOR[kk].chi2_min)*REDOR[kk].chi2_min;
                    chi2 = chi2 + curve_chi2[kk];
                }

            //If this structure has a lower Chi^2 value, then the shared variable chi2_min is replaced
            //We then also replace the top_thread variable to know which thread currently has the best structure
            if(chi2<chi2_min){
                chi2_min = chi2;
                top_thread=omp_get_thread_num();

                //We save the individual curve chi^2 values for this structure
                //These will be used to determine cutoff values later on.
                for(kk=0; kk<N_curves; kk++){
                    REDOR[kk].chi2_best = curve_chi2[kk];
                }

                //This structure overwrites xyz_best
                copy_structure(N_atoms,xyz_priv,xyz_best);
                printf("New Chi2 minimum at %lf\n",chi2_min);
                fprintf(log_file,"New Chi2 minimum at %lf\n",chi2_min);

                //Lastly, we save the distance and std indices for this structure, to write the best-fit curve
                for(kk=0; kk<N_curves; kk++){
                    REDOR[kk].d_index[0]=d_indices[kk];
                    REDOR[kk].std_index[0]=std_indices[kk];
                }

                //It is then saved as a mol2 file
                //each processor has its own "best" structure and these are compared at the end to determine
                //which one is the actual best-fit structure
                remove(min_chi2_output_filename);
                sprintf(min_chi2_output_filename, "%s_%d_best.mol2", filename_base,omp_get_thread_num());
                write_mol2(min_chi2_output_filename, N_bonds, atom_id, element, xyz_priv, atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);
            }

            //This this structure has the new lowest Chi^2 value for a given curve, the existing value is replaced.
        }
        else if(it==0){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nNOTE: The initial structure does not pass the collision distance constraints.\n");
            fclose(error_file);
        }
        }
        else if(it==0){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nNOTE: The initial structure does not pass the provided distance, angle, or dihedral constraints.\n");
            fclose(error_file);
        }
    }

    if (chi2_min == 5000000){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: No best-fit structure found with the given inputs\n");
        fclose(error_file);
        exit(1);
    }

    //We delete the individual thread best structures, leaving only the true best-fit structure mol2 file
    char best_filename[128];
    sprintf(best_filename,"%s_best.mol2", filename_base);
    remove(best_filename);

    #pragma omp parallel
    {
        char thread_filename[128];
        sprintf(thread_filename,"%s_%d_best.mol2", filename_base,omp_get_thread_num());

        if(omp_get_thread_num()==top_thread)
            rename(thread_filename,best_filename);

        remove(thread_filename);
    }

    for(i=0;i<N_curves;i++){
        copy_structure(N_atoms,xyz_best,xyz_it[i][0]);
    }

    //surface atoms are added to the mol2 file for plotting purposes.
    add_surface(best_filename);

    //the max_Chi2 function returns the highest allowable Chi^2 value given a minimum value and an error range
    //See numerical recipes in C for a description of the calculation.
    for(i=0; i<N_curves; i++){
        REDOR[i].chi2_max = max_Chi2(REDOR[i].chi2_min, threshold_accuracy/100.);
    }

    //If a cutoff value excludes the best fit structure, the cutoff is replaced with that from the best fit structure
    for(i=0; i<N_curves; i++){
        if(REDOR[i].chi2_max<REDOR[i].chi2_best){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nNOTE: Best-fit structure is not within experimental error for curve %d. Cutoff values were altered for this curve\n", i);
            fclose(error_file);
            REDOR[i].chi2_max = REDOR[i].chi2_best;
        }
    }

    //Next we find ALL the structures that fit within experimental error for all curves.
    //to prevent the writing of too many files a maximum number of structures is given
    //which will truncate the loop if too many structures are found

    int acceptable_structures = 0;
    int other_structures = 0;
    int max_accept_failsafe = 0;
    printf("_____________________________________________________________________________________________________\n");
    printf("\nFinding all structures that agree with experiment\n");
    printf("_____________________________________________________________________________________________________\n\n");

    fprintf(log_file,"_____________________________________________________________________________________________________\n");
    fprintf(log_file,"\nFinding all structures that agree with experiment\n");
    fprintf(log_file,"_____________________________________________________________________________________________________\n\n");

    #pragma omp parallel for
    for(long long int it=0;it<iterations;it++){
        if((acceptable_structures+other_structures)>max_acceptable_struct){
            //This is a safety feature that will limit the number of mol2 files that the program will save
            it = iterations;
            max_accept_failsafe = 1;
        }
        vector< vector<double> > xyz_priv;
        xyz_priv.resize(N_atoms, vector<double>(3,0.));
        int ii, jj, kk, bond_position[N_rotatable_bonds];
        int d_indices[N_curves], std_indices[N_curves];
        int check_chi2_threshold;
        double nominator, angle, deviation, R[4][3], bondvector[3], step;
        char min_chi2_output_filename[128];
        double curve_chi2[N_curves];
        copy_structure(N_atoms,xyz_ref,xyz_priv);

        int z_position = int(floor(it/z_mod));
        int y_position = int(floor((it-z_position*z_mod)/y_mod));
        int x_position = int(floor((it-z_position*z_mod-y_position*y_mod)/x_mod));

        for(jj=N_rotatable_bonds-1;jj>=0;jj--){
            nominator=it-z_position*z_mod-y_position*y_mod-x_position*x_mod;
            for(kk=jj+1; kk<N_rotatable_bonds;kk++){
                nominator=nominator - bond_position[kk] * bond[kk].mod;
            }
            bond_position[jj] = int(floor(nominator/bond[jj].mod));
        }

        translate_molecule_Z(N_atoms,xyz_priv,z_step*z_position+z_min);
        rotate_molecule_around_Y(N_atoms,xyz_priv,y_angle*y_position);
        rotate_molecule_around_X(N_atoms,xyz_priv,x_angle*x_position);

        for(jj=N_rotatable_bonds-1;jj>=0;jj--){
            if(bond[jj].type==0){//bond rotation
                angle=2.*Pi/bond[jj].N_steps*bond_position[jj];
                generate_bond_rot_matrix(R, xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2], angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                }
            }
            else if(bond[jj].type==1){//bond elongation
                get_internuclear_vector(bondvector,xyz_priv[bond[jj].atom1],xyz_priv[bond[jj].atom2]);
                step=bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin;
                bondvector[0]=bondvector[0]*step;
                bondvector[1]=bondvector[1]*step;
                bondvector[2]=bondvector[2]*step;
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    translate_atom_X(xyz_priv[bond[jj].affected_atom[ii]],bondvector[0]);
                    translate_atom_Y(xyz_priv[bond[jj].affected_atom[ii]],bondvector[1]);
                    translate_atom_Z(xyz_priv[bond[jj].affected_atom[ii]],bondvector[2]);
                }
            }
            else if(bond[jj].type==4){//bond elongation symmetric
                get_internuclear_vector(bondvector,xyz_priv[bond[jj].atom1],xyz_priv[bond[jj].atom2]);
                step=bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin;
                step=step/2.;
                bondvector[0]=bondvector[0]*step;
                bondvector[1]=bondvector[1]*step;
                bondvector[2]=bondvector[2]*step;
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    translate_atom_X(xyz_priv[bond[jj].affected_atom[ii]],bondvector[0]);
                    translate_atom_Y(xyz_priv[bond[jj].affected_atom[ii]],bondvector[1]);
                    translate_atom_Z(xyz_priv[bond[jj].affected_atom[ii]],bondvector[2]);
                }
                bondvector[0]=bondvector[0]*-1.;
                bondvector[1]=bondvector[1]*-1.;
                bondvector[2]=bondvector[2]*-1.;
                for(ii=0;ii<bond[jj].N_aff_atoms2; ii++){
                    translate_atom_X(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[0]);
                    translate_atom_Y(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[1]);
                    translate_atom_Z(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[2]);
                }
            }
            else if(bond[jj].type==2){//bond angle
                angle=(Pi/180.)*(bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin);
                generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                }
            }
            else{//bend_symmetric
                angle=(Pi/180.)*(bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin);
                angle=angle/2.;
                generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                }
                angle=-1.*angle;
                generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                //rotating all of atoms involved around the bond
                for(ii=0;ii<bond[jj].N_aff_atoms2; ii++){
                    rotate_around_bond2(xyz_priv[bond[jj].affected_atom2[ii]], R);
                }
            }
        }

        int check_constraints=0;
        double value;

        for(ii=0;ii<N_constraints;ii++){
            switch (constraint[ii].type){
                case 0 : //distance
                    value = distance_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2]);
                    break;
                case 1: //angle
                    value = angle_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2],xyz_priv[constraint[ii].atom3]);
                    break;
                case 2: //dihedral
                    value = dihedral_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2],xyz_priv[constraint[ii].atom3],xyz_priv[constraint[ii].atom4]);
                    break;
                case 3: //surface distance
                    value=xyz_priv[constraint[ii].atom1][2];
                    break;
            }//end switch

            if(value > constraint[ii].maximum){
                check_constraints=1;
                break;
            }
            if(value < constraint[ii].minimum){
                check_constraints=1;
                break;
            }
        }

        if(check_constraints == 0){
        if(!collisions(xyz_priv, neighbors, surface_collision_distance, interatomic_collision_distance)){
            check_chi2_threshold = 0;


            for(kk=0; kk<N_curves; kk++){
                d_indices[kk] = get_distance_index(REDOR[kk],xyz_priv);
                std_indices[kk] = get_STDEV_index(REDOR[kk],xyz_priv);

                if(metic[kk]==0)
                    curve_chi2[kk] = REDOR[kk].X2[d_indices[kk]][std_indices[kk]];
                else
                    curve_chi2[kk] = calculate_curve_Chi2(REDOR[kk],xyz_priv,REDORs);

                if(curve_chi2[kk]>REDOR[kk].chi2_max){
                    check_chi2_threshold++;
                    break;
                }
            }

            if(check_chi2_threshold == 0){
                //If a structure is found it is overlaid (reoriented) onto the best structure
                //The function also returns the RMSD between the two structures which can be used
                //to exclude certain structures from the overlay.
                deviation = overlay_structures(N_atoms, xyz_best, xyz_priv);

                //Each structure that is found is written either as blablabla_struct# or other_blablabla_struct#
                //These are later compiled into two overlay files.
                if(deviation<cutoff_RMSD){
                    acceptable_structures++;
                    sprintf(min_chi2_output_filename, "%s_struct%d.mol2", filename_base, acceptable_structures);
                    printf("(%d) Acceptable structure found\n", acceptable_structures);
                    fprintf(log_file,"(%d) Acceptable structure found\n", acceptable_structures);
                    write_mol2(min_chi2_output_filename, N_bonds, atom_id, element, xyz_priv, atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);

                    //The minimum and max distances are updated, if necessary, to print out the fitted curves ranges at the end.
                    for(kk=0; kk<N_curves; kk++){
                        if(d_indices[kk]<REDOR[kk].d_index[1]){
                            REDOR[kk].std_index[1] = std_indices[kk];
                            REDOR[kk].d_index[1] = d_indices[kk];
                            copy_structure(N_atoms,xyz_priv,xyz_it[kk][1]);
                        }
                        else if(d_indices[kk]>REDOR[kk].d_index[2]){
                            REDOR[kk].std_index[2] = std_indices[kk];
                            REDOR[kk].d_index[2] = d_indices[kk];
                            copy_structure(N_atoms,xyz_priv,xyz_it[kk][2]);
                        }

                        else if (d_indices[kk]==REDOR[kk].d_index[1]){
                            if(REDOR[kk].std_index[1]>std_indices[kk]){
                                REDOR[kk].std_index[1] = std_indices[kk];
                                copy_structure(N_atoms,xyz_priv,xyz_it[kk][1]);
                            }
                        }

                        else if (d_indices[kk]==REDOR[kk].d_index[2]){
                            if(REDOR[kk].std_index[2]<std_indices[kk]){
                                REDOR[kk].std_index[2] = std_indices[kk];
                                copy_structure(N_atoms,xyz_priv,xyz_it[kk][2]);
                            }
                        }
                    }
                }
                else{
                    other_structures++;
                    sprintf(min_chi2_output_filename, "other_%s_struct%d.mol2", filename_base, other_structures);
                    printf("(%d) Other acceptable structure found: rmsd with best = %lf\n", other_structures, deviation);
                    fprintf(log_file,"(%d) Other acceptable structure found: rmsd with best = %lf\n", other_structures, deviation);
                    write_mol2(min_chi2_output_filename, N_bonds, atom_id, element, xyz_priv, atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);
                }
            }
        }}
    }

    //Next we will write the structure overlay files
    char overlay_filename[128];
    sprintf(overlay_filename, "%s_overlay.mol2", filename_base);

    if((acceptable_structures+other_structures) <= 1){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Only one acceptable structure found within the given constraints\n");
        fclose(error_file);
        sprintf(overlay_filename, "%s_struct1.mol2", filename_base);
        remove(overlay_filename);
        exit(1);
    }

    printf("\nFound a total of %d acceptable structures,\n%d of which were within the requested RMSD from the best one\n", acceptable_structures+other_structures, acceptable_structures);
    fprintf(log_file,"\nFound a total of %d acceptable structures,\n%d of which were within the requested RMSD from the best one\n", acceptable_structures+other_structures, acceptable_structures);

    //if there are too many structures, erase the structures and exit the program
    if(max_accept_failsafe == 1){
        for(i=0;i<=acceptable_structures;i++){
            sprintf(overlay_filename, "%s_struct%d.mol2", filename_base,i);
            remove(overlay_filename);
        }
        for(i=0;i<=other_structures;i++){
            sprintf(overlay_filename, "other_%s_struct%d.mol2", filename_base,i);
            remove(overlay_filename);
        }
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Too many acceptable structures, exceeded limit of %d\n", max_acceptable_struct);
        fclose(error_file);
        return 0;
    }

    //Creating overlay and probability ellipsoids..
    printf("_____________________________________________________________________________________________________\n");
    printf("\nOverlaying structures\n");

    fprintf(log_file,"_____________________________________________________________________________________________________\n");
    fprintf(log_file,"\nOverlaying structures\n");

    //for primary structures
    compile_mol2_files(filename_base, acceptable_structures);
    printf("_____________________________________________________________________________________________________\n");
    fprintf(log_file,"_____________________________________________________________________________________________________\n");
    fclose(log_file);
    create_cif(overlay_filename, N_atoms+5);
    log_file=fopen("log.txt","a");

    //for all structures
    compile_all_mol2_files(filename_base, acceptable_structures, other_structures);

    //write out the fitted REDOR curve and ranges.
    printf("\n_____________________________________________________________________________________________________\n");
    printf("\nWriting the fitted RE(SP)DOR data to a file\n");

    fprintf(log_file,"\n_____________________________________________________________________________________________________\n");
    fprintf(log_file,"\nWriting the fitted RE(SP)DOR data to a file\n");
    if(meticulous==0)
        write_fits(filename_base,support,REDOR);
    else
        write_fits_meticulous(filename_base,support,REDOR,xyz_it);
    printf("_____________________________________________________________________________________________________\n");
    printf("\nStructure determination finished successfully\n");

    fprintf(log_file,"_____________________________________________________________________________________________________\n");
    fprintf(log_file,"\nStructure determination finished successfully\n");

    //Here starts the code for the minor structure search.
    //The code first calculates the dephasing levels form the best structure
    //It then performs a similar structure search in which calculated dephasing levels are the sum of those from the distinct species
    //For a minor species to be accepted, the calculated dephasing for each curve must be within error and the total dephasing much be
    //significantly low as to exclude the possibility of there being only the prior set of structures.
    if(found_structures>0){
        printf("_____________________________________________________________________________________________________\n");
        printf("\nBeginning a search for up to %d minor surface species\n", found_structures);

        fprintf(log_file,"_____________________________________________________________________________________________________\n");
        fprintf(log_file,"\nBeginning a search for up to %d minor surface species\n", found_structures);

        vector< vector< vector<double> > > xyz_minor;
        xyz_minor.resize(found_structures, vector< vector<double> > (N_atoms, vector<double>(3,0.)));
        vector<double> weights(1,1.);
        double new_weight;
        copy_structure(N_atoms,xyz_best,xyz_minor[0]);

        //This will be the total Chi2 threshold that the combined structure will have to beat
        double Chi2_limit, previous_Chi2=0., current_best_Chi2;
        for(i=0; i<N_curves; i++){
            previous_Chi2+= REDOR[i].chi2_best;
        }
        current_best_Chi2=previous_Chi2;
        Chi2_limit = max_Chi2_multi(previous_Chi2,minor_structures_CL/100.);

        for(i=0; i<N_curves; i++){
            int Npoints=REDOR[i].DSS0.size();
            REDOR[i].DSS0sim_prev.resize(Npoints, 0.);
            REDOR[i].DSS0sim_new.resize(Npoints, 0.);
            precalculate_dephasing(REDOR[i].DSS0sim_prev,REDOR[i],xyz_best,REDORs);
        }

        double current_CL=0.;
        do{//loop searching for minor surface species.
            xyz_minor.push_back(xyz_best);
            weights.push_back(0.);
            int found=0;
            new_weight=0.;
            printf("_____________________________________________________________________________________________________\n");
            printf("\nSearching for surface species number %d\n",found_structures+1);
            printf("Looking to reduce Chi2 from %lf to %lf\n",previous_Chi2,Chi2_limit);

            fprintf(log_file,"_____________________________________________________________________________________________________\n");
            fprintf(log_file,"\nSearching for surface species number %d\n",found_structures+1);
            fprintf(log_file,"Looking to reduce Chi2 from %lf to %lf\n",previous_Chi2,Chi2_limit);
            for(i=0;i<N_curves;i++){
                if(Chi2_limit>REDOR[i].chi2_max)
                    REDOR[i].chi2_max=Chi2_limit;
            }

            #pragma omp parallel for
            for(long long int it=0;it<iterations;it++){
                vector< vector<double> > xyz_priv;
                xyz_priv.resize(N_atoms, vector<double>(3,0.));
                int ii, jj, kk, bond_position[N_rotatable_bonds];
                int check_chi2_threshold;
                double nominator, angle, R[4][3], bondvector[3], step;
                double curve_chi2[N_curves][10];
                copy_structure(N_atoms,xyz_ref,xyz_priv);

                int z_position = int(floor(it/z_mod));
                int y_position = int(floor((it-z_position*z_mod)/y_mod));
                int x_position = int(floor((it-z_position*z_mod-y_position*y_mod)/x_mod));

                for(jj=N_rotatable_bonds-1;jj>=0;jj--){
                    nominator=it-z_position*z_mod-y_position*y_mod-x_position*x_mod;
                    for(kk=jj+1; kk<N_rotatable_bonds;kk++){
                        nominator=nominator - bond_position[kk] * bond[kk].mod;
                    }
                    bond_position[jj] = int(floor(nominator/bond[jj].mod));
                }

                translate_molecule_Z(N_atoms,xyz_priv,z_step*z_position+z_min);
                rotate_molecule_around_Y(N_atoms,xyz_priv,y_angle*y_position);
                rotate_molecule_around_X(N_atoms,xyz_priv,x_angle*x_position);

                for(jj=N_rotatable_bonds-1;jj>=0;jj--){
                    if(bond[jj].type==0){//bond rotation
                        angle=2.*Pi/bond[jj].N_steps*bond_position[jj];
                        generate_bond_rot_matrix(R, xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2], angle);
                        //rotating all of atoms involved around the bond
                        for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                            rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                        }
                    }
                    else if(bond[jj].type==1){//bond elongation
                        get_internuclear_vector(bondvector,xyz_priv[bond[jj].atom1],xyz_priv[bond[jj].atom2]);
                        step=bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin;
                        bondvector[0]=bondvector[0]*step;
                        bondvector[1]=bondvector[1]*step;
                        bondvector[2]=bondvector[2]*step;
                        for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                            translate_atom_X(xyz_priv[bond[jj].affected_atom[ii]],bondvector[0]);
                            translate_atom_Y(xyz_priv[bond[jj].affected_atom[ii]],bondvector[1]);
                            translate_atom_Z(xyz_priv[bond[jj].affected_atom[ii]],bondvector[2]);
                        }
                    }
                    else if(bond[jj].type==4){//bond elongation
                        get_internuclear_vector(bondvector,xyz_priv[bond[jj].atom1],xyz_priv[bond[jj].atom2]);
                        step=bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin;
                        step=step/2.;
                        bondvector[0]=bondvector[0]*step;
                        bondvector[1]=bondvector[1]*step;
                        bondvector[2]=bondvector[2]*step;
                        for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                            translate_atom_X(xyz_priv[bond[jj].affected_atom[ii]],bondvector[0]);
                            translate_atom_Y(xyz_priv[bond[jj].affected_atom[ii]],bondvector[1]);
                            translate_atom_Z(xyz_priv[bond[jj].affected_atom[ii]],bondvector[2]);
                        }
                        bondvector[0]=bondvector[0]*-1.;
                        bondvector[1]=bondvector[1]*-1.;
                        bondvector[2]=bondvector[2]*-1.;
                        for(ii=0;ii<bond[jj].N_aff_atoms2; ii++){
                            translate_atom_X(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[0]);
                            translate_atom_Y(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[1]);
                            translate_atom_Z(xyz_priv[bond[jj].affected_atom2[ii]],bondvector[2]);
                        }
                    }
                    else if(bond[jj].type==2){//bond angle
                        angle=(Pi/180.)*(bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin);
                        generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                        //rotating all of atoms involved around the bond
                        for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                            rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                        }
                    }
                    else{//bend_symmetric
                        angle=(Pi/180.)*(bond_position[jj]*(bond[jj].dmax-bond[jj].dmin)/(bond[jj].N_steps-1) + bond[jj].dmin);
                        angle=angle/2.;
                        generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                        //rotating all of atoms involved around the bond
                        for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                            rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
                        }
                        angle=-1.*angle;
                        generate_bond_angle_rot_matrix(R,xyz_priv[bond[jj].atom0], xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2],angle);
                        //rotating all of atoms involved around the bond
                        for(ii=0;ii<bond[jj].N_aff_atoms2; ii++){
                            rotate_around_bond2(xyz_priv[bond[jj].affected_atom2[ii]], R);
                        }
                    }
                }

                int check_constraints=0;
                double value;

                for(ii=0;ii<N_constraints;ii++){
                    switch (constraint[ii].type){
                        case 0 : //distance
                            value = distance_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2]);
                            break;
                        case 1: //angle
                            value = angle_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2],xyz_priv[constraint[ii].atom3]);
                            break;
                        case 2: //dihedral
                            value = dihedral_calc(xyz_priv[constraint[ii].atom1],xyz_priv[constraint[ii].atom2],xyz_priv[constraint[ii].atom3],xyz_priv[constraint[ii].atom4]);
                            break;
                        case 3: //surface distance
                            value=xyz_priv[constraint[ii].atom1][2];
                            break;
                    }//end switch

                    if(value > constraint[ii].maximum){
                        check_constraints=1;
                        break;
                    }
                    if(value < constraint[ii].minimum){
                        check_constraints=1;
                        break;
                    }
                }

                if(check_constraints == 0){
                if(!collisions(xyz_priv, neighbors, surface_collision_distance, interatomic_collision_distance)){
                    check_chi2_threshold = 0;

                    double chi2=0;
                    for(kk=0; kk<N_curves; kk++){
                        if(calculate_curve_Chi2_multi(curve_chi2[kk],REDOR[kk],xyz_priv,REDORs)>REDOR[kk].chi2_max){
                            check_chi2_threshold++;
                            break;
                        }
                    }

                    if(check_chi2_threshold == 0){
                        for(jj=0;jj<10;jj++){
                            chi2=0.;
                            check_chi2_threshold=0;
                            for(kk=0;kk<N_curves;kk++){
                                if(curve_chi2[kk][jj]>REDOR[kk].chi2_max){
                                    check_chi2_threshold++;
                                    break;
                                }
                                chi2 = chi2 + curve_chi2[kk][jj];
                            }
                            if(check_chi2_threshold == 0){
                                if(chi2<current_best_Chi2){
                                    current_best_Chi2=chi2;
                                    current_CL=return_CI(previous_Chi2,chi2)*100;
                                    printf("\nChi2 reduced to %lf; corresponding to a confidence level of %.1lf percent",chi2,current_CL);
                                    fprintf(log_file,"\nChi2 reduced to %lf; corresponding to a confidence level of %.1lf percent",chi2,current_CL);
                                }
                                if(chi2<Chi2_limit){
                                    Chi2_limit=chi2;
                                    copy_structure(N_atoms,xyz_priv,xyz_minor[found_structures]);
                                    found=1;
                                    new_weight=0.05+0.05*jj;
                                    printf("  (secondary structure found!)");
                                    fprintf(log_file,"  (secondary structure found!)");
                                }
                            }
                            else
                                check_chi2_threshold=0;
                        }
                    }
                }}
            }

            if(found !=0){
                for(i=0;i<found_structures;i++){
                    double deviation = overlay_structures(N_atoms, xyz_minor[i], xyz_minor[found_structures]);
                    if(deviation<cutoff_RMSD){
                        printf("\nMinor structure was found to be within error of the main structure.");
                        fprintf(log_file,"\nMinor structure was found to be within error of the main structure.");
                        found=0;
                        break;
            }}}

            found_structures++;
            weights[found_structures-1]=new_weight;
            for(i=0;i<found_structures-1;i++){
                weights[i]=weights[i]*(1.-new_weight);
            }

            if(found==0)
                break;

            for(i=0; i<N_curves; i++){
                int Npoints=find_Npoints(REDOR[i].filename);
                precalculate_dephasing(REDOR[i].DSS0sim_new,REDOR[i],xyz_minor[found_structures-1],REDORs);
                for(j=0;j<Npoints;j++){
                    REDOR[i].DSS0sim_prev[j]=(1.-new_weight)*REDOR[i].DSS0sim_prev[j] + new_weight*REDOR[i].DSS0sim_new[j];
                }
            }
            previous_Chi2=Chi2_limit;
            Chi2_limit=max_Chi2_multi(Chi2_limit,minor_structures_CL/100.);
            printf("\nIdentified surface species number %d\n", found_structures);
            fprintf(log_file,"\nIdentified surface species number %d\n", found_structures);
        }while(true);

        for(i=0;i<found_structures;i++){
            char conformer_filename[128];
            sprintf(conformer_filename, "%s_conformer_%d_%.2lf.mol2",filename_base,i+1,weights[i]);
            write_mol2(conformer_filename, N_bonds, atom_id, element, xyz_minor[i], atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);
            add_surface(conformer_filename);
        }
        printf("\n_____________________________________________________________________________________________________\n");
        printf("\nThe search for minor surface species has completed.\n");
        printf("_____________________________________________________________________________________________________\n");

        fprintf(log_file,"\n_____________________________________________________________________________________________________\n");
        fprintf(log_file,"\nThe search for minor surface species has completed.\n");
        fprintf(log_file,"_____________________________________________________________________________________________________\n");
        if(found_structures>1){
            printf("\nINTERFACES was able to improve the quality of the (%.1lf percent confidence interval)\n",current_CL);
            printf("REDOR fits with the inclusion of %d conformers. \n",found_structures);
            printf("These structures have been saved, inspect them to ensure that they are truly different.\n");
            printf("\nThe following are the weights of the %d conformers:",found_structures);

            fprintf(log_file,"\nINTERFACES was able to improve the quality of the (%.1lf percent confidence interval)\n",current_CL);
            fprintf(log_file,"REDOR fits with the inclusion of %d conformers. \n",found_structures);
            fprintf(log_file,"These structures have been saved, inspect them to ensure that they are truly different.\n");
            fprintf(log_file,"\nThe following are the weights of the %d conformers:",found_structures);
            for(i=0;i<found_structures;i++){
                printf("\nConformer %d has a weight of %.2lf",i+1,weights[i]);
                fprintf(log_file,"\nConformer %d has a weight of %.2lf",i+1,weights[i]);
            }

            //write out the fitted REDOR curve and ranges.
            printf("\n_____________________________________________________________________________________________________\n");
            printf("\nWriting the fitted RE(SP)DOR data to a file\n");
            fprintf(log_file,"\n_____________________________________________________________________________________________________\n");
            fprintf(log_file,"\nWriting the fitted RE(SP)DOR data to a file\n");
            write_fits_multi(found_structures,filename_base,REDOR,xyz_minor,weights);
        }
        else{
            printf("\nINTERFACES couldn't identify minor surface species within the specified %.1lf percent confidence interval\n",minor_structures_CL);
            fprintf(log_file,"\nINTERFACES couldn't identify minor surface species within the specified %.1lf percent confidence interval\n",minor_structures_CL);
        }

        printf("_____________________________________________________________________________________________________\n");
        printf("\nMinor structure determination finished successfully\n");
        printf("_____________________________________________________________________________________________________\n");

        fprintf(log_file,"_____________________________________________________________________________________________________\n");
        fprintf(log_file,"\nMinor structure determination finished successfully\n");
        fprintf(log_file,"_____________________________________________________________________________________________________\n");
    }//end of minor structure search

    fclose(log_file);
    return 0;
}//end int main
