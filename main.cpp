#include "affected_atoms.hpp"
#include "alterations.hpp"
#include "overlay_structures.hpp"
#include "probability_ellipsoids.hpp"

int main(){
    char input_filename[120], mol2_filename[120],  output_filename[128], error_filename[128], buffer[256], keyword[64], support[32];
    int I, i = 0, j = 0, line_Atoms, line_Bonds, N_atoms=0, N_bonds=0, N_curves=0, N_constraints = 0;
    int  N_steps_Z=1, N_steps_X=1, N_steps_Y=1, N_rotatable_bonds=0,max_acceptable_struct = 1000;
    double threshold_accuracy=90., z_min=0., z_max=0., cutoff_RMSD=2.5;
    double collision_distance = 1.5;
    vector<vector<int>> REDOR_det_index, REDOR_rec_index;
    FILE *input, *mol2_file, *out, *error_file;

    printf("**********************************************************\n");
    printf("**********************************************************\n");
    printf("***************Surface 3D Structure Solver****************\n");
    printf("****written by James Cunningham and Frederic A. Perras****\n");
    printf("******************Ames Laboratory, 2022*******************\n");
    printf("**********************************************************\n");
    printf("**********************************************************\n");


    printf ("\nEnter the name of the input file you wish to open.\n");
    scanf ("%s", &input_filename);

    //Opening the error file to print out any issues that come up
    sprintf(error_filename,"Errors.txt");
    remove(error_filename);
    error_file=fopen(error_filename,"w");
    fclose(error_file);

    //Here the provided input file is read
    input=fopen(input_filename,"r");
        if(input==NULL){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: Input file '%s' not found\n", input_filename);
            fclose(error_file);
            exit(1);
        }
    while(fgets(buffer, sizeof(buffer), input) != NULL){
        sscanf(buffer,"%s",&keyword);
        if(strcmp(keyword, "structure")==0){
            sscanf(buffer,"%s %s",&keyword,&mol2_filename);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "support")==0){
            sscanf(buffer,"%s %s",&keyword,&support);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "collision_distance")==0){
            sscanf(buffer,"%s %lf",&keyword,&collision_distance);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "z_distance")==0){
            sscanf(buffer,"%s %lf %lf %d",&keyword,&z_min, &z_max, &N_steps_Z);
            sprintf(keyword,"void");
        }

        else if(strcmp(keyword, "rotate_x")==0){
            sscanf(buffer,"%s %d",&keyword, &N_steps_X);
            sprintf(keyword,"void");
        }

        else if(strcmp(keyword, "rotate_y")==0){
            sscanf(buffer,"%s %d",&keyword, &N_steps_Y);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "confidence_level")==0){
            sscanf(buffer,"%s %lf",&keyword,&threshold_accuracy);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "max_structures")==0){
            sscanf(buffer,"%s %d",&keyword,&max_acceptable_struct);
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "cutoff_rmsd")==0){
            sscanf(buffer,"%s %lf",&keyword,&cutoff_RMSD);
            sprintf(keyword,"void");
        }
        //These calls are counted to determine the number of curves, bonds, and constraints, but they are not read yet.
        else if(strcmp(keyword, "surface-curve")==0){
            N_curves++;
            sprintf(keyword,"void");
        }

        else if(strcmp(keyword, "intramolecular-curve")==0){
            N_curves++;
            sprintf(keyword,"void");
        }

        else if(strcmp(keyword, "bond")==0){
            N_rotatable_bonds++;
            sprintf(keyword,"void");
        }
        else if(strcmp(keyword, "distance_constraint")==0 || strcmp(keyword, "angle_constraint")==0 || strcmp(keyword, "dihedral_constraint")==0){
            N_constraints++;
            sprintf(keyword,"void");
        }
    }
    fclose(input);
    i=0;
    j=0;
    int k=0;

    //Knowing the quantity of these variables the arrays to contain their values are created and the file is read a second time to extract the values.
    struct Bond bond[N_rotatable_bonds];
    struct Constraint constraint[N_constraints];
    int bond_index=0, Nconst=0, Nspins[N_curves], curve_type[N_curves], Nrecspins[N_curves];
    char curve_filename[N_curves][120];
    double scaling_factor[N_curves];
    N_curves = 0;

    input=fopen(input_filename,"r");

    while(fgets(buffer, sizeof(buffer), input) != NULL){
        sscanf(buffer," %s",&keyword);
            if(strcmp(keyword, "bond")==0){
                sscanf(buffer, "%s %d %d %d", &keyword, &bond[bond_index].atom1, &bond[bond_index].atom2, &bond[bond_index].N_steps);
                bond[bond_index].atom1=bond[bond_index].atom1-1;
                bond[bond_index].atom2=bond[bond_index].atom2-1;
                bond_index++;
                sprintf(keyword,"void");
            }


            else if(strcmp(keyword, "surface-curve")==0){
                REDOR_det_index.push_back(vector<int>());
                REDOR_rec_index.push_back(vector<int>());
                sscanf(buffer,"%s %s %lf",&keyword,&curve_filename[N_curves], &scaling_factor[N_curves]);
                curve_type[N_curves]=0;

                if((scaling_factor[N_curves]>1.)||(scaling_factor[N_curves]<=0.)){
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: Illegal REDOR curve scaling factor of %d in %s, set it to default (1.0) \n", scaling_factor[N_curves],curve_filename[N_curves]);
                    fclose(error_file);
                    scaling_factor[N_curves]=1.0;
                }

                N_curves++;
                j=0;

                fgets(buffer, sizeof(buffer), input);
                sscanf(buffer,"%s",&keyword);
                if(strcmp(keyword, "detected_spins")==0){
                    char *ptr = buffer, word[32];
                    sscanf(ptr,"%s",word);
                    ptr = strstr(ptr, word);
                    ptr += strlen(word); //skip the keyword
                    j=0;
                    while((sscanf(ptr,"%s",word)) == 1) {
                        REDOR_det_index[i].push_back(0);
                        REDOR_rec_index[i].push_back(0);
                        REDOR_det_index[i][j]=atoi(word)-1;
                        REDOR_rec_index[i][j]=0;
                        ptr = strstr(ptr, word);  // Find where the current word starts.
                        ptr += strlen(word); // Skip past the current word.
                        j++;
                    }
                    Nspins[i]=j;

                }
                else{
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", curve_filename[N_curves]);
                    fclose(error_file);
                    exit(1);
                }
                i++;
                sprintf(keyword,"void");
            }//surface_curve

            else if(strcmp(keyword, "intramolecular-curve")==0){
                REDOR_det_index.push_back(vector<int>());
                REDOR_rec_index.push_back(vector<int>());
                sscanf(buffer,"%s %s %lf",&keyword,&curve_filename[N_curves], &scaling_factor[N_curves]);
                curve_type[N_curves]=1;

                if((scaling_factor[N_curves]>1.)||(scaling_factor[N_curves]<0.)){
                    error_file=fopen(error_filename,"a");
                    fprintf(error_file, "\nERROR: Illegal REDOR curve scaling factor of %d in %s, set it to default (1.0) \n", scaling_factor[N_curves],curve_filename[N_curves]);
                    fclose(error_file);
                    scaling_factor[N_curves]=1.0;
                }

                N_curves++;
                j=0;

                    fgets(buffer, sizeof(buffer), input);
                    sscanf(buffer,"%s",&keyword);

                    if(strcmp(keyword, "detected_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR_det_index[i].push_back(0);
                            REDOR_det_index[i][j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }
                        Nspins[i]=j;

                        fgets(buffer, sizeof(buffer), input);
                        sscanf(buffer,"%s",&keyword);

                       if(strcmp(keyword, "recoupled_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR_rec_index[i].push_back(0);
                            REDOR_rec_index[i][j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }
                        Nrecspins[i]=j;
                    }
                       else{
                           error_file=fopen(error_filename,"a");
                           fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", curve_filename[N_curves]);
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
                            REDOR_rec_index[i].push_back(0);
                            REDOR_rec_index[i][j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }
                        Nrecspins[i]=j;

                        fgets(buffer, sizeof(buffer), input);
                        sscanf(buffer,"%s",&keyword);

                       if(strcmp(keyword, "detected_spins")==0){
                       char *ptr = buffer, word[32];
                       sscanf(ptr,"%s",word);
                       ptr = strstr(ptr, word);
                       ptr += strlen(word); //skip the keyword
                       j=0;
                        while((sscanf(ptr,"%s",word)) == 1) {
                            REDOR_det_index[i].push_back(0);
                            REDOR_det_index[i][j]=atoi(word)-1;
                            ptr = strstr(ptr, word);  // Find where the current word starts.
                            ptr += strlen(word); // Skip past the current word.
                            j++;
                        }
                        Nspins[i]=j;

                    }
                       else{
                           error_file=fopen(error_filename,"a");
                           fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", curve_filename[N_curves]);
                           fclose(error_file);
                           exit(1);
                        }
                    }

                    else{
                        error_file=fopen(error_filename,"a");
                        fprintf(error_file, "\nERROR: missing coupled spins for %s, exiting.\n", curve_filename[N_curves]);
                        fclose(error_file);
                        exit(1);
                    }
                    i++;
            }//intramolecular curve


            else if(strcmp(keyword, "distance_constraint")==0){
                sscanf(buffer,"%s %d %d %lf %lf",&keyword,&constraint[Nconst].atom1, &constraint[Nconst].atom2, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].atom2--;
                constraint[Nconst].type=0;

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].minimum=value;
                }

                Nconst++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "angle_constraint")==0){
                sscanf(buffer,"%s %d %d %d %lf %lf",&keyword,&constraint[Nconst].atom1, &constraint[Nconst].atom2, &constraint[Nconst].atom3, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].atom2--;
                constraint[Nconst].atom3--;
                constraint[Nconst].type=1;

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].minimum=value;
                }

                Nconst++;
                sprintf(keyword,"void");
            }
            else if(strcmp(keyword, "dihedral_constraint")==0){
                sscanf(buffer,"%s %d %d %d %d %lf %lf",&keyword,&constraint[Nconst].atom1, &constraint[Nconst].atom2,&constraint[Nconst].atom3,&constraint[Nconst].atom4, &constraint[Nconst].minimum, &constraint[Nconst].maximum);
                constraint[Nconst].atom1--;
                constraint[Nconst].atom2--;
                constraint[Nconst].atom3--;
                constraint[Nconst].atom4--;
                constraint[Nconst].type=2;

                if(constraint[Nconst].minimum>constraint[Nconst].maximum){
                    double value=constraint[Nconst].minimum;
                    constraint[Nconst].minimum=constraint[Nconst].maximum;
                    constraint[Nconst].minimum=value;
                }

                Nconst++;
                sprintf(keyword,"void");
            }
        }//end while
    fclose(input);

    //Next we read the provided starting mol2 file to extract the atomic coordinates and bonding connectivities.
    //The bonding connectivities are used when determining which atoms are affected by a given rotation.
    mol2_file=fopen(mol2_filename, "r");

    if(mol2_file==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: mol2 file '%s' not found\n", mol2_filename);
        fclose(error_file);
        exit(1);
    }

    i=0;
    j=0;

    //First the program extracts the line numbers where the atom and bond tables begin, in addition to the number of atoms and bonds.
    while(fgets(buffer, sizeof(buffer), mol2_file) != NULL){
        sscanf(buffer,"%s",&keyword);
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
    double xyz[N_atoms][3];
    char element[N_atoms][3], atom_type[N_atoms][8];
    int atom_id[N_atoms];
    //bond variables
    int bond_id[N_bonds], ori_atom_id[N_bonds], tar_atom_id[N_bonds];
    char bond_type[N_bonds][20];

    k = 1;

    //The mol2 file is opened for a second time to extract the coordinates and bond connections.
    mol2_file=fopen(mol2_filename,"r");
    while(fgets(buffer, sizeof(buffer), mol2_file) != NULL){
        if(k==line_Atoms+1){
            for(i=0; i<N_atoms; i++){
                fgets(buffer, sizeof(buffer), mol2_file);
                sscanf(buffer,"%d %s %lf %lf %lf %s", &atom_id[i], &element[i],&xyz[i][0],&xyz[i][1],&xyz[i][2], atom_type[i]);
                k++;
            }}
        else if(k==line_Bonds){
            for(i=0; i<N_bonds; i++){
                fgets(buffer, sizeof(buffer), mol2_file);
                sscanf(buffer,"%d %d %d %s", &bond_id[i], &ori_atom_id[i], &tar_atom_id[i], &bond_type[i]);
                k++;
            }}
        else{
            k++;
        }
    }
    fclose(mol2_file);

    center_structure(N_atoms,xyz);

    double NMR_z_coord_sum[N_curves], NMR_z_avg[N_curves];
    k = 0;

    //Here the program created tables of Chi^2 values as a function of distance and standard deviation of distance
    //between the atom and the surface plane.  There is one table per atom.
    vector<vector<vector<double>>> X2;
    X2.resize(N_curves, vector<vector<double>>(200,vector<double>(101,0.)));
    for(i=0; i<N_curves; i++){
        create_X2_table(curve_filename[i], support, element[REDOR_det_index[i][0]],element[REDOR_rec_index[i][0]], X2[i],scaling_factor[i], Nspins[i], curve_type[i]);
    }


    //This function will use the bond list from the mol2 file to determine what atoms will be affected by
    //the rotation about a given bond.
    get_aff_atoms(N_bonds, N_rotatable_bonds, bond, ori_atom_id, tar_atom_id);

    k=1;
    int fnamer_len=strlen(mol2_filename)-5;
    char base_iteration_filename[fnamer_len+1];
    sprintf(base_iteration_filename,"%.*s",fnamer_len,mol2_filename);

    //We calculate the total number of structural iterations that will be done.
    //In order to make the program parallel it was necessary to have a single loop
    //over the various structural variables. To achieve this the program loops over
    //a generic iterator: it and uses modulo operations to determine this iteration's
    //bond angles, distances, etc. So each variable has an increment and a modulo (mod) variable
    int iterations=1;
    for(i=0;i<N_rotatable_bonds;i++){
        bond[i].mod=iterations;
        iterations=iterations*bond[i].N_steps;
    }

    int x_mod=iterations;
    iterations = iterations * N_steps_X;
    double x_angle = 2*Pi/N_steps_X;

    int y_mod=iterations;
    iterations = iterations * N_steps_Y;
    double y_angle = 2*Pi/N_steps_Y;

    int z_mod=iterations;
    iterations = iterations * N_steps_Z;
    double z_step = (z_max-z_min)/(N_steps_Z-1);
    if(N_steps_Z==1)
        z_step=0.;

    double xyz_ref[N_atoms][3];
    copy_structure(N_atoms,xyz,xyz_ref);
    double xyz_best[N_atoms][3];

    //set chi2 values to a very high starting values to ensure a fit is found
    double chi2_min = 5000000;
    double curve_chi2[N_curves], curve_chi2_min[N_curves], curve_chi2_max[N_curves], best_struct_curve_chi2[N_curves];
    int d_indices_min[N_curves], std_indices_min[N_curves];

    for(i=0; i<N_curves; i++){
        curve_chi2_min[i]=5000000;
        best_struct_curve_chi2[i]=5000000;
        d_indices_min[i]=5000000;
        std_indices_min[i]=5000000;

    }

    //These arrays to store the average  and stdev distances of the (0)best fit structure, (1)Smallest distance and (2)Largest distance from the surface
    double distance_BSL[N_curves][3], stdev_BSL[N_curves][3];
    int d_indices_range[N_curves][3], std_indices_range[N_curves][3];
    for(i=0;i<N_curves;i++){
        d_indices_range[i][1]=1000;
        d_indices_range[i][2]=0;
        std_indices_range[i][1]=1000;
        std_indices_range[i][2]=0;
    }

    //A loop begins below where all variables are sampled to find the best-fit structure
    //We also store the distance and std indices for this structure in order to supply its
    //simulated curves.
    printf("\nWill perform a search over a total of %d conformations\n",iterations);

    printf("\nSearching for the best-fit structure\n");
    int top_thread;

    #pragma omp parallel for
    for(int it=0;it<iterations;it++){
        //processor-specific variables

        int ii, jj, kk, bond_position[N_rotatable_bonds];
        int d_indices[N_curves], std_indices[N_curves];
        double xyz_priv[N_atoms][3], nominator, angle, R[4][3];
        char min_chi2_output_filename[128];
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
            angle=2*Pi/bond[jj].N_steps*bond_position[jj];
            generate_bond_rot_matrix(R, xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2], angle);
            //rotating all of atoms involved around the bond
            for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
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
        if(!collisions(xyz_priv, N_atoms, collision_distance)){
        if(check_constraints == 0){

            //Here we loop over the different curves to calculate the distance and STD index of that structure
            //in the Chi^2 table
            for(kk=0; kk<N_curves; kk++){
                d_indices[kk] = get_distance_index(REDOR_det_index,REDOR_rec_index, xyz_priv,kk,curve_type[kk]);
                std_indices[kk] = get_STDEV_index(REDOR_det_index,REDOR_rec_index, xyz_priv, kk,curve_type[kk]);
                curve_chi2[kk] = X2[kk][d_indices[kk]][std_indices[kk]];
                curve_chi2_min[kk] = (curve_chi2[kk]<curve_chi2_min[kk])*curve_chi2[kk] + (curve_chi2[kk]>=curve_chi2_min[kk])*curve_chi2_min[kk] ;
            }
            //We then calculate the total Chi^2 of that structure
            double chi2 = calc_chi2(X2, std_indices, d_indices, N_curves);

            //If this structure has a lower Chi^2 value, then the shared variable chi2_min is replaced
            //We then also replace the top_thread variable to know which thread currently has the best structure
            if(chi2<chi2_min){
                chi2_min = chi2;
                top_thread=omp_get_thread_num();

                //We save the individual curve chi^2 values for this structure
                //These will be used to determine cutoff values later on.
                for(kk=0; kk<N_curves; kk++){
                    best_struct_curve_chi2[kk] = X2[kk][d_indices[kk]][std_indices[kk]];
                }

                //This structure overwrites xyz_best
                copy_structure(N_atoms,xyz_priv,xyz_best);
                printf("New Chi2 minimum at %lf\n",chi2_min);

                //Lastly, we save the distance and std indices for this structure, to write the best-fit curve
                for(kk=0; kk<N_curves; kk++){
                    d_indices_range[kk][0]=d_indices[kk];
                    std_indices_range[kk][0]=std_indices[kk];
                }

                //It is then saved as a mol2 file
                //each processor has its own "best" structure and these are compared at the end to determine
                //which one is the actual best-fit structure
                remove(min_chi2_output_filename);
                sprintf(min_chi2_output_filename, "%s_%d_best.mol2", base_iteration_filename,omp_get_thread_num());
                write_mol2(min_chi2_output_filename, N_atoms, N_bonds, atom_id, element, xyz_priv, atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);
            }

            //This this structure has the new lowest Chi^2 value for a given curve, the existing value is replaced.
        }}
    }

    if (chi2_min == 5000000){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: No best-fit structure found with the given inputs\n");
        fclose(error_file);
        exit(1);
    }

    //We delete the individual thread best structures, leaving only the true best-fir structure mol2 file
    char best_filename[128];
    sprintf(best_filename,"%s_best.mol2", base_iteration_filename);
    #pragma omp parallel
    {
        char thread_filename[128];
        sprintf(thread_filename,"%s_%d_best.mol2", base_iteration_filename,omp_get_thread_num());

        if(omp_get_thread_num()==top_thread)
            rename(thread_filename,best_filename);

        remove(thread_filename);
    }

    //surface atoms are added to the mol2 file for plotting purposes.
    add_surface(best_filename);

    //store the distance information for the best-fit structure
    for(i=0;i<N_curves;i++){
        distance_BSL[i][0] = get_average_distance(REDOR_det_index,REDOR_rec_index, xyz_best, i,curve_type[i]);
        stdev_BSL[i][0] = get_STDEV(REDOR_det_index,REDOR_rec_index, xyz_best, i,curve_type[i]);
    }

    //the max_Chi2 function returns the highest allowable Chi^2 value given a minimum value and an error range
    //See numerical recipes in C for a description of the calculation.
    for(i=0; i<N_curves; i++){
        curve_chi2_max[i] = max_Chi2(curve_chi2_min[i], threshold_accuracy/100.);
    }

    //If a cutoff value excludes the best fit structure, the cutoff is replaced with that from the best fit structure
    for(i=0; i<N_curves; i++){
        if(curve_chi2_max[i]<best_struct_curve_chi2[i]){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nNOTE: Best-fit structure is not within experimental error for curve %d. Cutoff values were altered for this curve\n", i);
            fclose(error_file);
            curve_chi2_max[i] = best_struct_curve_chi2[i];
        }
    }

    //Next we find ALL the structures that fit within experimental error for all curves.
    //to prevent the writing of too many files a maximum number of structures is given
    //which will truncate the loop if too many structures are found

    int acceptable_structures = 0;
    int other_structures = 0;
    int max_accept_failsafe = 0;
    printf("\nFinding all structures that agree with experiment\n");

    #pragma omp parallel for
    for(int it=0;it<iterations;it++){
        if((acceptable_structures+other_structures)>max_acceptable_struct){
            it = iterations;
            max_accept_failsafe = 1;
        }
        int ii, jj, kk, bond_position[N_rotatable_bonds];
        int d_indices[N_curves], std_indices[N_curves];
        int check_chi2_threshold;
        double xyz_priv[N_atoms][3], nominator, angle, deviation, R[4][3];
        char min_chi2_output_filename[128];

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
            angle=2*Pi/bond[jj].N_steps*bond_position[jj];
            generate_bond_rot_matrix(R, xyz_priv[bond[jj].atom1], xyz_priv[bond[jj].atom2], angle);
            //rotating all of atoms involved around the bond
            for(ii=0;ii<bond[jj].N_aff_atoms; ii++){
                rotate_around_bond2(xyz_priv[bond[jj].affected_atom[ii]], R);
            }
        }

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

        if(!collisions(xyz_priv, N_atoms, collision_distance)){
        if(check_constraints == 0){
            check_chi2_threshold = 0;
            for(kk=0; kk<N_curves; kk++){
                d_indices[kk] = get_distance_index(REDOR_det_index,REDOR_rec_index, xyz_priv, kk, curve_type[kk]);
                std_indices[kk] = get_STDEV_index(REDOR_det_index,REDOR_rec_index, xyz_priv, kk, curve_type[kk]);
                curve_chi2[kk] = X2[kk][d_indices[kk]][std_indices[kk]];
                    if(curve_chi2[kk]>curve_chi2_max[kk]){
                        check_chi2_threshold++;
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
                    sprintf(min_chi2_output_filename, "%s_struct%d.mol2", base_iteration_filename, acceptable_structures);
                    printf("(%d) Acceptable structure found\n", acceptable_structures);
                    write_mol2(min_chi2_output_filename, N_atoms, N_bonds, atom_id, element, xyz_priv, atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);

                    //The minimum and max distances are updated, if necessary, to print out the fitted curves ranges at the end.
                    for(kk=0; kk<N_curves; kk++){
                        std_indices_range[kk][1] = (d_indices_range[kk][1] <= d_indices[kk])*std_indices_range[kk][1] + (d_indices_range[kk][1] > d_indices[kk])*std_indices[kk];
                        std_indices_range[kk][2] = (d_indices_range[kk][2] >= d_indices[kk])*std_indices_range[kk][2] + (d_indices_range[kk][2] < d_indices[kk])*std_indices[kk];
                        d_indices_range[kk][1] = (d_indices_range[kk][1] <= d_indices[kk])*d_indices_range[kk][1] + (d_indices_range[kk][1] > d_indices[kk])*d_indices[kk];
                        d_indices_range[kk][2] = (d_indices_range[kk][2] >= d_indices[kk])*d_indices_range[kk][2] + (d_indices_range[kk][2] < d_indices[kk])*d_indices[kk];
                    }

                }
                else{
                    other_structures++;
                    sprintf(min_chi2_output_filename, "other_%s_struct%d.mol2", base_iteration_filename, other_structures);
                    printf("(%d) Other acceptable structure found: rmsd with best = %lf\n", other_structures, deviation);
                    write_mol2(min_chi2_output_filename, N_atoms, N_bonds, atom_id, element, xyz_priv, atom_type, bond_id, ori_atom_id, tar_atom_id, bond_type);
                }
            }
        }}
    }



    if(acceptable_structures <= 1){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Only one acceptable structure found with the given constraints\n");
        fclose(error_file);
        exit(1);
    }

    printf("\nFound a total of %d acceptable structures,\n%d of which were within the requested RMSD from the best one\n", acceptable_structures+other_structures, acceptable_structures);
    fclose(mol2_file);

    //Next we will write the structure overlay files
    char overlay_filename[128];
    sprintf(overlay_filename, "%s_overlay.mol2", base_iteration_filename);

    //if there are too many structures, erase the structures and exit the program
    if(max_accept_failsafe == 1){
        for(i=0;i<acceptable_structures;i++){
            sprintf(overlay_filename, "%s_struct%d.mol2", base_iteration_filename,i);
            remove(overlay_filename);
        }
        for(i=0;i<other_structures;i++){
            sprintf(overlay_filename, "other_%s_struct%d.mol2", base_iteration_filename,i);
            remove(overlay_filename);
        }
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Too many acceptable structures, exceeded limit of %d\n", max_acceptable_struct);
        fclose(error_file);
        return 0;
    }

    //Creating overlay and probability ellipsoids..
    printf("\nOverlaying Structures\n");

    //for primary structures
    compile_mol2_files(base_iteration_filename, acceptable_structures);
    create_cif(overlay_filename, N_atoms+5);

    //for all structures
    compile_all_mol2_files(base_iteration_filename, acceptable_structures, other_structures);

    //write out the fitted REDOR curve and ranges.
    write_fits(d_indices_range, std_indices_range, base_iteration_filename, N_curves, support, REDOR_det_index, REDOR_rec_index, element,curve_filename, scaling_factor, Nspins, curve_type);
}//end int main