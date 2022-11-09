#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
using namespace std;

void generate_rotation_matrix(double (*R)[3],double phi, double theta, double psi){
    //This function generates a rotation matrix, stored at the pointer *R
    //see JACS 2006, 128, 11860-11871 (SI).
    double ct=cos(theta);
    double st=sin(theta);
    double cp=cos(phi);
    double sp=sin(phi);
    double cg=cos(psi);
    double sg=sin(psi);

    R[0][0]=-sg*ct*sp+cg*cp;
    R[0][1]= sg*ct*cp+cg*sp;
    R[0][2]= sg*st;

    R[1][0]=-cg*ct*sp-sg*cp;
    R[1][1]= cg*ct*cp-sg*sp;
    R[1][2]= cg*st;

    R[2][0]= st*sp;
    R[2][1]=-st*cp;
    R[2][2]= ct;
}

void create_cif( char *compiled_mol2_filename, int N_atoms){
   //This function will read a mol2 file containing multiple structures in an overlay
    //and generate a cif file containing the average structure with error ellipsoids
    char error_filename[128];
    sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    double Pi = 3.1415926535897932384626433;
    char buffer[128], crap[3], element[N_atoms][3];
    int N_structures,N_total_atoms,i,j,k,line_Atoms,junk, atom_id[N_atoms];
    int compiled_mol2_filename_len=strlen(compiled_mol2_filename);

    char cif_filename[compiled_mol2_filename_len+1];
    sprintf(cif_filename,"%.*s.cif",compiled_mol2_filename_len-5,compiled_mol2_filename);

    remove(cif_filename);

    FILE *fp, *cif;

    fp=fopen(compiled_mol2_filename,"r");
    if(fp==NULL){
        error_file=fopen(error_filename,"a");
        fprintf(error_file, "\nERROR: Compiled mol2 file '%s' not found when trying to create .cif file\n", compiled_mol2_filename);
        fclose(error_file);
        exit(1);
    }

    i=1;

    while(fgets(buffer, sizeof(buffer), fp) != NULL){
            if(strcmp(buffer, "@<TRIPOS>MOLECULE\n")==0){
                fgets(buffer, sizeof(buffer), fp);
                fgets(buffer, sizeof(buffer), fp);
                sscanf(buffer,"%d",&N_total_atoms);
                N_structures=N_total_atoms/N_atoms;
                i++;
                i++;
            }
            if(strcmp(buffer, "@<TRIPOS>ATOM\n")==0){
                line_Atoms = i;
                }
            i++;
    }
    fclose(fp);

    vector< vector< vector<double> > > coordinates;
    coordinates.resize(N_atoms, vector< vector<double> >(N_structures,vector<double>(3,0.)));
    vector< vector<double> >  ave_coordinates;
    ave_coordinates.resize(N_atoms, vector<double>(3,0.));

    fp=fopen(compiled_mol2_filename,"r");
    i=1;
    while(fgets(buffer, sizeof(buffer), fp) != NULL){
        if(i==line_Atoms+1){
            i++;
            for(j=0; j<N_atoms; j++){
                sscanf(buffer,"%d %s %lf %lf %lf", &atom_id[j], element[j],&coordinates[j][0][0],&coordinates[j][0][1],&coordinates[j][0][2]);
                fgets(buffer, sizeof(buffer), fp);
            }
            for(k=1;k<N_structures;k++){
                for(j=0; j<N_atoms; j++){
                    sscanf(buffer,"%d %s %lf %lf %lf", &junk, crap,&coordinates[j][k][0],&coordinates[j][k][1],&coordinates[j][k][2]);
                    fgets(buffer, sizeof(buffer), fp);
            }}}
        else
            i++;
    }

    fclose(fp);

    for(i=0;i<N_atoms;i++){
        for(j=0;j<N_structures;j++){
            for(k=0;k<3;k++){
                ave_coordinates[i][k]=ave_coordinates[i][k]+coordinates[i][j][k]/N_structures;
            }
        }
    }

    cif=fopen(cif_filename,"w");
    fprintf(cif,"data_%s\n",cif_filename);
    fprintf(cif,"_symmetry_space_group_name_H-M    'P1'\n");
    fprintf(cif,"_symmetry_Int_Tables_number       1\n");
    fprintf(cif,"_symmetry_cell_setting            triclinic\n");
    fprintf(cif,"loop_\n");
    fprintf(cif,"_symmetry_equiv_pos_as_xyz\n");
    fprintf(cif,"  x,y,z\n");
    fprintf(cif,"_cell_length_a                    100.0\n");
    fprintf(cif,"_cell_length_b                    100.0\n");
    fprintf(cif,"_cell_length_c                    100.0\n");
    fprintf(cif,"_cell_angle_alpha                 90.0000\n");
    fprintf(cif,"_cell_angle_beta                  90.0000\n");
    fprintf(cif,"_cell_angle_gamma                 90.0000\n");
    fprintf(cif,"loop_\n");
    fprintf(cif,"_atom_site_label\n");
    fprintf(cif,"_atom_site_type_symbol\n");
    fprintf(cif,"_atom_site_fract_x\n");
    fprintf(cif,"_atom_site_fract_y\n");
    fprintf(cif,"_atom_site_fract_z\n");
    fprintf(cif,"_atom_site_U_iso_or_equiv\n");
    fprintf(cif,"_atom_site_adp_type\n");
    fprintf(cif,"_atom_site_occupancy\n");

    for(i=0;i<N_atoms;i++){
        fprintf(cif,"%s%d    %s    %lf    %lf    %lf    0.100000    Uani    1.00\n",element[i],atom_id[i],element[i],ave_coordinates[i][0]/100.+.5,ave_coordinates[i][1]/100.+.5,ave_coordinates[i][2]/100.);
    }
    fprintf(cif,"\nloop_\n");
    fprintf(cif,"_atom_site_aniso_label\n");
    fprintf(cif,"_atom_site_aniso_U_11\n");
    fprintf(cif,"_atom_site_aniso_U_22\n");
    fprintf(cif,"_atom_site_aniso_U_33\n");
    fprintf(cif,"_atom_site_aniso_U_23\n");
    fprintf(cif,"_atom_site_aniso_U_13\n");
    fprintf(cif,"_atom_site_aniso_U_12\n");

    fclose(cif);
    cif=fopen(cif_filename,"a");

    //shifting all atoms to the origin to prepare the calculation of U tensors
    for(i=0;i<N_atoms;i++){
        for(j=0;j<N_structures;j++){
            for(k=0;k<3;k++){
                coordinates[i][j][k]=coordinates[i][j][k]-ave_coordinates[i][k];
            }
        }
    }

    double U[3][3], E[3][3];
    double rms_max[N_atoms][3], phi_max[N_atoms], theta_max[N_atoms], psi_max[N_atoms];

    //finding the theta and psi Euler angles along with the Z RMS values
    int counter=0;
    #pragma omp parallel for
    for(i=0;i<N_atoms;i++){
        double R[3][3], rms;
        int jj;
        rms_max[i][0]=rms_max[i][1]=rms_max[i][2]=0.;
        phi_max[i]=theta_max[i]=psi_max[i]=0.;
        for(double theta=0.;theta<Pi;theta=theta+Pi/180.){
            for(double psi=0.;psi<2.*Pi;psi=psi+Pi/180.){
                generate_rotation_matrix(R,0.,theta,psi);
                rms=0.;
                for(jj=0;jj<N_structures;jj++){
                    rms=rms+pow((R[0][2]*coordinates[i][jj][0]+R[1][2]*coordinates[i][jj][1]+R[2][2]*coordinates[i][jj][2]),2.)/N_structures;
                }
                rms=sqrt(rms);
                if(rms>rms_max[i][2]){
                    rms_max[i][2]= rms;
                    theta_max[i]= theta;
                    psi_max[i]= psi;
                }
            }
        }
        counter++;
        printf("\nCalculating probability ellipsoids for atom %d of %d",counter,N_atoms);
    }

    //finding the value of phi and the X and Y RMS values
    #pragma omp parallel for
    for(i=0;i<N_atoms;i++){
        double phi, theta, psi, R[3][3], rms;
        int jj;
        for(phi=0.;phi<2.*Pi;phi=phi+Pi/180.){
            generate_rotation_matrix(R,phi,theta_max[i],psi_max[i]);
            rms=0.;
            for(jj=0;jj<N_structures;jj++){
                rms=rms+pow((R[0][1]*coordinates[i][jj][0]+R[1][1]*coordinates[i][jj][1]+R[2][1]*coordinates[i][jj][2]),2.)/N_structures;
            }

            if(rms>rms_max[i][1]){
                rms_max[i][1]= rms;
                phi_max[i]= phi;
            }
        }
        generate_rotation_matrix(R,phi_max[i],theta_max[i],psi_max[i]);
        for(jj=0;jj<N_structures;jj++){
            rms_max[i][0]=rms_max[i][0]+pow((R[0][0]*coordinates[i][jj][0]+R[1][0]*coordinates[i][jj][1]+R[2][0]*coordinates[i][jj][2]),2.)/N_structures;
        }
    }

    //calculating the U tensor and printing it to the cif file
    for(i=0;i<N_atoms;i++){
        double R[3][3];
        rms_max[i][0]=sqrt(rms_max[i][0]);
        generate_rotation_matrix(R,phi_max[i],theta_max[i],psi_max[i]);

        //giving the RMS values a minimum of 0.1 A
        rms_max[i][0]=rms_max[i][0]+0.1;
        rms_max[i][1]=rms_max[i][1]+0.1;
        rms_max[i][2]=rms_max[i][2]+0.1;


        U[0][0]=R[0][0]*rms_max[i][0]*R[0][0] + R[0][1]*rms_max[i][1]*R[0][1] + R[0][2]*rms_max[i][2]*R[0][2];
        U[0][1]=R[0][0]*rms_max[i][0]*R[1][0] + R[0][1]*rms_max[i][1]*R[1][1] + R[0][2]*rms_max[i][2]*R[1][2];
        U[0][2]=R[0][0]*rms_max[i][0]*R[2][0] + R[0][1]*rms_max[i][1]*R[2][1] + R[0][2]*rms_max[i][2]*R[2][2];

        U[1][0]=R[1][0]*rms_max[i][0]*R[0][0] + R[1][1]*rms_max[i][1]*R[0][1] + R[1][2]*rms_max[i][2]*R[0][2];
        U[1][1]=R[1][0]*rms_max[i][0]*R[1][0] + R[1][1]*rms_max[i][1]*R[1][1] + R[1][2]*rms_max[i][2]*R[1][2];
        U[1][2]=R[1][0]*rms_max[i][0]*R[2][0] + R[1][1]*rms_max[i][1]*R[2][1] + R[1][2]*rms_max[i][2]*R[2][2];

        U[2][0]=R[2][0]*rms_max[i][0]*R[0][0] + R[2][1]*rms_max[i][1]*R[0][1] + R[2][2]*rms_max[i][2]*R[0][2];
        U[2][1]=R[2][0]*rms_max[i][0]*R[1][0] + R[2][1]*rms_max[i][1]*R[1][1] + R[2][2]*rms_max[i][2]*R[1][2];
        U[2][2]=R[2][0]*rms_max[i][0]*R[2][0] + R[2][1]*rms_max[i][1]*R[2][1] + R[2][2]*rms_max[i][2]*R[2][2];

        E[0][0]=U[0][0]*U[0][0]+U[1][0]*U[1][0]+U[2][0]*U[2][0];
        E[1][1]=U[0][1]*U[0][1]+U[1][1]*U[1][1]+U[2][1]*U[2][1];
        E[2][2]=U[0][2]*U[0][2]+U[1][2]*U[1][2]+U[2][2]*U[2][2];

        E[0][1]=U[0][0]*U[0][1]+U[1][0]*U[1][1]+U[2][0]*U[2][1];
        E[0][2]=U[0][0]*U[0][2]+U[1][0]*U[1][2]+U[2][0]*U[2][2];
        E[1][2]=U[0][1]*U[0][2]+U[1][1]*U[1][2]+U[2][1]*U[2][2];


        fprintf(cif,"%s%d    %lf    %lf    %lf    %lf    %lf    %lf\n",element[i],atom_id[i],E[0][0],E[1][1],E[2][2],E[1][2],E[0][2],E[0][1]);

    }
    fclose(cif);
}
