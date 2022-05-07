#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
using namespace std;

int translate_atom_Z(double *xyz, double distance){
    //Translates xyz coordinates of one atom up by distance along the Z direction
    xyz[2]=xyz[2]+distance;
    return 0;
}

void translate_molecule_Z(int N_atoms, double (*xyz)[3], double distance){
   //Translates xyz coordinates of entire molecule up by distance along the Z direction
    int i;

    for(i=0; i<N_atoms; i++){
    translate_atom_Z(xyz[i],distance);
    }
}

int translate_atom_X(double *xyz, double distance){
    //Translates xyz coordinates of one atom up by distance along the X direction
    xyz[0]=xyz[0]+distance;
    return 0;
}

int translate_molecule_X(int N_atoms, double (*xyz)[3], double distance){
    //Translates xyz coordinates of entire molecule up by distance along the X direction
    int i;

    for(i=0; i<N_atoms; i++){
    translate_atom_X(xyz[i],distance);
    }
    return 0;
}

int translate_atom_Y(double *xyz, double distance){
    //Translates xyz coordinates of one atom up by distance along the Y direction

    xyz[1]=xyz[1]+distance;
    return 0;
}

int translate_molecule_Y(int N_atoms, double (*xyz)[3], double distance){
    //Translates xyz coordinates of entire molecule up by distance along the Y direction
    int i;

    for(i=0; i<N_atoms; i++){
    translate_atom_Y(xyz[i],distance);
    }
    return 0;
}

void rotate_around_bond(double *atom, double *atom_rot, double *bond1, double *bond2, double angle){
    //This function rotates atom around the bond formed by bond1 and bond2 by an angle of angle.
    //The rotated atom is stored in atom_rot
    double bond_axis[3], bond_center[3], xyz_rot[3], R_mat[3][3];
    double ct=cos(angle);
    double st=sin(angle);

    //Normalized vector oriented with the bond
    bond_axis[0]=bond1[0]-bond2[0];
    bond_axis[1]=bond1[1]-bond2[1];
    bond_axis[2]=bond1[2]-bond2[2];
    double len=sqrt(bond_axis[0]*bond_axis[0]+bond_axis[1]*bond_axis[1]+bond_axis[2]*bond_axis[2]);
    bond_axis[0]=bond_axis[0]/len;
    bond_axis[1]=bond_axis[1]/len;
    bond_axis[2]=bond_axis[2]/len;

    //center of the bond
    bond_center[0]=(bond1[0]+bond2[0])/2.;
    bond_center[1]=(bond1[1]+bond2[1])/2.;
    bond_center[2]=(bond1[2]+bond2[2])/2.;

    //rotation matrix
    R_mat[0][0]=ct+bond_axis[0]*bond_axis[0]*(1.-ct);
    R_mat[0][1]=bond_axis[0]*bond_axis[1]*(1.-ct)-bond_axis[2]*st;
    R_mat[0][2]=bond_axis[0]*bond_axis[2]*(1.-ct)+bond_axis[1]*st;

    R_mat[1][0]=bond_axis[0]*bond_axis[1]*(1.-ct)+bond_axis[2]*st;
    R_mat[1][1]=ct+bond_axis[1]*bond_axis[1]*(1.-ct);
    R_mat[1][2]=bond_axis[1]*bond_axis[2]*(1.-ct)-bond_axis[0]*st;

    R_mat[2][0]=bond_axis[0]*bond_axis[2]*(1.-ct)-bond_axis[1]*st;
    R_mat[2][1]=bond_axis[1]*bond_axis[2]*(1.-ct)+bond_axis[0]*st;
    R_mat[2][2]=ct+bond_axis[2]*bond_axis[2]*(1.-ct);

    //calculating new coordinates for atom 2
    atom_rot[0]=atom[0]-bond_center[0];
    atom_rot[1]=atom[1]-bond_center[1];
    atom_rot[2]=atom[2]-bond_center[2];

    xyz_rot[0]=atom_rot[0]*R_mat[0][0]  +atom_rot[1]*R_mat[0][1]  +atom_rot[2]*R_mat[0][2];
    xyz_rot[1]=atom_rot[0]*R_mat[1][0]  +atom_rot[1]*R_mat[1][1]  +atom_rot[2]*R_mat[1][2];
    xyz_rot[2]=atom_rot[0]*R_mat[2][0]  +atom_rot[1]*R_mat[2][1]  +atom_rot[2]*R_mat[2][2];

    atom_rot[0]=xyz_rot[0]+bond_center[0];
    atom_rot[1]=xyz_rot[1]+bond_center[1];
    atom_rot[2]=xyz_rot[2]+bond_center[2];
}

void rotate_around_bond2(double *atom, double *bond1, double *bond2, double angle){
    //same as rotate_around_bond() but overwrites the existing coordinates
    double dummy[3];
    rotate_around_bond(atom,dummy,bond1,bond2,angle);
    atom[0]=dummy[0];
    atom[1]=dummy[1];
    atom[2]=dummy[2];
}

void rotate_molecule_around_Z(int N_atoms, double (*xyz)[3], double angle){
    //Rotates an entire molecule around the origin in the Z direction
    int i;
    double origin[3]={0.,0.,0.}, z_vector[3]={0.,0.,1.};

    for(i=0; i<N_atoms; i++){
        rotate_around_bond2(xyz[i], origin, z_vector, angle);
    }
}

void rotate_molecule_around_X(int N_atoms, double (*xyz)[3], double angle){
    //Rotates an entire molecule around its center in the X direction
    int i,j;
    double X_vec[3],origin[3] = {0.,0.,0.};

    for(i=0;i<N_atoms;i++){
        for(j=0;j<3;j++){
            origin[j] = origin[j] + xyz[i][j]/N_atoms;
    }}

    for(j=0;j<3;j++){
        X_vec[j] = origin[j];
    }
    X_vec[0]=origin[0]+1.;


    for(i=0;i<N_atoms;i++){
        rotate_around_bond2(xyz[i],origin,X_vec,angle);
    }
}

void rotate_molecule_around_Y(int N_atoms, double (*xyz)[3], double angle){
    //Rotates an entire molecule around its center in the Y direction
    int i,j;
    double Y_vec[3],origin[3] = {0.,0.,0.};

    for(i=0;i<N_atoms;i++){
        for(j=0;j<3;j++){
            origin[j] = origin[j] + xyz[i][j]/N_atoms;
    }}

    for(j=0;j<3;j++){
        Y_vec[j] = origin[j];
    }
    Y_vec[1]=origin[1]+1.;


    for(i=0;i<N_atoms;i++){
        rotate_around_bond2(xyz[i],origin,Y_vec,angle);
    }
}

void copy_structure(int N_atoms, double (*xyz)[3], double (*xyz_copy)[3]){
    //Makes a copy of xyz in xyz_copy
    int i;
    for(i=0;i<N_atoms;i++){
        xyz_copy[i][0]=xyz[i][0];
        xyz_copy[i][1]=xyz[i][1];
        xyz_copy[i][2]=xyz[i][2];
    }
}



