#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <time.h>
#include <algorithm>
static double Pi = 3.1415926535897932384626433;
using namespace std;

struct Constraint{
    int atom1;
    int atom2;
    int atom3;
    int atom4;
    double minimum;
    double maximum;
    int type; //0 for distance, 1 for angle, 2 for dihedral, 3 for surface_distance
};

double distance_sq(double *atom1, double *atom2){
    //function that returns the squared distance.
    //In enables the calculation of collisions without having to call a sqrt
    double x=(atom1[0]-atom2[0]);
    double y=(atom1[1]-atom2[1]);
    double z=(atom1[2]-atom2[2]);
    return x*x+y*y+z*z;
}

double distance_calc(double *atom1, double *atom2){
    //returns the distance between two atoms
    double x=(atom1[0]-atom2[0]);
    double y=(atom1[1]-atom2[1]);
    double z=(atom1[2]-atom2[2]);
    return sqrt(x*x+y*y+z*z);
}

double angle_calc(double *atom1, double *atom2, double *atom3){
    //returns the bond angle between 3 atoms;
    double vector1[3], vector2[3], length, dot=0.;
    int i;

    length=distance_calc(atom1,atom2);
    for(i=0;i<3;i++){
        vector1[i]=atom1[i]-atom2[i];
        vector1[i]=vector1[i]/length;
    }

    length=distance_calc(atom2,atom3);
    for(i=0;i<3;i++){
        vector2[i]=atom3[i]-atom2[i];
        vector2[i]=vector2[i]/length;
    }

    for(i=0;i<3;i++){
        dot=dot+vector1[i]*vector2[i];
    }
    return acos(dot)*180./Pi;
}

double dihedral_calc(double *atom1, double *atom2, double *atom3, double *atom4){
    //returns the dihedral angle between 4 atoms.
    double vector1[3], vector2[3], vector3[3], normal1[3], normal2[3], origin[3], dihedral, triple=0;
    int i;

    for(i=0;i<3;i++){
        vector1[i]=atom1[i]-atom2[i];
        vector2[i]=atom2[i]-atom3[i];
        vector3[i]=atom3[i]-atom4[i];
        origin[i]=0.;
    }

    normal1[0]=vector1[1]*vector2[2]-vector1[2]*vector2[1];
    normal1[1]=vector1[2]*vector2[0]-vector1[0]*vector2[2];
    normal1[2]=vector1[0]*vector2[1]-vector1[1]*vector2[0];

    normal2[0]=vector2[1]*vector3[2]-vector2[2]*vector3[1];
    normal2[1]=vector2[2]*vector3[0]-vector2[0]*vector3[2];
    normal2[2]=vector2[0]*vector3[1]-vector2[1]*vector3[0];

    for(i=0;i<3;i++){
        triple=triple+vector1[i]*normal2[i];
    }

    dihedral=angle_calc(normal1,origin,normal2);

    return ((triple<0.) - (triple>=0.))*dihedral;
}

int not_bonded(vector< vector<int> > neighbors, int atom1, int atom2){
    int i, bonded=1;
    for(i=0; i<neighbors[atom1].size(); i++){
        bonded=bonded*(neighbors[atom1][i] != atom2);
    }
    return bonded;
}

int collisions(double (*xyz)[3], vector< vector<int> > neighbors, int N_atoms, double surface_collision_distance, double interatomic_collision_distance){
    //This function returns true if there is either an interatomic (< 1A) or surface-atom collision.
    //A surface-atom distance of exactly zero is ignored as it is assumed to be a surface atom.
    int i,j;
    double d,x,y,z;

    for(i=0;i<N_atoms;i++){
        if(xyz[i][2]>surface_collision_distance)
            continue;
        if(xyz[i][2]!=0.0)
            return 1;
        }

    for(i=0;i<N_atoms;i++){
        for(j=i+1;j<N_atoms; j++){
            x=xyz[i][0]-xyz[j][0];
            if(x>interatomic_collision_distance)
                continue;

            y=xyz[i][1]-xyz[j][1];
            if(y>interatomic_collision_distance)
                continue;

            z=xyz[i][2]-xyz[j][2];
            if(z>interatomic_collision_distance)
                continue;

            d=x*x+y*y+z*z;
            if(d> interatomic_collision_distance*interatomic_collision_distance)
                continue;

            //check if the two atoms are bonded
            if(find(neighbors[i].begin(),neighbors[i].end(),j)==neighbors[i].end())
                return 1;
    }
}
    return 0;
}

int collisions_old(double (*xyz)[3], vector< vector<int> > neighbors, int N_atoms, double surface_collision_distance, double interatomic_collision_distance){
    //This function returns true if there is either an interatomic (< 1A) or surface-atom collision.
    //A surface-atom distance of exactly zero is ignored as it is assumed to be a surface atom.
    int i,j;

    for(i=0;i<N_atoms;i++){
        if(xyz[i][2]<surface_collision_distance){
            if(xyz[i][2]!=0.0)
            return 1;
        }
        for(j=i+1;j<N_atoms; j++){
            if(distance_calc(xyz[i],xyz[j])<1.0)
                return 1;
    }}
    return 0;
}
