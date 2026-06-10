using namespace std;
#include <vector>

struct fragments{
    int N_fragments;
    vector< vector<int> > frag_indices;  //indices of each atoms in a given fragment
    vector< vector<double> > mol_xyz; //xyz coordinates of the center of each fragment
    vector< vector<double> > mol_frac; //fractional coordinates of the center of each fragment
};

void xyz_to_frac(int N_atoms, vector< vector<double> > &xyz, vector< vector<double> > &frac, vector< vector<double> > &cell){
    int i;
    for(i=0;i<N_atoms;i++){
        frac[i][1] = xyz[i][1] / cell[1][1];
        frac[i][0] = (xyz[i][0] - cell[1][0]*frac[i][1]) / cell[0][0];
        frac[i][2] = xyz[i][2] / cell[2][2];
    }
}
void frac_to_xyz(int N_atoms, vector< vector<double> > &xyz, vector< vector<double> > &frac, vector< vector<double> > &cell){
    int i;
    for(i=0;i<N_atoms;i++){
        xyz[i][0] = cell[0][0]*frac[i][0];
        xyz[i][1] = cell[1][0]*frac[i][0] + cell[1][1]*frac[i][1];
        xyz[i][2] = cell[2][0]*frac[i][0] + cell[2][1]*frac[i][1]+ cell[2][2]*frac[i][2];
    }
}

void calc_cell_matrix(vector< vector<double> > &cell, vector<double> &unit_cell){
    //Function that takes the unit cell parameters and calculated the cell cartesian matrix
    cell[0][1]=cell[0][2]=cell[1][2]=0.;
    cell[0][0]=unit_cell[0]; //a

    cell[1][0]=unit_cell[1]*cos(unit_cell[5]); //b*cos(gamma)
    cell[1][1]=unit_cell[1]*sin(unit_cell[5]); //b*sin(gamma)

    cell[2][0]= unit_cell[2]*cos(unit_cell[4]); //c*cos(beta)
    cell[2][1]= unit_cell[2]*(cos(unit_cell[3])-cos(unit_cell[4])*cos(unit_cell[5]))/sin(unit_cell[5]); //c*(cos(alpha)-cos(beta)*cos(gamma))/sin(gamma)
    cell[2][2]= sqrt(unit_cell[2]*unit_cell[2]-cell[2][0]*cell[2][0]-cell[2][1]*cell[2][1]);
}

void calc_cell_dimensions(vector< vector<double> > &cell, vector<double> &unit_cell){
    //Function that takes the cell cartesian matrix and calculates the unit cell parameters
    unit_cell[0]=cell[0][0];//a
    unit_cell[1]=sqrt(cell[1][0]*cell[1][0]+cell[1][1]*cell[1][1]);//b
    unit_cell[2]=sqrt(cell[2][0]*cell[2][0]+cell[2][1]*cell[2][1]+cell[2][2]*cell[2][2]);//c

    unit_cell[5]=acos(cell[1][0]/unit_cell[1]);//gamma
    unit_cell[4]=acos(cell[2][0]/unit_cell[2]);//beta
    unit_cell[3]=acos(cell[2][1]/unit_cell[2]*sin(unit_cell[5])+cos(unit_cell[4])*cos(unit_cell[5]));
}

double calc_cell_volume(vector<double> &unit_cell){
    double ca=cos(unit_cell[3]);
    double cb=cos(unit_cell[4]);
    double cg=cos(unit_cell[5]);

    return unit_cell[0]*unit_cell[1]*unit_cell[2]*sqrt(1.+2.*ca*cb*cg-ca*ca-cb*cb-cg*cg);
}

void get_fragments(int N_atoms, vector< vector<int> > &neighbors, fragments *fragment, vector< vector<double> > &xyz, vector< vector<double> > &cell){
    //Function that counts and tags the fragments in the supplied mol2 file
    //Data are stored as a struct fragments;
    int i,j,k,l,track[N_atoms], lowest, check, start;
    fragment->N_fragments=0;

    for(i=0;i<N_atoms;i++){
        track[i]=i+1;
    }

    do{
        for(i=0;i<N_atoms;i++){
            lowest=0;
            if(track[i]!=0){
                lowest=1;
                fragment->frag_indices.push_back(vector<int>());
                fragment->frag_indices[fragment->N_fragments].push_back(i);
                fragment->N_fragments++;
                track[i]=0;
                break;
            }
        }
        if(lowest==0)
            break;

        start=0;
        for(i=start;i<fragment->frag_indices[fragment->N_fragments-1].size();i++){

            for(j=1;j<neighbors[fragment->frag_indices[fragment->N_fragments-1][i]].size();j++){
                check=0;
                for(k=0;k<fragment->frag_indices[fragment->N_fragments-1].size();k++){
                    if(neighbors[fragment->frag_indices[fragment->N_fragments-1][i]][j]==fragment->frag_indices[fragment->N_fragments-1][k]){
                        check=1;
                        break;
                    }
                }
                if(check==0){
                    fragment->frag_indices[fragment->N_fragments-1].push_back(neighbors[fragment->frag_indices[fragment->N_fragments-1][i]][j]);
                    track[neighbors[fragment->frag_indices[fragment->N_fragments-1][i]][j]]=0;
                    track[neighbors[fragment->frag_indices[fragment->N_fragments-1][i]][0]]=0;
                    for(l=0;l<N_atoms;l++){
                    }
                    start++;
                }
            }//neighbors of the atom in the fragment
        }//atoms in fragment
    }while(true);

    //Now that all fragments have been found and labeled,
    //we will determine the fractional coordinates of the center of each fragment
    //First we find the mean xyz coordinates of each fragment
    fragment->mol_xyz.resize(fragment->N_fragments,vector<double>(3,0.));
    fragment->mol_frac.resize(fragment->N_fragments,vector<double>(3,0.));

    for(i=0;i<fragment->N_fragments;i++){
        for(j=0;j<fragment->frag_indices[i].size();j++){
            fragment->mol_xyz[i][0] += xyz[fragment->frag_indices[i][j]][0] / fragment->frag_indices[i].size();
            fragment->mol_xyz[i][1] += xyz[fragment->frag_indices[i][j]][1] / fragment->frag_indices[i].size();
            fragment->mol_xyz[i][2] += xyz[fragment->frag_indices[i][j]][2] / fragment->frag_indices[i].size();
        }
    }

    //Using the cell information, the mean xyz coordinates are converted to fractiona ones
    xyz_to_frac(fragment->N_fragments,fragment->mol_xyz,fragment->mol_frac,cell);
}
