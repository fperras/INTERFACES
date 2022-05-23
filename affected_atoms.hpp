#include "mol2_files.hpp"

struct Bond{
    int atom1;
    int atom2;
    int atom0;
    int N_aff_atoms;
    int N_steps;
    vector<int> affected_atom;
    int mod;
    int type; //0=revolve 1=stretch 2=bend
    double dmin;
    double dmax;
};

void get_affected_atoms(int N_rotatable_bonds, struct Bond *bond, vector< vector<int> > neighbors){
    //Same function of get_aff_atoms() but uses neighbors and should be faster.

    int found,start,check=0,i,j,k,l;

    for(i=0;i<N_rotatable_bonds;i++){
        bond[i].N_aff_atoms=1;
        bond[i].affected_atom.resize(bond[i].N_aff_atoms);
        bond[i].affected_atom[0]=bond[i].atom2;
    }

    for(k=0;k<N_rotatable_bonds;k++){
        start=0;
        do{
            found=1;
            for(i=start;i<bond[k].N_aff_atoms;i++){
                for(j=1;j<neighbors[bond[k].affected_atom[i]].size();j++){
                    for(l=0;l<bond[k].affected_atom.size();l++){
                        check=0;
                        if(neighbors[bond[k].affected_atom[i]][j]==bond[k].affected_atom[l]){
                            check=1;
                            break;
                        }
                    }//loop over current affected atoms to check
                    if(check==0){
                        bond[k].affected_atom.push_back(neighbors[bond[k].affected_atom[i]][j]);
                        found=1;
                        bond[k].N_aff_atoms++;
                        start++;
                    }
                }//loop over neighbors of affected atom
            }//loop over last shell of affected_atoms
        }while(found==0);
    }//looping over rotatable bonds
}

void get_aff_atoms(int Nbonds, int N_rotatable_bonds, struct Bond *bond, int *ori_atom_id, int *tar_atom_id){
    //This function is used to find all the atoms that are affected by the rotation of a particular bond
    //It looks down the chain in the direction of atom1->atom2 to find all the atoms down the chain and its branches.

    int found, check,i,j,k,l;

    for(i=0;i<N_rotatable_bonds;i++){
        bond[i].N_aff_atoms=1;
        bond[i].affected_atom.resize(bond[i].N_aff_atoms);
        bond[i].affected_atom[0]=bond[i].atom2;
    }

    for(k=0;k<N_rotatable_bonds;k++){
    do{
        found=1;
        for(i=0;i<bond[k].N_aff_atoms;i++){
            for(j=0;j<Nbonds;j++){
                if((ori_atom_id[j]-1 == bond[k].affected_atom[i]) && (tar_atom_id[j]-1 != bond[k].atom1)){
                   check=1;
                   for(l=0;l<bond[k].N_aff_atoms;l++){
                        if(tar_atom_id[j]-1 == bond[k].affected_atom[l])
                            check=0;
                   }
                    if(check!=0){
                        found=0;
                        bond[k].affected_atom.push_back(1);
                        bond[k].affected_atom[bond[k].N_aff_atoms]=tar_atom_id[j]-1;
                        bond[k].N_aff_atoms++;
                    }
                }

                else if((tar_atom_id[j]-1 == bond[k].affected_atom[i]) && (ori_atom_id[j]-1 != bond[k].atom1)){
                   check=1;
                   for(l=0;l<bond[k].N_aff_atoms;l++){
                        if(ori_atom_id[j]-1 == bond[k].affected_atom[l])
                            check=0;
                   }
                    if(check!=0){
                        found=0;
                        bond[k].affected_atom.push_back(1);
                        bond[k].affected_atom[bond[k].N_aff_atoms]=ori_atom_id[j]-1;
                        bond[k].N_aff_atoms++;
                    }
                }
            }
        }
    }while(found==0);
    }
}

void get_internuclear_vector(double *bondvector, double *atom1, double *atom2){
    bondvector[0]=atom2[0]-atom1[0];
    bondvector[1]=atom2[1]-atom1[1];
    bondvector[2]=atom2[2]-atom1[2];
    double len=sqrt(bondvector[0]*bondvector[0]+bondvector[1]*bondvector[1]+bondvector[2]*bondvector[2]);
    bondvector[0]=bondvector[0]/len;
    bondvector[1]=bondvector[1]/len;
    bondvector[2]=bondvector[2]/len;
}

