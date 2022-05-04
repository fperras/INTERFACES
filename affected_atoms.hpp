#include "mol2_files.hpp"

struct Bond{
    int atom1;
    int atom2;
    int N_aff_atoms;
    int N_steps;
    vector<int> affected_atom;
    int mod;
};

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



