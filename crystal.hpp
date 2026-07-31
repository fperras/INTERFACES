using namespace std;
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <string>

struct fragments{
    int N_fragments;
    vector< vector<int> > frag_indices;  //indices of each atoms in a given fragment
    vector< vector<double> > mol_xyz; //xyz coordinates of the center of each fragment
    vector< vector<double> > mol_frac; //fractional coordinates of the center of each fragment
};

void xyz_to_frac(int N_atoms,const std::vector<std::vector<double>>& xyz,std::vector<std::vector<double>>& frac,const std::vector<std::vector<double>>& cell)
{
    for (int i = 0; i < N_atoms; i++) {
        frac[i][2] = xyz[i][2] / cell[2][2];
        frac[i][1] = xyz[i][1] / cell[1][1];
        frac[i][0] = (xyz[i][0] - cell[2][0] * frac[i][2]) / cell[0][0];
    }
}

void frac_to_xyz(int N_atoms, vector< vector<double> > &xyz, vector< vector<double> > &frac, vector< vector<double> > &cell){
    int i;

    for (int i = 0; i < N_atoms; i++) {
        xyz[i][0] = cell[0][0]*frac[i][0] + cell[1][0]*frac[i][1] + cell[2][0]*frac[i][2];
        xyz[i][1] = cell[0][1]*frac[i][0] + cell[1][1]*frac[i][1] + cell[2][1]*frac[i][2];
        xyz[i][2] = cell[0][2]*frac[i][0] + cell[1][2]*frac[i][1] + cell[2][2]*frac[i][2];
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

struct Nvec{
    //variable used for sorting by distance with compare_Nvec()
    int nx,ny,nz;
    int index;
    double distance;
} ;

int compare_Nvec(const void *a, const void *b)
{
    //sorting function used in find_images() to sort symmetry-generated atoms by distance to a central detected atom
    const Nvec *A = (const Nvec *)a;
    const Nvec *B = (const Nvec *)b;

    if (A->distance < B->distance) return -1;
    if (A->distance > B->distance) return 1;
    return 0;
}

void find_extra_atoms(vector<int> &extra_atoms, vector< vector<double> > &xyz, vector< vector<double> > &cell){
    //This function creates a list of the atom indices from symmetyr-generated atoms that were included in the mol2 file
    //This list is used to prevent double-counting of some atoms with sym_gen or density calculations.
    int i,j;

    //conversion to fractional coordinates
    vector< vector<double> > frac;
    frac.resize(xyz.size(), vector<double>(3,0.));
    xyz_to_frac(xyz.size(),xyz,frac,cell);

    //finding symmetry-related atoms on cell edges
    for(i=0;i<xyz.size();i++){
        for(j=0;j<xyz.size();j++){
            if(fabs(frac[i][0]-(frac[j][0]-1.))<1e-6){
                if((fabs(frac[i][1]-frac[j][1])<1e-6)&&(fabs(frac[i][2]-frac[j][2])<1e-6)){
                    extra_atoms.push_back(i);
                    break;
                }
            }
            else if(fabs(frac[i][1]-(frac[j][1]-1.))<1e-6){
                if((fabs(frac[i][0]-frac[j][0])<1e-6)&&(fabs(frac[i][2]-frac[j][2])<1e-6)){
                    extra_atoms.push_back(i);
                    break;
                }
            }
            else if(fabs(frac[i][2]-(frac[j][2]-1.))<1e-6){
                if((fabs(frac[i][0]-frac[j][0])<1e-6)&&(fabs(frac[i][1]-frac[j][1])<1e-6)){
                    extra_atoms.push_back(i);
                    break;
                }
            }
        }
    }
    sort(extra_atoms.rbegin(),extra_atoms.rend());
}

double atomic_mass(const char* element){
    //This function returns the atomic mass of an element.
    char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    static const std::unordered_map<std::string,double> mass = {
        {"H", 1.008}, {"He", 4.003},
        {"Li", 6.94}, {"Be", 9.012}, {"B", 10.81}, {"C", 12.01}, {"N", 14.01},
        {"O", 16.00}, {"F", 19.00}, {"Ne", 20.18},

        {"Na", 22.99}, {"Mg", 24.31}, {"Al", 26.98}, {"Si", 28.09}, {"P", 30.97},
        {"S", 32.06}, {"Cl", 35.45}, {"Ar", 39.95},

        {"K", 39.10}, {"Ca", 40.08}, {"Sc", 44.96}, {"Ti", 47.87}, {"V", 50.94},
        {"Cr", 52.00}, {"Mn", 54.94}, {"Fe", 55.85}, {"Co", 58.93}, {"Ni", 58.69},
        {"Cu", 63.55}, {"Zn", 65.38},

        {"Ga", 69.72}, {"Ge", 72.63}, {"As", 74.92}, {"Se", 78.97}, {"Br", 79.90},
        {"Kr", 83.80},

        {"Rb", 85.47}, {"Sr", 87.62}, {"Y", 88.91}, {"Zr", 91.22}, {"Nb", 92.91},
        {"Mo", 95.95}, {"Tc", 98.00}, {"Ru", 101.07}, {"Rh", 102.91}, {"Pd", 106.42},
        {"Ag", 107.87}, {"Cd", 112.41},

        {"In", 114.82}, {"Sn", 118.71}, {"Sb", 121.76}, {"Te", 127.60}, {"I", 126.90},
        {"Xe", 131.29},

        {"Cs", 132.91}, {"Ba", 137.33}, {"La", 138.91}, {"Ce", 140.12}, {"Pr", 140.91},
        {"Nd", 144.24}, {"Pm", 145.00}, {"Sm", 150.36}, {"Eu", 151.96}, {"Gd", 157.25},
        {"Tb", 158.93}, {"Dy", 162.50}, {"Ho", 164.93}, {"Er", 167.26}, {"Tm", 168.93},
        {"Yb", 173.05}, {"Lu", 174.97},

        {"Hf", 178.49}, {"Ta", 180.95}, {"W", 183.84}, {"Re", 186.21}, {"Os", 190.23},
        {"Ir", 192.22}, {"Pt", 195.08}, {"Au", 196.97}, {"Hg", 200.59},

        {"Tl", 204.38}, {"Pb", 207.2}, {"Bi", 208.98}, {"Po", 209.00}, {"At", 210.00},
        {"Rn", 222.00},

        {"Fr", 223.00}, {"Ra", 226.00}, {"Ac", 227.00}, {"Th", 232.04}, {"Pa", 231.04},
        {"U", 238.03}, {"Np", 237.00}, {"Pu", 244.00}, {"Am", 243.00}, {"Cm", 247.00},
        {"Bk", 247.00}, {"Cf", 251.00}, {"Es", 252.00}, {"Fm", 257.00}, {"Md", 258.00},
        {"No", 259.00}, {"Lr", 262.00},

        {"Rf", 267.00}, {"Db", 270.00}, {"Sg", 271.00}, {"Bh", 270.00}, {"Hs", 277.00},
        {"Mt", 276.00}, {"Ds", 281.00}, {"Rg", 280.00}, {"Cn", 285.00}, {"Nh", 284.00},
        {"Fl", 289.00}, {"Mc", 288.00}, {"Lv", 293.00}, {"Ts", 294.00}, {"Og", 294.00}
    };

    auto it = mass.find(element);
    if(it != mass.end())
        return it->second;

    error_file=fopen(error_filename,"a");
    fprintf(error_file, "\nERROR: Mass for '%s' is unknown\n", element);
    fclose(error_file);
    exit(1);
}

double density(int N_atoms,char (*element)[3], vector<double> &unit_cell, const vector<int> &extra_atoms){
    //This function returns the density of a crystalline solid using the atomic masses and
    //the cell volume. symmetry-related extra atoms are removed from the calculation.
    int i,j;
    double cell_mass=0.;
    for(int i = 0; i < N_atoms; i++){
        if(std::find(extra_atoms.begin(), extra_atoms.end(), i) == extra_atoms.end()){
            cell_mass += atomic_mass(element[i]);
        }
    }

    return 1.660539*cell_mass/calc_cell_volume(unit_cell);
}

