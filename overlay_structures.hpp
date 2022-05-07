#include "REDOR_data.hpp"

void center_structure(int N_atoms, double (*xyz)[3]){
//This function will recenter the structure at the origin on the surface to enable optimal viewing of the overlays
    int i;
    double x_ave=0., y_ave=0.;

    for(i=0;i<N_atoms;i++){
        x_ave=x_ave+xyz[i][0];
        y_ave=y_ave+xyz[i][0];
    }
    x_ave=x_ave/N_atoms;
    y_ave=y_ave/N_atoms;

    for(i=0;i<N_atoms;i++){
        xyz[i][0]=xyz[i][0]-x_ave;
        xyz[i][0]=xyz[i][0]-x_ave;
    }
}

double rmsd(int N_atoms, double (*xyz_ref)[3], double (*xyz)[3]){
    //returns the rmsd between two structures
    //used to overlay two structures

    int i,j;
    double RMSD=0., dist;

    for(i=0;i<N_atoms;i++){
        dist=distance_calc(xyz[i],xyz_ref[i]);
        RMSD=RMSD + dist*dist;
    }
    RMSD=sqrt(RMSD/N_atoms);
    return RMSD;
}

double overlay_structures(int N_atoms, double (*xyz_ref)[3], double (*xyz)[3]){
    //This function will reorient structure xyz so as to minimize the RMSD with xyz_ref
    //The function also returns this RMSD which can be used to exclude certain structured
    //from an overlaid mol2/cif file.
    double xyz_temp[N_atoms][3];
    double RMSDm, RMSDp, RMSDo;
    int i,check_x=1, check_y=1, check_z=1, count=0;
    double RMSD_start, RMSD_end, RMSD_threshold=0.001;

    //The function is a nested do-while loop wherein the structure is translated along x, then y, and then rotated along z, in that order.
    //When a minimum is found in a given dimension it moves to the next dimension
    //Once the minimum RMSD is found the outer, Z, while loop exits
    do{//rotate along z
        RMSDo=rmsd(N_atoms,xyz_ref,xyz);
        RMSD_start=RMSDo;
        rotate_molecule_around_Z(N_atoms, xyz,  1.*Pi/180.);
        RMSDp=rmsd(N_atoms,xyz_ref,xyz);
        rotate_molecule_around_Z(N_atoms, xyz, -2.*Pi/180.);
        RMSDm=rmsd(N_atoms,xyz_ref,xyz);

        if(RMSDp<RMSDo){
            rotate_molecule_around_Z(N_atoms, xyz, 2.*Pi/180.);
            check_x=check_y=1;
        }
        else if(RMSDo<RMSDm){
            rotate_molecule_around_Z(N_atoms, xyz, 1.*Pi/180.);
            check_z=0;
        }
        else
            check_x=check_y=1;
        do{ //shift along y
            if(count>10){
                count=0;
                break;
            }
            RMSDo=rmsd(N_atoms,xyz_ref,xyz);
            translate_molecule_Y(N_atoms, xyz, 0.05);
            RMSDp=rmsd(N_atoms,xyz_ref,xyz);
            translate_molecule_Y(N_atoms, xyz, -0.10);
            RMSDm=rmsd(N_atoms,xyz_ref,xyz);
            count++;



            if(RMSDp<RMSDo){
                translate_molecule_Y(N_atoms, xyz, 0.10);
                check_x=check_z=1;
            }
            else if(RMSDo<RMSDm){
                translate_molecule_Y(N_atoms, xyz, 0.05);
                check_y=0;
            }
            else
                check_x=check_z=1;
            do{//shift along x
                if(count>10){
                    count=0;
                    break;
                }
                RMSDo=rmsd(N_atoms,xyz_ref,xyz);
                translate_molecule_X(N_atoms, xyz, 0.05);
                RMSDp=rmsd(N_atoms,xyz_ref,xyz);
                translate_molecule_X(N_atoms, xyz, -0.10);
                RMSDm=rmsd(N_atoms,xyz_ref,xyz);
                count++;

                if(RMSDp<RMSDo){
                    translate_molecule_X(N_atoms, xyz, 0.10);
                    check_y=check_z=1;
                }
                else if(RMSDo<RMSDm){
                    check_x=0;
                    translate_molecule_X(N_atoms, xyz, 0.05);
                }
                else
                    check_y=check_z=1;
            }while(check_x==1);
        }while(check_y==1);
        RMSD_end=rmsd(N_atoms,xyz_ref,xyz);
        if((RMSD_start-RMSD_end)<RMSD_threshold)
            break;

    }while(check_z==1);

    return RMSD_end;
}
