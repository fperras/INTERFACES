#include "REDOR_data.hpp"

double DSS0_full(double time, double scaling_factor, double order_parameter, vector<vector<int> > &REDOR_det_index, double (*xyz)[3], int curve_index, vector< vector< vector<double> > > DSS0_lib){
    //This function returns the dephasing form a group of detected spins as an average
    //for surface-to-atom REDOR

    int time_index = round((time/0.0002)*order_parameter - 1.);
	time_index = (time_index>249)*249 + (time_index<=249)*time_index;
	time_index = (time_index>-1)*time_index;
	double DSS0=0.;
	double distance;
	int i,d_index, ndet=REDOR_det_index[curve_index].size();

    for(i=0; i<ndet; i++){
        distance = xyz[REDOR_det_index[curve_index][i]][2];
        d_index=round(distance*10.+0.1001-1.);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0 = DSS0 + DSS0_lib[curve_index][time_index][d_index]/ndet;
    }

    return DSS0*scaling_factor;

}

double REDOR_full(double time, double scaling_factor, double order_parameter,double RDD1A, int spin, vector<vector<int> > &REDOR_det_index,vector<vector<int> > &REDOR_rec_index, double (*xyz)[3], int curve_index, vector< vector<double> > &REDORs){
    //This function returns the dephasing form a group of detected spins as an average
    //for intramolecular REDOR
	double DSS0=0.;
	double distance;
	int i, ndet=REDOR_det_index[curve_index].size();

    for(i=0; i<ndet; i++){
        distance = get_effective_distance(REDOR_det_index,REDOR_rec_index,xyz,curve_index,i);
        DSS0 = DSS0+  REDOR_DSS0(RDD(RDD1A,distance),time,order_parameter,spin,REDORs)/ndet;
    }

    return DSS0*scaling_factor;
}

double calculate_curve_Chi2(double scaling_factor, double order_parameter, vector<vector<int> > &REDOR_det_index,vector<vector<int> > &REDOR_rec_index, double (*xyz)[3], char (*element)[3], int curve_index, int curve_type, vector< vector<double> > &REDORs, vector< vector< vector<double> > > DSS0_lib, vector<double> &DSS0,vector<double> tmix){
    //This function returns the Chi2 contribution from a single curve.


    int i, Npoints = DSS0.size();
    double X2=0, DSS0_calc;
    double RDD1A=RDD_1A(element[REDOR_det_index[curve_index][0]],element[REDOR_rec_index[curve_index][0]]);
    int S=spin(element[REDOR_rec_index[curve_index][0]]);


    if(curve_type==0){//surface-to-atom
        for(i=0; i<Npoints; i++){
            DSS0_calc=DSS0_full(tmix[i],scaling_factor,order_parameter,REDOR_det_index,xyz,curve_index,DSS0_lib);
            X2 = X2 + pow((DSS0[i]-DSS0_calc),2.) / pow((DSS0[i]),2.);
        }
        return X2;
    }

    for(i=0; i<Npoints; i++){
        DSS0_calc=REDOR_full(tmix[i],scaling_factor,order_parameter,RDD1A,S,REDOR_det_index,REDOR_rec_index,xyz,curve_index,REDORs);
        X2 = X2 + pow((DSS0[i]-DSS0_calc),2.) / pow((DSS0[i]),2.);
    }
    return X2;
}

void write_fits_meticulous(int N_atoms,char *base_filename, int N_curves, const char *support_name, vector<vector<int> > &REDOR_det_index, vector<vector<int> > &REDOR_rec_index, double *xyz_it, char (*element)[3], char (*curve_filename)[120], double *scaling_factor, double *order_parameter, int *Nspins, int *curve_type){
        //This function writes out the REDOR curves form the best-fit structure as well as the range of dephasing for each curve
    //Data is stored in a CSV file
    //In the arrays listing the best-fin, minimum, and maximum distances and STDEV the indices have the following meaning:
    //0=Best Fit, 1=Smallest distance and 2=Largest distance from the surface
    FILE *error_file;
    char fits_filename[128];
    int i, j, k, l;
    double DSS0;
    FILE *out;
    double xyz[N_curves][3][N_atoms][3];

    for(i=0;i<N_curves;i++){
        for(j=0;j<3;j++){
            for(k=0;k<N_atoms;k++){
                for(l=0;l<3;l++){
                    xyz[i][j][k][l]=xyz_it[l+k*3+j*3*N_atoms+i*9*N_atoms];
                }
            }
        }
    }

    sprintf(fits_filename, "%s_REDOR_fits.csv", base_filename);

    remove(fits_filename);

    out=fopen(fits_filename,"w");

    if(out==NULL){
        error_file=fopen("Errors.txt","a");
        fprintf(error_file, "\nERROR: Could not write REDOR file '%s' for best fit structure\n", fits_filename);
        fclose(error_file);
    }

    vector< vector< vector<double> > > DSS0_lib;
    DSS0_lib.resize(N_curves, vector<vector<double> > (250, vector<double>(200,0.)));

    int all_types=0;
    for(j=0;j<N_curves;j++){
            all_types=all_types+curve_type[j];
            if(curve_type[j]==0)
                load_simulations(support_name, element[REDOR_det_index[j][0]], DSS0_lib[j]);
    }

    vector< vector<double> > REDORs;
    REDORs.resize(10000, vector<double>(9,0.));

    if(all_types!=0)
        generate_REDORs(REDORs);

    //writing the header to the CSV file
    fprintf(out,",");
    for(i=0;i<N_curves;i++){
        fprintf(out,",%s,,,",curve_filename[i]);
    }
    fprintf(out,"\ntime (s),,");
    for(i=0;i<N_curves;i++){
        fprintf(out,"best-fit,max,min,,");
    }
    fprintf(out,"\n");

    //writing out the dephasing values to the CSV file
    for(i=1;i<=250;i++){
        double time=i*0.0002;
        fprintf(out,"%lf,,",time);

        for(j=0;j<N_curves;j++){
            if(curve_type[j]==0){//surface curves
                DSS0 = DSS0_full(time,scaling_factor[j],order_parameter[j],REDOR_det_index,xyz[j][0],j,DSS0_lib);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_full(time,scaling_factor[j],order_parameter[j],REDOR_det_index,xyz[j][1],j,DSS0_lib);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_full(time,scaling_factor[j],order_parameter[j],REDOR_det_index,xyz[j][2],j,DSS0_lib);
                fprintf(out,"%lf,,",DSS0);
            }
            else{
                double RDD1A=RDD_1A(element[REDOR_det_index[j][0]],element[REDOR_rec_index[j][0]]);
                int S=spin(element[REDOR_rec_index[j][0]]);
                DSS0=REDOR_full(time,scaling_factor[j],order_parameter[j],RDD1A,S,REDOR_det_index,REDOR_rec_index,xyz[j][0],j,REDORs);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,scaling_factor[j],order_parameter[j],RDD1A,S,REDOR_det_index,REDOR_rec_index,xyz[j][1],j,REDORs);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,scaling_factor[j],order_parameter[j],RDD1A,S,REDOR_det_index,REDOR_rec_index,xyz[j][2],j,REDORs);
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }

    fclose(out);
}
