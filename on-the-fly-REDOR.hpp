#include "REDOR_data.hpp"

double DSS0_full(double time, REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the dephasing form a group of detected spins as an average
    //for surface-to-atom REDOR

    int time_index = round((time/0.0002)*REDOR.order_parameter - 1.);
	time_index = (time_index>249)*249 + (time_index<=249)*time_index;
	time_index = (time_index>-1)*time_index;
	double DSS0=0.;
	double distance;
	int i,d_index, ndet=REDOR.detected.size();

    for(i=0; i<ndet; i++){
        distance = xyz[REDOR.detected[i]][2];
        d_index=round(distance*10.+0.1001-1.);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0 = DSS0 + REDOR.DSS0_lib[time_index][d_index]/ndet;
    }

    return DSS0*REDOR.scaling_factor;

}

double REDOR_full(double time, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &REDORs){
    //This function returns the dephasing form a group of detected spins as an average
    //for intramolecular REDOR
	double DSS0=0.;
	double distance;
	int i, ndet=REDOR.detected.size();

    for(i=0; i<ndet; i++){
        distance = get_effective_distance(REDOR,xyz,i);
        DSS0 = DSS0+  REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin,REDORs)/ndet;
    }

    return DSS0*REDOR.scaling_factor;
}

double calculate_curve_Chi2(REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &REDORs){
    //This function returns the Chi2 contribution from a single curve.

    int i, Npoints = REDOR.DSS0.size();
    double X2=0, DSS0_calc;

    if(REDOR.type==0){//surface-to-atom
        for(i=0; i<Npoints; i++){
            DSS0_calc=DSS0_full(REDOR.tmix[i],REDOR,xyz);
            X2 = X2 + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]),2.);
        }
        return X2;
    }

    for(i=0; i<Npoints; i++){
        DSS0_calc=REDOR_full(REDOR.tmix[i],REDOR,xyz,REDORs);
        X2 = X2 + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]),2.);
    }
    return X2;
}

void write_fits_meticulous(char *base_filename, const char *support_name, vector< REDOR_dataset > &REDOR, vector< vector< vector< vector<double> > > > &xyz){
    //This function writes out the REDOR curves form the best-fit structure as well as the range of dephasing for each curve
    //Data is stored in a CSV file
    //In the arrays listing the best-fin, minimum, and maximum distances and STDEV the indices have the following meaning:
    //0=Best Fit, 1=Smallest distance and 2=Largest distance from the surface
    FILE *error_file;
    char fits_filename[128];
    int i, j, k, l;
    double DSS0;
    FILE *out;

    sprintf(fits_filename, "%s_REDOR_fits.csv", base_filename);
    remove(fits_filename);
    out=fopen(fits_filename,"w");

    if(out==NULL){
        error_file=fopen("Errors.txt","a");
        fprintf(error_file, "\nERROR: Could not write REDOR file '%s' for best fit structure\n", fits_filename);
        fclose(error_file);
    }

    vector< vector<double> > REDORs;
    REDORs.resize(10000, vector<double>(9,0.));
    generate_REDORs(REDORs);

    //writing the header to the CSV file
    fprintf(out,",");
    for(i=0;i<REDOR.size();i++){
        fprintf(out,",%s,,,",REDOR[i].filename);
    }
    fprintf(out,"\ntime (s),,");
    for(i=0;i<REDOR.size();i++){
        fprintf(out,"best-fit,max,min,,");
    }
    fprintf(out,"\n");

    //writing out the dephasing values to the CSV file
    for(i=1;i<=250;i++){
        double time=i*0.0002;
        fprintf(out,"%lf,,",time);

        for(j=0;j<REDOR.size();j++){
            if(REDOR[j].type==0){//surface curves
                DSS0 = DSS0_full(time,REDOR[j],xyz[j][0]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_full(time,REDOR[j],xyz[j][1]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_full(time,REDOR[j],xyz[j][2]);
                fprintf(out,"%lf,,",DSS0);
            }
            else{
                DSS0=REDOR_full(time,REDOR[j],xyz[j][0],REDORs);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,REDOR[j],xyz[j][1],REDORs);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,REDOR[j],xyz[j][2],REDORs);
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

void precalculate_dephasing(vector<double> &calc_DSS0, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &REDORs){
    //This function calculated the dephasing levels for a given set of xyz coordinates at the experimentally-specified recoupling times
    //Designed to be used with the calculate_curve_Chi2_multi() function.

    int i,j, Npoints = REDOR.DSS0.size();

    if(REDOR.type==0){//surface-to-atom
        for(i=0; i<Npoints; i++){
            calc_DSS0[i]=DSS0_full(REDOR.tmix[i],REDOR,xyz);
        }
        return;
    }

    for(i=0; i<Npoints; i++){
        calc_DSS0[i]=REDOR_full(REDOR.tmix[i],REDOR,xyz,REDORs);
    }
    return;
}

double calculate_curve_Chi2_multi(double *curve_chi2, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &REDORs){
    //This function returns the Chi2 contributions from a single curve in a system containing multiple molecules
    //The function loops over site populations from 5% to 50 % in 5% increments and stores the curve_chi2 values
    //in an array of the same name. These values can then be compared to determine whether the site is statistically
    //significant or not. The dephasing levels determined from a prior structure determination are imported as
    //base_DSS0[Npoints].

    int i,j, Npoints = REDOR.DSS0.size();
    double DSS0_calc, weight, minimum=100000000.;
    vector<double> new_DSS0;
    new_DSS0.resize(Npoints, 0.);
    precalculate_dephasing(new_DSS0,REDOR,xyz,REDORs);

    if(REDOR.type==0){//surface-to-atom
        for(j=0;j<10;j++){
            weight = 0.05 + 0.05*j;
            curve_chi2[j]=0.0;
            for(i=0; i<Npoints; i++){
                DSS0_calc = weight*new_DSS0[i] + (1.-weight)*REDOR.DSS0sim_prev[i];
                curve_chi2[j] = curve_chi2[j] + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]),2.);
            }
            minimum = (minimum < curve_chi2[j])*minimum + (minimum > curve_chi2[j])*curve_chi2[j];
        }
        return minimum;
    }

    for(j=0;j<10;j++){
        weight = 0.05 + 0.05*j;
        curve_chi2[j]=0.0;
        for(i=0; i<Npoints; i++){
            DSS0_calc = weight*new_DSS0[i] + (1.-weight)*REDOR.DSS0sim_prev[i];
            curve_chi2[j] = curve_chi2[j] + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]),2.);
        }
        minimum = (minimum < curve_chi2[j])*minimum + (minimum > curve_chi2[j])*curve_chi2[j];
    }
    return minimum;
}

void write_fits_multi(int found_structures, char *base_filename, vector< REDOR_dataset > &REDOR, vector< vector< vector<double> > > &xyz, vector<double> &weights){
    //This function writes out the REDOR curves for the multi-model fit
    //Data is stored in a CSV file
    //In the arrays listing the best-fin, minimum, and maximum distances and STDEV the indices have the following meaning:
    //0=Best Fit, 1=Smallest distance and 2=Largest distance from the surface
    FILE *error_file;
    char fits_filename[128];
    int i, j, k, l;
    double DSS0;
    FILE *out;

    sprintf(fits_filename, "%s_multi-site_REDOR_fits.csv", base_filename);
    remove(fits_filename);
    out=fopen(fits_filename,"w");

    if(out==NULL){
        error_file=fopen("Errors.txt","a");
        fprintf(error_file, "\nERROR: Could not write REDOR file '%s' for best fit structure\n", fits_filename);
        fclose(error_file);
    }

    vector< vector<double> > REDORs;
    REDORs.resize(10000, vector<double>(9,0.));

    generate_REDORs(REDORs);

    //writing the header to the CSV file
    fprintf(out,",");
    for(i=0;i<REDOR.size();i++){
        fprintf(out,",%s,",REDOR[i].filename);
    }
    fprintf(out,"\ntime (s),,");
    for(i=0;i<REDOR.size();i++){
        fprintf(out,"best-fit,,");
    }
    fprintf(out,"\n");

    //writing out the dephasing values to the CSV file
    for(i=1;i<=250;i++){
        double time=i*0.0002;
        fprintf(out,"%lf,,",time);

        for(j=0;j<REDOR.size();j++){
            if(REDOR[j].type==0){//surface curves
                DSS0=0.;
                for(k=0;k<found_structures;k++){
                    DSS0 += weights[k]*DSS0_full(time,REDOR[j],xyz[k]);
                }
                fprintf(out,"%lf,,",DSS0);
            }
            else{

                DSS0=0.;
                for(k=0;k<found_structures;k++){
                    DSS0 += weights[k]*REDOR_full(time,REDOR[j],xyz[k],REDORs);
                }
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

