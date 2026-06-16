#include "REDOR_data.hpp"
#include "3ZCW.h"

double DSS0_full(double time, REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the dephasing form a group of detected spins as an average
    //for surface-to-atom REDOR

    int time_index = round((time/0.0001)*REDOR.order_parameter - 1.);
	time_index = (time_index>499)*499 + (time_index<=499)*time_index;
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

    return DSS0*REDOR.scaling_factor*REDOR.NA;

}

double SEDOR_coswDt(double sa, double ca, double sb, double cb, double sg, double cg, double time, double RDD, double sa2, double ca2, double sb2, double cb2, int spin) {
    //similar to the next function but for SEDOR data instead of REDOR
    //similarly assumes saturation of all transitions unless it is a spin-1/2 for the CT_sat keyword is used in the input
    double x = sb2*ca2;
    double y = sb2*sa2;
    double z = cb2;

    double cb3 = sb*ca*x+sb*sa*y+cb*z;

    double w = Pi*RDD*(3.*cb3*cb3-1.)*time;

    switch(spin){
        case 1 :
        return -cos(w);
        break;

        case 2 :
        return -1./3. -4./9.*cos(w)-2./9.*cos(2*w);
        break;

        case 3 :
        return -1./4. -3./8.*cos(w)-1./4.*cos(2*w)-1./8.*cos(3*w);
        break;

        case 4 :
        return -1./5. -8./25.*cos(w)-6./25.*cos(2*w)-4./25.*cos(3*w)-2./25.*cos(4*w);
        break;

        case 5 :
        return -1./6. - 5./18.*cos(w) - 2./9.*cos(2*w) - 1./6.*cos(3*w) - 1./9.*cos(4*w) - 1./18.*cos(5*w);
        break;

        case 6 :
        return -1./7. - 12./49.*cos(w) - 10./49.*cos(2*w) - 8./49.*cos(3*w) - 6./49.*cos(4*w) - 4./49.*cos(5*w) - 2./49.*cos(6*w);
        break;

        case 7 :
        return -1./8. - 14./64.*cos(w) - 12./64.*cos(2*w) - 10./64.*cos(3*w) - 8./64.*cos(4*w) - 6./64.*cos(5*w) - 4./64.*cos(6*w) - 2./64.*cos(7*w);
        break;

        case 8 :
        return -1./9. - 16./81.*cos(w) - 14./81.*cos(2*w) - 12./81.*cos(3*w) - 10./81.*cos(4*w) - 8./81.*cos(5*w) - 6./81.*cos(6*w) - 4./81.*cos(7*w) - 2./81.*cos(8*w);
        break;

        case 9 :
        return -1./10. - 9./50.*cos(w) - 8./50.*cos(2*w) - 7./50.*cos(3*w) - 6./50.*cos(4*w) - 5./50.*cos(5*w) - 4./50.*cos(6*w) - 3./50.*cos(7*w) - 2./50.*cos(8*w) - 1./50.*cos(9*w);
        break;
	   }

    return 0.0;
}

double coswDt(double sa, double ca, double sb, double cb, double sg, double cg, double time, double RDD, double sa2, double ca2, double sb2, double cb2, int spin) {
    //This is a function to calculate the cosine of the dipolar frequency
    //equations taken from JMR 127, 147-154 (1997) and PCCP 12, 9395-9405 (2010)

       double x = sb2*ca2;
       double y = sb2*sa2;
       double z = cb2;

       double niz = sb*ca*x + sb*sa*y + cb*z;
       double niy = (-cg*sa - cb*ca*sg)*x + (cg*ca - cb*sa*sg)*y + sg*sb*z;
       double s2bca = 2*niz*niy;

       double w = 2.828427125*time*(RDD*s2bca);

       switch(spin){
			case 1 :
			return -cos(w);
			break;

			case 2 :
			return -1./3. -4./9.*cos(w)-2./9.*cos(2*w);
			break;

			case 3 :
			return -1./4. -3./8.*cos(w)-1./4.*cos(2*w)-1./8.*cos(3*w);
			break;

			case 4 :
			return -1./5. -8./25.*cos(w)-6./25.*cos(2*w)-4./25.*cos(3*w)-2./25.*cos(4*w);
			break;

			case 5 :
			return -1./6. - 5./18.*cos(w) - 2./9.*cos(2*w) - 1./6.*cos(3*w) - 1./9.*cos(4*w) - 1./18.*cos(5*w);
			break;

			case 6 :
			return -1./7. - 12./49.*cos(w) - 10./49.*cos(2*w) - 8./49.*cos(3*w) - 6./49.*cos(4*w) - 4./49.*cos(5*w) - 2./49.*cos(6*w);
			break;

			case 7 :
			return -1./8. - 14./64.*cos(w) - 12./64.*cos(2*w) - 10./64.*cos(3*w) - 8./64.*cos(4*w) - 6./64.*cos(5*w) - 4./64.*cos(6*w) - 2./64.*cos(7*w);
			break;

			case 8 :
			return -1./9. - 16./81.*cos(w) - 14./81.*cos(2*w) - 12./81.*cos(3*w) - 10./81.*cos(4*w) - 8./81.*cos(5*w) - 6./81.*cos(6*w) - 4./81.*cos(7*w) - 2./81.*cos(8*w);
			break;

			case 9 :
			return -1./10. - 9./50.*cos(w) - 8./50.*cos(2*w) - 7./50.*cos(3*w) - 6./50.*cos(4*w) - 5./50.*cos(5*w) - 4./50.*cos(6*w) - 3./50.*cos(7*w) - 2./50.*cos(8*w) - 1./50.*cos(9*w);
			break;
	   }
    return 0.0;
}

double REDOR_full_mean(double time, REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the dephasing form a group of detected spins as an average
    //for intramolecular REDOR
	double DSS0=0.;
	double distance;
	int i, ndet=REDOR.detected.size();

	if(REDOR.type==1){//REDOR
    for(i=0; i<ndet; i++){
        distance = get_effective_distance(REDOR,xyz,i);
        DSS0 = DSS0+  REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin)/ndet;
    }
	}

	else if(REDOR.type==2){//SEDOR
    for(i=0; i<ndet; i++){
        distance = get_effective_distance(REDOR,xyz,i);
        DSS0 = DSS0+  SEDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin)/ndet;
    }
	}

    return DSS0*REDOR.scaling_factor*REDOR.NA;
}

double REDOR_full(double time, REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the dephasing form a group of detected spins as an average
    //for intramolecular REDOR
	double DSS0=0.;
	double distance;
	int i,j,k, ndet=REDOR.detected.size(), nrec=REDOR.recoupled.size();

	if(nrec==1){//uses the faster Bessel approach
        if(REDOR.type==1){//REDOR
            for(i=0; i<ndet; i++){
                distance = get_effective_distance(REDOR,xyz,i);
                DSS0 = DSS0+  REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin)/ndet;
            }
            return DSS0*REDOR.scaling_factor*REDOR.NA;
        }
        else if(REDOR.type==2){//SEDOR
            for(i=0; i<ndet; i++){
                distance = get_effective_distance(REDOR,xyz,i);
                DSS0 = DSS0+  SEDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin)/ndet;
            }
            return DSS0*REDOR.scaling_factor*REDOR.NA;
        }
	}

	//multi-spin numerical integration
	//3-angle ZCW grids
    int N_orient=REDOR.ZCWg;
    vector<double> sa,ca,sb,cb,sg,cg, intensity;
    sa.resize(N_orient,0.);
    ca.resize(N_orient,0.);
    sb.resize(N_orient,0.);
    cb.resize(N_orient,0.);
    sg.resize(N_orient,0.);
    cg.resize(N_orient,0.);
    intensity.resize(N_orient,0.);
    ZCWt(sa,ca,sb,cb,sg,cg,intensity,N_orient);
    double x[nrec], y[nrec], z[nrec], xy[nrec], alphaD, betaD, D[nrec], DS;
    double saD[nrec],caD[nrec],sbD[nrec],cbD[nrec];

    //calculation of the multispin REDOR or RESPDOR datapoint intensity
	for(i=0; i<ndet; i++){
        //calculating polar angles
        for(j=0;j<nrec;j++){
            x[j] = xyz[REDOR.detected[i]][0]-xyz[REDOR.recoupled[j]][0];
            y[j] = xyz[REDOR.detected[i]][1]-xyz[REDOR.recoupled[j]][1];
            z[j] = xyz[REDOR.detected[i]][2]-xyz[REDOR.recoupled[j]][2];
            xy[j] = sqrt(x[j]*x[j]+y[j]*y[j]);
            z[j] = fabs(z[j]);
            if(x[j]>0.)
                alphaD=atan(y[j]/x[j]);
            else if(x[j]<0.)
                alphaD=atan(y[j]/x[j])+Pi;
            else
                alphaD=0.0;
            betaD=atan(xy[j]/(z[j]));
            D[j] = RDD(REDOR.RDD1A,sqrt(x[j]*x[j]+y[j]*y[j]+z[j]*z[j]))*REDOR.order_parameter;
            sincos(alphaD,&saD[j],&caD[j]);
            sincos(betaD,&sbD[j],&cbD[j]);
        }
        //dephasing calculation
        if(REDOR.type==1){//REDOR
            for(k=0; k<N_orient; k++){
                DS=1.;
                for(j=0;j<nrec;j++){
                    DS=DS*(1.-REDOR.NA*(1.-    (-coswDt(sa[k], ca[k], sb[k],cb[k], sg[k], cg[k], time, D[j], saD[j], caD[j], sbD[j], cbD[j], REDOR.spin))));
                }
                DSS0 += (1.-DS)*intensity[k]/ndet;
            }
        }
        if(REDOR.type==2){//SEDOR
            for(k=0; k<N_orient; k++){
                DS=1.;
                for(j=0;j<nrec;j++){
                    DS=DS*(1.-REDOR.NA*(1.-    (-SEDOR_coswDt(sa[k], ca[k], sb[k],cb[k], sg[k], cg[k], time, D[j], saD[j], caD[j], sbD[j], cbD[j], REDOR.spin))));
                }
                DSS0 += (1.-DS)*intensity[k]/ndet;

            }
        }
	}

    return DSS0*REDOR.scaling_factor;
}

double calculate_curve_Chi2(REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the Chi2 contribution from a single curve.

    int i, Npoints = REDOR.DSS0.size();
    double X2=0, DSS0_calc;

    if(REDOR.type==0){//surface-to-atom
        for(i=0; i<Npoints; i++){
            DSS0_calc=DSS0_full(REDOR.tmix[i],REDOR,xyz);
            X2 = X2 + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]+REDOR.scaling_factor/10.),2.);
        }
        return X2;
    }

    for(i=0; i<Npoints; i++){
        DSS0_calc=REDOR_full(REDOR.tmix[i],REDOR,xyz);
        X2 = X2 + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]+REDOR.scaling_factor/10.),2.);
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
    for(i=1;i<=500;i++){
        double time=i*0.0001;
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
                DSS0=REDOR_full(time,REDOR[j],xyz[j][0]);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,REDOR[j],xyz[j][1]);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,REDOR[j],xyz[j][2]);
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

void precalculate_dephasing(vector<double> &calc_DSS0, REDOR_dataset &REDOR, vector< vector<double> > &xyz){
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
        calc_DSS0[i]=REDOR_full(REDOR.tmix[i],REDOR,xyz);
    }
    return;
}

double calculate_curve_Chi2_multi(double *curve_chi2, REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the Chi2 contributions from a single curve in a system containing multiple molecules
    //The function loops over site populations from 5% to 50 % in 5% increments and stores the curve_chi2 values
    //in an array of the same name. These values can then be compared to determine whether the site is statistically
    //significant or not. The dephasing levels determined from a prior structure determination are imported as
    //base_DSS0[Npoints].

    int i,j, Npoints = REDOR.DSS0.size();
    double DSS0_calc, weight, minimum=100000000.;
    vector<double> new_DSS0;
    new_DSS0.resize(Npoints, 0.);
    precalculate_dephasing(new_DSS0,REDOR,xyz);

    if(REDOR.type==0){//surface-to-atom
        for(j=0;j<10;j++){
            weight = 0.05 + 0.05*j;
            curve_chi2[j]=0.0;
            for(i=0; i<Npoints; i++){
                DSS0_calc = weight*new_DSS0[i] + (1.-weight)*REDOR.DSS0sim_prev[i];
                curve_chi2[j] = curve_chi2[j] + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]+REDOR.scaling_factor/10.),2.);
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
            curve_chi2[j] = curve_chi2[j] + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]+REDOR.scaling_factor/10.),2.);
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
    for(i=1;i<=500;i++){
        double time=i*0.0001;
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
                    DSS0 += weights[k]*REDOR_full(time,REDOR[j],xyz[k]);
                }
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

