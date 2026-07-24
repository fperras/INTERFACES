#include "REDOR_data.hpp"
#include "3ZCW.h"

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

void find_images(int detected_index,int max_nrec, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &cell){
    //Function to find all the required symmetry-generated atoms from a unit cell for a given detected spin
    //Truncates the list to 97% of the sum of dipolar coupling squared, which gives a final 99% accuracy
    //Number of this is then increased by the inverse cubic root of the abundance of the recoupled spins
    int i,j,k,l;
    double abundance = REDOR.NA;

    //converting cell matrix to a cartesian form
    vector<double> unit_cell;// a,b,c, alpha, beta, gamma
    unit_cell.resize(6,1.);
    calc_cell_dimensions(cell,unit_cell);

    //Predicting an overestimate of the distance cutoff to build a supercell
    double distance = 10000000.;
    for(i=0;i<REDOR.recoupled.size();i++){
        double d_temp = sqrt(pow(xyz[REDOR.detected[detected_index]][0]-xyz[REDOR.recoupled[i]][0],2.) + pow(xyz[REDOR.detected[detected_index]][1]-xyz[REDOR.recoupled[i]][1],2.) + pow(xyz[REDOR.detected[detected_index]][2]-xyz[REDOR.recoupled[i]][2],2.));
        if(d_temp<distance)
            distance=d_temp;
    }
    double cutoff = distance * 3.5/cbrt(abundance);

    //Finding the bounds of the supercell that incorporates the cutoff calculated as above
    int nmax[3];
    for(i=0;i<3;i++){
        nmax[i]= (int)(cutoff/fabs(unit_cell[i])+1);
    }

    //variable declarations
    int images=0;
    distance=0.;
    vector< vector<double> > lattice_shift;
    lattice_shift.resize(1, vector<double>(3,0.));
    vector< vector<double> > n_vector;
    n_vector.resize(1, vector<double>(3,0.));

    //The vector dists will contain the structs Nvec which contain the lattice n multipliers, the atom index, and distance from detected atom
    vector< Nvec > dists;
    dists.resize((2*nmax[0]+1)*(2*nmax[1]+1)*(2*nmax[2]+1)*REDOR.recoupled.size());
    double D2_final=0.;

    //Here we populate the dists vector and calculate the final limit for the dipolar coupling squared (D2)
    for(i=-nmax[0];i<=nmax[0];i++){
        for(j=-nmax[1];j<=nmax[1];j++){
            for(k=-nmax[2];k<=nmax[2];k++){
                n_vector[0][0]=(double)i;
                n_vector[0][1]=(double)j;
                n_vector[0][2]=(double)k;
                frac_to_xyz(1,lattice_shift,n_vector,cell);

                for(l=0;l<REDOR.recoupled.size();l++){
                    distance = sqrt(pow(xyz[REDOR.detected[detected_index]][0]-xyz[REDOR.recoupled[l]][0]+lattice_shift[0][0],2.) + pow(xyz[REDOR.detected[detected_index]][1]-xyz[REDOR.recoupled[l]][1]+lattice_shift[0][1],2.) + pow(xyz[REDOR.detected[detected_index]][2]-xyz[REDOR.recoupled[l]][2]+lattice_shift[0][2],2.));
                    dists[images].nx=i;
                    dists[images].ny=j;
                    dists[images].nz=k;
                    dists[images].index=l;
                    dists[images].distance=distance;
                    images++;
                    D2_final+=pow(distance,-6.);
    }}}}

    //Sorting all images from closest to the detected atom to farthest
    qsort(dists.data(),dists.size(),sizeof(Nvec),compare_Nvec);

    //looping over images to find how many are needed to be 95% from teh final D2 value
    images=0;
    double D2=0.;
    double D2_cutoff=0.97;
    int keep=dists.size();
    for(i=0;i<dists.size();i++){
        images++;
        D2+=pow(dists[i].distance,-6.);
        if((D2/D2_final)>D2_cutoff){
            keep=images/abundance;
            break;
            }
    }

    if(keep>max_nrec)
        keep=max_nrec;

    //removing the excess recoupled atom images from the vector
    dists.resize(keep);

    //storing the results with the REDOR data
    REDOR.recoupled_index[detected_index].resize(dists.size());
    REDOR.nx[detected_index].resize(dists.size());
    REDOR.ny[detected_index].resize(dists.size());
    REDOR.nz[detected_index].resize(dists.size());

    for(i=0;i<dists.size();i++){
        REDOR.recoupled_index[detected_index][i]=dists[i].index;
        REDOR.nx[detected_index][i]=dists[i].nx;
        REDOR.ny[detected_index][i]=dists[i].ny;
        REDOR.nz[detected_index][i]=dists[i].nz;
    }
}

void find_all_images(int max_nrec, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &cell){
    //Finds all symmetry-generated atom positions required to converge a given REDOR dataset
    //from a unit cell specification, structure, and abundance.
    for(int i=0;i<REDOR.detected.size();i++){
        find_images(i,max_nrec,REDOR,xyz,cell);
    }
}

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

double REDOR_full(double time, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &cell){
    //This function returns the dephasing form a group of detected spins as an average
    //for intramolecular REDOR
	double DSS0=0.;
	double distance;
	int i,j,k, ndet=REDOR.detected.size(), nrec;//, nrec=REDOR.recoupled.size();

	//for finding symmetry images
	vector< vector<double> > lattice_shift;
    lattice_shift.resize(1, vector<double>(3,0.));
    vector< vector<double> > n_vector;
    n_vector.resize(1, vector<double>(3,0.));

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
    double x, y, z, xy, alphaD, betaD, DS;


    //calculation of the multispin REDOR or RESPDOR datapoint intensity
	for(i=0; i<ndet; i++){
	    nrec=REDOR.recoupled_index[i].size();
	    double D[nrec],saD[nrec],caD[nrec],sbD[nrec],cbD[nrec];
	    int rec;
        if(nrec==1){//uses the faster Bessel approach

            n_vector[0][0]=(double) REDOR.nx[i][0];
            n_vector[0][1]=(double) REDOR.ny[i][0];
            n_vector[0][2]=(double) REDOR.nz[i][0];
            frac_to_xyz(1,lattice_shift,n_vector,cell);

            rec=REDOR.recoupled[REDOR.recoupled_index[i][0]];

            x = xyz[REDOR.detected[i]][0]-xyz[rec][0]+lattice_shift[0][0];
            y = xyz[REDOR.detected[i]][1]-xyz[rec][1]+lattice_shift[0][1];
            z = xyz[REDOR.detected[i]][2]-xyz[rec][2]+lattice_shift[0][2];

            distance = sqrt(x*x+y*y+z*z);

            if(REDOR.type==1)//REDOR
                DSS0 = DSS0+  REDOR.NA*REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin)/ndet;

            else if(REDOR.type==2)//SEDOR
                    DSS0 = DSS0+  REDOR.NA*SEDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin)/ndet;

        }

        else{//multispin (nrec>1)
        //calculating polar angles
        for(j=0;j<nrec;j++){

            n_vector[0][0]=(double) REDOR.nx[i][j];
            n_vector[0][1]=(double) REDOR.ny[i][j];
            n_vector[0][2]=(double) REDOR.nz[i][j];

            frac_to_xyz(1,lattice_shift,n_vector,cell);

            rec=REDOR.recoupled[REDOR.recoupled_index[i][j]];

            x = xyz[REDOR.detected[i]][0]-xyz[rec][0]+lattice_shift[0][0];
            y = xyz[REDOR.detected[i]][1]-xyz[rec][1]+lattice_shift[0][1];
            z = xyz[REDOR.detected[i]][2]-xyz[rec][2]+lattice_shift[0][2];

            xy = sqrt(x*x+y*y);
            z = fabs(z);
            if(x>0.)
                alphaD=atan(y/x);
            else if(x<0.)
                alphaD=atan(y/x)+Pi;
            else
                alphaD=0.0;
            betaD=atan(xy/(z));
            D[j] = RDD(REDOR.RDD1A,sqrt(x*x+y*y+z*z))*REDOR.order_parameter;
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

	}

    return DSS0*REDOR.scaling_factor;
}

double calculate_curve_Chi2(REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &cell){
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
        DSS0_calc=REDOR_full(REDOR.tmix[i],REDOR,xyz,cell);
        X2 = X2 + pow((REDOR.DSS0[i]-DSS0_calc),2.) / pow((REDOR.DSS0[i]+REDOR.scaling_factor/10.),2.);
    }

    return X2;
}

void write_fits_meticulous(char *base_filename, const char *support_name, vector< REDOR_dataset > &REDOR, vector< vector< vector< vector<double> > > > &xyz, vector< vector<double> > &cell){
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
                DSS0=REDOR_full(time,REDOR[j],xyz[j][0],cell);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,REDOR[j],xyz[j][1],cell);
                fprintf(out,"%lf,",DSS0);
                DSS0=REDOR_full(time,REDOR[j],xyz[j][2],cell);
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

void precalculate_dephasing(vector<double> &calc_DSS0, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &cell){
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
        calc_DSS0[i]=REDOR_full(REDOR.tmix[i],REDOR,xyz,cell);
    }
    return;
}

double calculate_curve_Chi2_multi(double *curve_chi2, REDOR_dataset &REDOR, vector< vector<double> > &xyz, vector< vector<double> > &cell){
    //This function returns the Chi2 contributions from a single curve in a system containing multiple molecules
    //The function loops over site populations from 5% to 50 % in 5% increments and stores the curve_chi2 values
    //in an array of the same name. These values can then be compared to determine whether the site is statistically
    //significant or not. The dephasing levels determined from a prior structure determination are imported as
    //base_DSS0[Npoints].

    int i,j, Npoints = REDOR.DSS0.size();
    double DSS0_calc, weight, minimum=100000000.;
    vector<double> new_DSS0;
    new_DSS0.resize(Npoints, 0.);
    precalculate_dephasing(new_DSS0,REDOR,xyz,cell);

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

void write_fits_multi(int found_structures, char *base_filename, vector< REDOR_dataset > &REDOR, vector< vector< vector<double> > > &xyz, vector<double> &weights, vector< vector<double> > &cell){
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
                    DSS0 += weights[k]*REDOR_full(time,REDOR[j],xyz[k],cell);
                }
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }
    fclose(out);
}

