#include "statistics.hpp"
#include "constraints.hpp"

struct REDOR_dataset{
    vector<int> detected;
    vector<int> recoupled;
    int type;
    int spin;
    char rec_element[3];
    char det_element[3];
    double RDD1A;
    double scaling_factor;
    double order_parameter;
    char filename[120];
    double chi2_min;
    double chi2_max;
    double chi2_best;

    //These arrays to store the average and stdev distances of the (0)best fit structure, (1)Smallest distance and (2)Largest distance from the surface
    int d_index[3];
    int std_index[3];

    vector<double> DSS0;
    vector<double> tmix;
    vector<double> DSS0sim_prev;
    vector<double> DSS0sim_new;

    vector< vector<double> > X2;
    vector< vector<double> > DSS0_lib;
};

double get_effective_distance(REDOR_dataset &REDOR, vector< vector<double> > &xyz, int det_index){
    int i;
    double Deff=0., dist;

    for(i=0; i<REDOR.recoupled.size(); i++){
        dist = distance_calc(xyz[REDOR.detected[det_index]],xyz[REDOR.recoupled[i]]);
        Deff=Deff+1./pow(dist,6.);
    }
    return pow(1./Deff,1./6.);
}

double get_average_distance(REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the average vertical distance for the atoms that contribute to a given curve (curve_index)

    double d_sum=0.;
    int j;

    if(REDOR.type==0){

        for(j=0; j<REDOR.detected.size(); j++){
            d_sum = d_sum + xyz[REDOR.detected[j]][2];
        }
        return d_sum/REDOR.detected.size();
    }

    //will return the distance corresponding to the effective recoupled dipolar coupling in multispin case
    for(j=0; j<REDOR.detected.size(); j++){
        d_sum = d_sum + get_effective_distance(REDOR,xyz,j);
    }
    return d_sum/REDOR.detected.size();
}

double get_STDEV(REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the standard deviation of the vertical distance for the atoms that contribute to a given curve (curve_index)
    double d_STDEV_numerator=0.;
    int j;
    double d_avg = get_average_distance(REDOR,xyz);

    if(REDOR.type==0){
        for(j=0; j<REDOR.detected.size(); j++){
            d_STDEV_numerator = d_STDEV_numerator + pow(xyz[REDOR.detected[j]][2] - d_avg,2.);
        }
        return sqrt(d_STDEV_numerator/REDOR.detected.size());
    }

    for(j=0; j<REDOR.detected.size(); j++){
        d_STDEV_numerator = d_STDEV_numerator + pow(get_effective_distance(REDOR,xyz,j) - d_avg,2.);
    }
    return sqrt(d_STDEV_numerator/REDOR.detected.size());

}

int get_distance_index(REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the average vertical distance for the atoms that contribute to a given curve (curve_index)
    //result is given as an index in the chi squared table (distance*10)
        double distance = get_average_distance(REDOR, xyz);
        int d_index=round(distance*10.+0.1001-1.);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        return d_index;
}

int get_STDEV_index(REDOR_dataset &REDOR, vector< vector<double> > &xyz){
    //This function returns the standard deviation of the vertical distance for the atoms that contribute to a given curve (curve_index)
    //result is given as an index in the chi squared table (STDEV*10)
    double d_std = get_STDEV(REDOR, xyz);
    int std_index = round(d_std*10+0.001);
	std_index = (std_index>100)*100 + (std_index<=100)*std_index;
	return std_index;
}

int load_simulations(const char *support_name, REDOR_dataset &REDOR){
    //This function loads the simulated REDOR results for a given support-nuclear combination
	//Data are stored as DSS0_lib[250][200] with the first index corresponding to time and the second distance
	//DSS0_lib is organized as DSS0_lib[time][distance] with the distance indices running from 0.1 to 20 A (indices of 0 to 199)
    //and the time time running from 0.0002 s to 0.05 s in 0.0002 s increments.
	char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

	char filename[32];
	FILE *fp;
	int i,j;
	double time;

	sprintf(filename,"%s_%s.txt",support_name,REDOR.det_element);
	fp=fopen(filename,"r");

	if(fp == NULL){
		error_file=fopen(error_filename,"a");
		fprintf(error_file, "\nERROR: Reference file '%s' could not be found\n", filename);
		fclose(error_file);
        exit(1);
	}

	for(i=0; i<250;i++){
		fscanf(fp,"%lf",&time);
		for(j=0;j<200;j++){
			fscanf(fp,"%lf",&REDOR.DSS0_lib[i][j]);
		}
	}
	fclose(fp);
	return 0;
}

double DSS0_pred(double time, double d_mean, double d_std, vector< vector<double> > &gaussians, REDOR_dataset &REDOR){
	//This function is used to calculate the REDOR dephasing given a gaussian distribution of distances with d_mean and d_std
	//DSS0_lib is organized as DSS0_lib[time][distance] with the distance indices running from 0.1 to 20 A (indices of 0 to 199)
	//and the time time running from 0.0002 s to 0.05 s in 0.0002 s increments.

	int i, time_index, d_index;
	double DSS0 = 0., distance;

	time_index = round((time/0.0002)*REDOR.order_parameter - 1.);
	time_index = (time_index>249)*249 + (time_index<=249)*time_index;
	time_index = (time_index>-1)*time_index;

	if(REDOR.detected.size()==2){
        distance=d_mean-d_std;
        d_index=round(distance*10.+0.1001-1.);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0= DSS0 + 0.5*REDOR.DSS0_lib[time_index][d_index];

        distance=d_mean+d_std;
        d_index=round(distance*10.+0.1001-1.);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0= DSS0 + 0.5*REDOR.DSS0_lib[time_index][d_index];
	}

	else{
        for(i=0;i<200;i++){
            distance=0.1+i*0.1;
            DSS0= DSS0 + d_prob(distance, d_mean, d_std, gaussians)*REDOR.DSS0_lib[time_index][i];
        }
	}
	return DSS0*REDOR.scaling_factor;
}

int find_Npoints(const char *filename){
    //Returns the number of points in a REDOR curve file by counting the non-empty lines
	char error_filename[128], buffer[256];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file, *fp;
	int Npoints=0;

	fp=fopen(filename,"r");

	if(fp == NULL){
        error_file=fopen(error_filename,"a");
		fprintf(error_file, "\nERROR: could not find Npoints in file %s\n", filename);
        fclose(error_file);
        exit(1);
	}

	while(fgets(buffer,256,fp)!=NULL){
        if(strlen(buffer) > 2)
            Npoints++;
	}

	fclose(fp);
	return Npoints;
}

int load_exp_curve(REDOR_dataset &REDOR){
    //This function is used to load a single experimental REDOR curve.
	//experimental curves are a vector of  the form DSS0[t_index]
	//these files have the format "time(s)   DS/S0\n"
	//a similar vector exists to save the recoupling times (tmix).
	char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

	char  buffer[64];
	FILE *fp;
	int i,Npoints;

	Npoints=find_Npoints(REDOR.filename);
	fp=fopen(REDOR.filename,"r");

	if(fp == NULL){
		error_file=fopen(error_filename,"a");
		fprintf(error_file, "\nERROR: could not find exp file '%s'\n", REDOR.filename);
		fclose(error_file);
        exit(1);
	}

	for(i=0;i<Npoints;i++){
		if(fgets(buffer,64,fp)!=NULL){
		sscanf(buffer,"%lf %lf",&REDOR.tmix[i],&REDOR.DSS0[i]);
		if(REDOR.DSS0[i]<=0.){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: illegal REDOR intensity of %lf found in %s\n",REDOR.DSS0[i], REDOR.filename);
            fclose(error_file);
            exit(1);
		}
		else if(REDOR.DSS0[i]<0.05){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nWARNING: Very low dephasing value of %lf in %s can lead to over fitting\nYou may obtain a better fit by eliminating this data point.\n",REDOR.DSS0[i], REDOR.filename);
            fclose(error_file);
		}
	}}

	fclose(fp);
	return 0;
}

int spin(const char *element){
    //Returns the spin of the nuclide, multiplied by 2 as an integer
    //Assumes that the preferred isotope is used.
    char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

        if(strcmp(element,"H")==0)
            return 1;
	else if(strcmp(element,"D")==0)
            return 2;
        else if(strcmp(element,"Li")==0)
            return 3;
        else if(strcmp(element,"Be")==0)
            return 3;
         else if(strcmp(element,"B")==0)
            return 3;
         else if(strcmp(element,"C")==0)
            return 1;
         else if(strcmp(element,"N")==0)
            return 1;
         else if(strcmp(element,"O")==0)
            return 5;
         else if(strcmp(element,"F")==0)
            return 1;
         else if(strcmp(element,"Na")==0)
            return 3;
         else if(strcmp(element,"Mg")==0)
            return 5;
         else if(strcmp(element,"Al")==0)
            return 5;
         else if(strcmp(element,"Si")==0)
            return 1;
         else if(strcmp(element,"P")==0)
            return 1;
         else if(strcmp(element,"Cl")==0)
            return 3;
         else if(strcmp(element,"Sc")==0)
            return 7;
         else if(strcmp(element,"V")==0)
            return 7;
         else if(strcmp(element,"Mn")==0)
            return 5;
         else if(strcmp(element,"Co")==0)
            return 7;
         else if(strcmp(element,"Cu")==0)
            return 3;
         else if(strcmp(element,"Ga")==0)
            return 3;
         else if(strcmp(element,"As")==0)
            return 3;
         else if(strcmp(element,"Se")==0)
            return 1;
         else if(strcmp(element,"Br")==0)
            return 3;
         else if(strcmp(element,"Rb")==0)
            return 5;
         else if(strcmp(element,"Y")==0)
            return 1;
         else if(strcmp(element,"Nb")==0)
            return 9;
         else if(strcmp(element,"Mo")==0)
            return 5;
         else if(strcmp(element,"Tc")==0)
            return 9;
         else if(strcmp(element,"Rh")==0)
            return 1;
         else if(strcmp(element,"Ag")==0)
            return 1;
         else if(strcmp(element,"Cd")==0)
            return 1;
         else if(strcmp(element,"In")==0)
            return 9;
         else if(strcmp(element,"Sn")==0)
            return 1;
         else if(strcmp(element,"Sb")==0)
            return 5;
         else if(strcmp(element,"Te")==0)
            return 1;
         else if(strcmp(element,"I")==0)
            return 5;
         else if(strcmp(element,"Cs")==0)
            return 7;
         else if(strcmp(element,"La")==0)
            return 7;
         else if(strcmp(element,"W")==0)
            return 1;
         else if(strcmp(element,"Pb")==0)
            return 1;
         else if(strcmp(element,"Tl")==0)
            return 1;
         else{
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: Gyromagnetic ratio for '%s' is unknown\n", element);
            fclose(error_file);
            exit(1);
         }
}

double RDD_1A(const char *element1, const char *element2){
    //This function return the dipolar coupling constant between two nuclei that are at a 1A distance from each other
    //Assumes that the preferred isotope is used.
    //If it is 1H or 19F it is further scaled down to correspond to SREDOR
    char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    char nuclei[2][3];
    sprintf(nuclei[0],"%s",element1); sprintf(nuclei[1],"%s",element2);
    int i;
    double gamma[2], SREDOR=1.;

    for(i=0;i<2;i++){
        if(strcmp(nuclei[i],"H")==0){
			gamma[i]=346.5883;
			SREDOR = 0.5554;
		}
		else if(strcmp(nuclei[i],"D")==0)
            gamma[i]=53.203;
        else if(strcmp(nuclei[i],"Li")==0)
            gamma[i]=134.7074;
        else if(strcmp(nuclei[i],"Be")==0)
            gamma[i]=48.7054;
         else if(strcmp(nuclei[i],"B")==0)
            gamma[i]=111.2191;
         else if(strcmp(nuclei[i],"C")==0)
            gamma[i]=87.1683;
         else if(strcmp(nuclei[i],"N")==0)
            gamma[i]=35.1433;
         else if(strcmp(nuclei[i],"O")==0)
            gamma[i]=47.0036;
         else if(strcmp(nuclei[i],"F")==0){
            gamma[i]=326.2387;
			SREDOR = 0.5554;
		 }
         else if(strcmp(nuclei[i],"Na")==0)
            gamma[i]=91.7360;
         else if(strcmp(nuclei[i],"Mg")==0)
            gamma[i]=21.2324;
         else if(strcmp(nuclei[i],"Al")==0)
            gamma[i]=90.3811;
         else if(strcmp(nuclei[i],"Si")==0)
            gamma[i]=68.9103;
         else if(strcmp(nuclei[i],"P")==0)
            gamma[i]=140.4299;
         else if(strcmp(nuclei[i],"Cl")==0)
            gamma[i]=33.9978;
         else if(strcmp(nuclei[i],"Sc")==0)
            gamma[i]=84.3247;
         else if(strcmp(nuclei[i],"V")==0)
            gamma[i]=91.2781;
         else if(strcmp(nuclei[i],"Mn")==0)
            gamma[i]=86.0926;
         else if(strcmp(nuclei[i],"Co")==0)
            gamma[i]=82.0342;
         else if(strcmp(nuclei[i],"Cu")==0)
            gamma[i]=92.1368;
         else if(strcmp(nuclei[i],"Ga")==0)
            gamma[i]=105.9912;
         else if(strcmp(nuclei[i],"As")==0)
            gamma[i]=95.5456;
         else if(strcmp(nuclei[i],"Se")==0)
            gamma[i]=66.4019;
         else if(strcmp(nuclei[i],"Br")==0)
            gamma[i]=93.9245;
         else if(strcmp(nuclei[i],"Rb")==0)
            gamma[i]=33.5898;
         else if(strcmp(nuclei[i],"Y")==0)
            gamma[i]=17.0531;
         else if(strcmp(nuclei[i],"Nb")==0)
            gamma[i]=85.0839;
         else if(strcmp(nuclei[i],"Mo")==0)
            gamma[i]=22.6851;
         else if(strcmp(nuclei[i],"Tc")==0)
            gamma[i]=78.3291;
         else if(strcmp(nuclei[i],"Rh")==0)
            gamma[i]=10.9707;
         else if(strcmp(nuclei[i],"Ag")==0)
            gamma[i]=16.2185;
         else if(strcmp(nuclei[i],"Cd")==0)
            gamma[i]=77.2267;
         else if(strcmp(nuclei[i],"In")==0)
            gamma[i]=76.4012;
         else if(strcmp(nuclei[i],"Sn")==0)
            gamma[i]=129.9657;
         else if(strcmp(nuclei[i],"Sb")==0)
            gamma[i]=83.4788;
         else if(strcmp(nuclei[i],"Te")==0)
            gamma[i]=110.2622;
         else if(strcmp(nuclei[i],"I")==0)
            gamma[i]=69.8246;
         else if(strcmp(nuclei[i],"Cs")==0)
            gamma[i]=45.7750;
         else if(strcmp(nuclei[i],"La")==0)
            gamma[i]=49.3388;
         else if(strcmp(nuclei[i],"W")==0)
            gamma[i]=14.6169;
         else if(strcmp(nuclei[i],"Pb")==0)
            gamma[i]=73.5387;
         else if(strcmp(nuclei[i],"Tl")==0)
            gamma[i]=203.3;
         else{
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: Gyromagnetic ratio for '%s' is unknown\n", nuclei[i]);
            fclose(error_file);
            exit(1);
         }
    }

    return gamma[0]*gamma[1]*SREDOR;
}

double RDD(double RDD_1A, double distance){
    //This function returns the dipolar coupling constant between two nuclei.
    return RDD_1A*pow(distance,-3.0);
}

void generate_REDORs(vector< vector<double> > &REDORs){
    //This function generates a REDOR/RESPDOR library of curves assuming a 100 Hz RDD
    //in time increments of 10us up to 100ms
    //data are organized as REDORs[(time-10us)/10us][spin*2]=DS/S0
    //dimensions are fixed at REDORs[1000][9] and the data are stored in the REDOR_library.txt file
    FILE *fp, *error_file;
    int i;
    char buffer[128];

    fp=fopen("REDOR_library.txt","r");
    if(fp==NULL){
        error_file=fopen("Errors.txt","a");
        fprintf(error_file, "\nERROR: REDOR_library.txt not found\n");
        fclose(error_file);
        exit(1);
    }

    for(i=0;i<10000;i++){
        if(fgets(buffer, sizeof(buffer), fp)!=NULL){
            sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf %lf",&REDORs[i][0],&REDORs[i][1],&REDORs[i][2],&REDORs[i][3],&REDORs[i][4],&REDORs[i][5],&REDORs[i][6],&REDORs[i][7],&REDORs[i][8]);
        }
        else{
            error_file=fopen("Errors.txt","a");
            fprintf(error_file, "\nERROR: REDOR_library.txt file is missing data\n");
            fclose(error_file);
            exit(1);
        }
    }
    fclose(fp);
}

double REDOR_DSS0(double RDD, double time, double order_parameter, int spin, vector< vector<double> > &REDORs){
   // int time_index = round(RDD*time/0.1 - 1.); //compensates for different RDD than used in the library
    int time_index = round(1000*time*RDD*order_parameter - 1.);
    time_index = (time_index>9999)*9999 + (time_index<=9999)*time_index;
    time_index = (time_index>-1)*time_index;
    return REDORs[time_index][spin-1];
}

double REDOR_pred(double time, double d_mean, double d_std, vector< vector<double> > &gaussians, vector< vector<double> > &REDORs, REDOR_dataset &REDOR){
	//This function is used to calculate the REDOR dephasing given a gaussian distribution of distances with d_mean and d_std
	//DSS0_lib is organized as DSS0_lib[time][distance] with the distance indices running from 0.1 to 20 A (indices of 0 to 199)
	//and the time time running from 0.0002 s to 0.05 s in 0.0002 s increments.
	int i;
	double DSS0 = 0., distance;

	if(REDOR.detected.size()==2){
        distance=d_mean-d_std*d_std;
        DSS0= DSS0 + 0.5*REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin,REDORs);
        distance=d_mean+d_std*d_std;
        DSS0= DSS0 + 0.5*REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin,REDORs);
	}

	else{
        for(i=0;i<200;i++){
            distance=0.1+i*0.1;
            DSS0= DSS0 + d_prob(distance, d_mean, d_std, gaussians)*REDOR_DSS0(RDD(REDOR.RDD1A,distance),time,REDOR.order_parameter,REDOR.spin,REDORs);
        }
	}
	return DSS0*REDOR.scaling_factor;
}

int create_X2_table(REDOR_dataset &REDOR, const char *support_name){
	//This function is used to create the Chi^2 table for a single curve
	//X2 is ordered as X2[d_index][std_index]
	//X2 has dimensions of X2[200][101]

	printf("Creating Chi2 table for %s\n",REDOR.filename);
	int Npoints=REDOR.tmix.size();

	int i, d_index;
    for(d_index=0;d_index<200;d_index++){
        for(i=0;i<101; i++){
            REDOR.X2[d_index][i] = 0.;
        }
    }

    vector< vector<double> > gaussians;
    gaussians.resize(401, vector<double>(101,0.));
    generate_gaussians(gaussians);

    if(REDOR.type==0){
        if(REDOR.detected.size()==1){
            #pragma omp parallel for
            for(d_index=0;d_index<200;d_index++){
                double d_mean=0.1*d_index+0.1;
                double DSS0_calc, d_std=0.;
                int std_index=0, ii;

                for(ii=0;ii<Npoints;ii++){
                    DSS0_calc = DSS0_pred(REDOR.tmix[ii],d_mean,d_std,gaussians,REDOR);
                    REDOR.X2[d_index][std_index]=REDOR.X2[d_index][std_index] + pow((REDOR.DSS0[ii]-DSS0_calc),2.) / pow((REDOR.DSS0[ii]),2.);
                }
            }
        }

        else{
            #pragma omp parallel for
            for(d_index=0;d_index<200;d_index++){
                double d_mean=0.1*d_index+0.1;
                double DSS0_calc;
                int std_index, ii;

                for(std_index=0;std_index<101; std_index++){
                    double d_std=0.1*std_index;

                    for(ii=0;ii<Npoints;ii++){
                        DSS0_calc = DSS0_pred(REDOR.tmix[ii], d_mean, d_std, gaussians, REDOR);
                        REDOR.X2[d_index][std_index]=REDOR.X2[d_index][std_index] + pow((REDOR.DSS0[ii]-DSS0_calc),2.) / pow((REDOR.DSS0[ii]),2.);
                    }
                }
            }
        }
    }

    else{//intramolecular REDOR
        vector< vector<double> > REDORs;
        REDORs.resize(10000, vector<double>(9,0.));
        generate_REDORs(REDORs);

        if(REDOR.detected.size()==1){
            #pragma omp parallel for
            for(d_index=0;d_index<200;d_index++){
                double d_mean=0.1*d_index+0.1;
                double DSS0_calc, d_std=0.;
                int std_index=0, ii;

                for(ii=0;ii<Npoints;ii++){
                    DSS0_calc = REDOR_pred(REDOR.tmix[ii], d_mean, d_std, gaussians, REDORs, REDOR);
                    REDOR.X2[d_index][std_index]=REDOR.X2[d_index][std_index] + pow((REDOR.DSS0[ii]-DSS0_calc),2.) / pow((REDOR.DSS0[ii]),2.);
                }
            }
        }

        else{
            #pragma omp parallel for
            for(d_index=0;d_index<200;d_index++){
                double d_mean=0.1*d_index+0.1;
                double DSS0_calc;
                int std_index, ii;

                for(std_index=0;std_index<101; std_index++){
                    double d_std=0.1*std_index;

                    for(ii=0;ii<Npoints;ii++){
                        DSS0_calc = REDOR_pred(REDOR.tmix[ii], d_mean, d_std, gaussians, REDORs,REDOR);
                        REDOR.X2[d_index][std_index]=REDOR.X2[d_index][std_index] + pow((REDOR.DSS0[ii]-DSS0_calc),2.) / pow((REDOR.DSS0[ii]),2.);
                    }
                }
            }
        }
    }
	return 0;
}

void write_fits(char *base_filename, const char *support_name, vector< REDOR_dataset > &REDOR){
    //This function writes out the REDOR curves form the best-fit structure as well as the range of dephasing for each curve
    //Data is stored in a CSV file
    //In the arrays listing the best-fin, minimum, and maximum distances and STDEV the indices have the following meaning:
    //0=Best Fit, 1=Smallest distance and 2=Largest distance from the surface
    FILE *error_file;
    char fits_filename[128];
    int i, j;
    double DSS0;
    FILE *out;
    int N_curves=REDOR.size();

    sprintf(fits_filename, "%s_REDOR_fits.csv", base_filename);
    remove(fits_filename);
    out=fopen(fits_filename,"w");

    if(out==NULL){
        error_file=fopen("Errors.txt","a");
        fprintf(error_file, "\nERROR: Could not write REDOR file '%s' for best fit structure\n", fits_filename);
        fclose(error_file);
    }

    vector< vector<double> > gaussians;
    gaussians.resize(401, vector<double>(101,0.));
    generate_gaussians(gaussians);

    vector< vector<double> > REDORs;
    REDORs.resize(10000, vector<double>(9,0.));

    generate_REDORs(REDORs);

    //writing the header to the CSV file
    fprintf(out,",");
    for(i=0;i<REDOR.size();i++){
        fprintf(out,"%s,,,,,",REDOR[i].filename);
    }
    fprintf(out,"\n");
    for(i=0;i<REDOR.size();i++){
        fprintf(out,"time (s),best-fit,max,min,,");
    }
    fprintf(out,"\n");

    //writing out the dephasing values to the CSV file
    for(i=1;i<=250;i++){
        double time=i*0.0002;

        for(j=0;j<REDOR.size();j++){
            fprintf(out,"%lf,",time/REDOR[j].order_parameter);
            if(REDOR[j].type==0){//surface curves
                DSS0 = DSS0_pred(time/REDOR[j].order_parameter,REDOR[j].d_index[0]/10., REDOR[j].std_index[0]/10., gaussians, REDOR[j]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_pred(time/REDOR[j].order_parameter, REDOR[j].d_index[1]/10., REDOR[j].std_index[1]/10., gaussians, REDOR[j]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_pred(time/REDOR[j].order_parameter, REDOR[j].d_index[2]/10., REDOR[j].std_index[2]/10., gaussians, REDOR[j]);
                fprintf(out,"%lf,,",DSS0);
            }
            else{
                DSS0 = REDOR_pred(time/REDOR[j].order_parameter, REDOR[j].d_index[0]/10., REDOR[j].std_index[0]/10., gaussians, REDORs, REDOR[j]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = REDOR_pred(time/REDOR[j].order_parameter, REDOR[j].d_index[1]/10., REDOR[j].std_index[1]/10., gaussians, REDORs, REDOR[j]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = REDOR_pred(time/REDOR[j].order_parameter, REDOR[j].d_index[2]/10., REDOR[j].std_index[2]/10., gaussians, REDORs, REDOR[j]);
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }

    fclose(out);
}
