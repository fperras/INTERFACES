#include "statistics.hpp"
#include "constraints.hpp"
#include "REDOR_library.hpp"

int get_Nspins(char *buffer){
    //function returns the number of nuclei in a given recoupled_spins or detected_spins line
    int word_count=0;
    int len=strlen(buffer);
    int i=1;

    while(buffer[i] != '\0'){
        if(((buffer[i]==' ')||(buffer[i]=='\n'))&&(buffer[i-1]!=' ')){
            word_count++;
            }
        i++;
    }
    return word_count-1;
}

double get_effective_distance(vector<vector<int>> &REDOR_det_index,vector<vector<int>> &REDOR_rec_index, double (*xyz)[3], int curve_index, int det_index){
    int i;
    double Deff=0., dist;

    for(i=0; i<REDOR_rec_index[curve_index].size(); i++){
        dist = distance_calc(xyz[REDOR_det_index[curve_index][det_index]],xyz[REDOR_rec_index[curve_index][i]]);
        Deff=Deff+1./pow(dist,6.);
    }
    return pow(1./Deff,1./6.);
}

double get_average_distance(vector<vector<int>> &REDOR_det_index,vector<vector<int>> &REDOR_rec_index, double (*xyz)[3], int curve_index, int curve_type){
    //This function returns the average vertical distance for the atoms that contribute to a given curve (curve_index)

    double d_sum=0.;
    int j;

    if(curve_type==0){

        for(j=0; j<REDOR_det_index[curve_index].size(); j++){
            d_sum = d_sum + xyz[REDOR_det_index[curve_index][j]][2];
        }
        return d_sum/REDOR_det_index[curve_index].size();
    }

    //will return the distance corresponding to the effective recoupled dipolar coupling in multispin case
    for(j=0; j<REDOR_det_index[curve_index].size(); j++){
        d_sum = d_sum + get_effective_distance(REDOR_det_index,REDOR_rec_index,xyz,curve_index,j);
    }
    return d_sum/REDOR_det_index[curve_index].size();
}

double get_STDEV(vector<vector<int>> &REDOR_det_index,vector<vector<int>> &REDOR_rec_index, double (*xyz)[3], int curve_index, int curve_type){
    //This function returns the standard deviation of the vertical distance for the atoms that contribute to a given curve (curve_index)
    double d_STDEV_numerator=0.;
    int j;
    double d_avg = get_average_distance(REDOR_det_index,REDOR_rec_index, xyz, curve_index,curve_type);

    if(curve_type==0){
        for(j=0; j<REDOR_det_index[curve_index].size(); j++){
            d_STDEV_numerator = d_STDEV_numerator + pow(xyz[REDOR_det_index[curve_index][j]][2] - d_avg,2.);
        }
        return sqrt(d_STDEV_numerator/REDOR_det_index[curve_index].size());
    }

    for(j=0; j<REDOR_det_index[curve_index].size(); j++){
        d_STDEV_numerator = d_STDEV_numerator + pow(get_effective_distance(REDOR_det_index,REDOR_rec_index,xyz,curve_index,j) - d_avg,2.);
    }
    return sqrt(d_STDEV_numerator/REDOR_det_index[curve_index].size());

}

int get_distance_index(vector<vector<int>> &REDOR_det_index, vector<vector<int>> &REDOR_rec_index, double (*xyz)[3], int curve_index, int curve_type){
    //This function returns the average vertical distance for the atoms that contribute to a given curve (curve_index)
    //result is given as an index in the chi squared table (distance*10)
    int i,j;

        double distance = get_average_distance(REDOR_det_index,REDOR_rec_index, xyz, curve_index,curve_type);
        int d_index=round(distance*10.+0.1001);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        return d_index;
}

int get_STDEV_index(vector<vector<int>> &REDOR_det_index,vector<vector<int>> &REDOR_rec_index, double (*xyz)[3], int curve_index, int curve_type){
    //This function returns the standard deviation of the vertical distance for the atoms that contribute to a given curve (curve_index)
    //result is given as an index in the chi squared table (STDEV*10)
    double d_std = get_STDEV(REDOR_det_index,REDOR_rec_index, xyz, curve_index,curve_type);
    int std_index = round(d_std*10+0.001);
	std_index = (std_index>100)*100 + (std_index<=100)*std_index;
	return std_index;
}

int load_simulations(const char *support_name, const char *element, vector< vector<double> > &DSS0_lib){
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

	sprintf(filename,"%s_%s.txt",support_name,element);
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
			fscanf(fp,"%lf",&DSS0_lib[i][j]);
		}
	}
	fclose(fp);
	return 0;
}

double DSS0_pred(double time, double d_mean, double d_std, vector< vector<double> > &gaussians, vector< vector<double> > &DSS0_lib, double scaling_factor, int Nspins){
	//This function is used to calculate the REDOR dephasing given a gaussian distribution of distances with d_mean and d_std
	//DSS0_lib is organized as DSS0_lib[time][distance] with the distance indices running from 0.1 to 20 A (indices of 0 to 199)
	//and the time time running from 0.0002 s to 0.05 s in 0.0002 s increments.

	int i, time_index, d_index;
	double DSS0 = 0., distance;

	time_index = round(time/0.0002 - 1.);
	time_index = (time_index>249)*249 + (time_index<=249)*time_index;
	time_index = (time_index>-1)*time_index;

	if(Nspins==2){
        distance=d_mean-d_std*d_std;
        d_index=round(distance*10.+0.1001);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0= DSS0 + 0.5*DSS0_lib[time_index][d_index];

        distance=d_mean+d_std*d_std;
        d_index=round(distance*10.+0.1001);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0= DSS0 + 0.5*DSS0_lib[time_index][d_index];
	}

	else{
        for(i=0;i<200;i++){
            distance=0.1+i*0.1;
            DSS0= DSS0 + d_prob(distance, d_mean, d_std, gaussians)*DSS0_lib[time_index][i];
        }
	}
	return DSS0*scaling_factor;
}

int find_Npoints(const char *filename){
    //Returns the number of points in a REDOR curve file
	char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

	FILE *fp;
	int Npoints=1;
	char ch;

	fp=fopen(filename,"r");

	if(fp == NULL){
        error_file=fopen(error_filename,"a");
		fprintf(error_file, "\nERROR: could not find Npoints in file %s\n", filename);
        fclose(error_file);
        exit(1);
	}

	while(!feof(fp)){
		ch = fgetc(fp);
		if(ch == '\n')
		Npoints++;
	}

	fclose(fp);
	return Npoints;
}

int load_exp_curve(const char *filename, vector<double> &DSS0, vector<double> &tmix){
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

	Npoints=find_Npoints(filename);
	fp=fopen(filename,"r");

	if(fp == NULL){
		error_file=fopen(error_filename,"a");
		fprintf(error_file, "\nERROR: could not find exp file '%s'\n", filename);
		fclose(error_file);
        exit(1);
	}

	for(i=0;i<Npoints;i++){
		fgets(buffer,64,fp);
		sscanf(buffer,"%lf %lf",&tmix[i],&DSS0[i]);
		if(DSS0[i]<=0.){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nERROR: illegal REDOR intensity of %lf found in %s\n",DSS0[i], filename);
            fclose(error_file);
            exit(1);
		}
		else if(DSS0[i]<0.05){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nWARNING: Very low dephasing value of %lf in %s can lead to over fitting\n You may obtain a better fit by eliminating this data point.\n",DSS0[i], filename);
            fclose(error_file);
		}
	}

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
    double gamma[2];

    for(i=0;i<2;i++){
        if(strcmp(nuclei[i],"H")==0)
            gamma[i]=(i)*192.4814-(i-1)*346.5883;
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
         else if(strcmp(nuclei[i],"F")==0)
            gamma[i]=(i)*181.18-(i-1)*326.2387;
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

    return gamma[0]*gamma[1];
}

double RDD(double RDD_1A, double distance){
    //This function returns the dipolar coupling constant between two nuclei.
    char error_filename[128];
	sprintf(error_filename,"Errors.txt");
    FILE *error_file;

    return RDD_1A*pow(distance,-3.0);
}

double REDOR_DSS0(double RDD,double time,int spin, vector< vector<double> > &REDORs){
   // int time_index = round(RDD*time/0.1 - 1.); //compensates for different RDD than used in the library
    int time_index = round(1000*time*RDD - 1.);
    time_index = (time_index>9999)*9999 + (time_index<=9999)*time_index;
    time_index = (time_index>-1)*time_index;
    return REDORs[time_index][spin-1];
}

double REDOR_pred(double time, double d_mean, double d_std, double RDD1A, vector< vector<double> > &gaussians, vector< vector<double> > &REDORs, double scaling_factor, int Nspins, int spin){
	//This function is used to calculate the REDOR dephasing given a gaussian distribution of distances with d_mean and d_std
	//DSS0_lib is organized as DSS0_lib[time][distance] with the distance indices running from 0.1 to 20 A (indices of 0 to 199)
	//and the time time running from 0.0002 s to 0.05 s in 0.0002 s increments.
	int i, time_index, d_index;
	double DSS0 = 0., distance;

	if(Nspins==2){
        distance=d_mean-d_std*d_std;
        DSS0= DSS0 + 0.5*REDOR_DSS0(RDD(RDD1A,distance),time,spin,REDORs);
        distance=d_mean+d_std*d_std;
        DSS0= DSS0 + 0.5*REDOR_DSS0(RDD(RDD1A,distance),time,spin,REDORs);
	}

	else{
        for(i=0;i<200;i++){
            distance=0.1+i*0.1;
            DSS0= DSS0 + d_prob(distance, d_mean, d_std, gaussians)*REDOR_DSS0(RDD(RDD1A,distance),time,spin,REDORs);
        }
	}
	return DSS0*scaling_factor;
}



int create_X2_table(const char *filename, const char *support_name, const char *element, const char *rec_element, vector< vector<double> > &X2, double scaling_factor, int Nspins, int curve_type){
	//This function is used to create the Chi^2 table for a single curve
	//X2 is ordered as X2[d_index][std_index]
	//X2 has dimensions of X2[200][101]

	printf("creating Chi2 table for %s\n",filename);
	int Npoints=find_Npoints(filename);

	vector<double> DSS0(Npoints,0.);
	vector<double> tmix(Npoints,0.);
	load_exp_curve(filename, DSS0, tmix);



	int i, d_index;
    for(d_index=0;d_index<200;d_index++){
        for(i=0;i<101; i++){
            X2[d_index][i] = 0.;
        }
    }

    vector< vector<double> > gaussians;
    gaussians.resize(401, vector<double>(101,0.));
    generate_gaussians(gaussians);

    if(curve_type==0){
        vector< vector<double> > DSS0_lib;
        DSS0_lib.resize(250, vector<double>(200,0.));
        load_simulations(support_name, element, DSS0_lib);

        if(Nspins==1){
            #pragma omp parallel for
            for(d_index=0;d_index<200;d_index++){
                double d_mean=0.1*d_index+0.1;
                double DSS0_calc, d_std=0.;
                int std_index=0, ii;

                for(ii=0;ii<Npoints;ii++){
                    DSS0_calc = DSS0_pred(tmix[ii], d_mean, d_std, gaussians, DSS0_lib,scaling_factor, Nspins);
                    X2[d_index][std_index]=X2[d_index][std_index] + pow((DSS0[ii]-DSS0_calc),2.) / pow((DSS0[ii]),2.);
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
                        DSS0_calc = DSS0_pred(tmix[ii], d_mean, d_std, gaussians, DSS0_lib,scaling_factor,Nspins);
                        X2[d_index][std_index]=X2[d_index][std_index] + pow((DSS0[ii]-DSS0_calc),2.) / pow((DSS0[ii]),2.);
                    }
                }
            }
        }
    }

    else{//intramolecular REDOR
        vector< vector<double> > REDORs;
        REDORs.resize(10000, vector<double>(9,0.));
        generate_REDORs(REDORs);

        double RDD1A=RDD_1A(element,rec_element);
        int S=spin(rec_element);

        if(Nspins==1){
            #pragma omp parallel for
            for(d_index=0;d_index<200;d_index++){
                double d_mean=0.1*d_index+0.1;
                double DSS0_calc, d_std=0.;
                int std_index=0, ii;

                for(ii=0;ii<Npoints;ii++){
                    DSS0_calc = REDOR_pred(tmix[ii], d_mean, d_std, RDD1A, gaussians, REDORs,scaling_factor, Nspins, S);
                    X2[d_index][std_index]=X2[d_index][std_index] + pow((DSS0[ii]-DSS0_calc),2.) / pow((DSS0[ii]),2.);
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
                        DSS0_calc = REDOR_pred(tmix[ii], d_mean, d_std, RDD1A, gaussians, REDORs,scaling_factor, Nspins, S);
                        X2[d_index][std_index]=X2[d_index][std_index] + pow((DSS0[ii]-DSS0_calc),2.) / pow((DSS0[ii]),2.);
                    }
                }
            }
        }
    }
	return 0;
}

int write_fits(int (*distances)[3], int (*stdevs)[3], char *base_filename, int N_curves, const char *support_name, vector<vector<int>> &REDOR_det_index, vector<vector<int>> &REDOR_rec_index, char (*element)[3], char (*curve_filename)[120], double *scaling_factor, int *Nspins, int *curve_type){
    //This function writes out the REDOR curves form the best-fit structure as well as the range of dephasing for each curve
    //Data is stored in a CSV file
    //In the arrays listing the best-fin, minimum, and maximum distances and STDEV the indices have the following meaning:
    //0=Best Fit, 1=Smallest distance and 2=Largest distance from the surface
    FILE *error_file;
    char fits_filename[128];
    int i, j, k;
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

    vector< vector< vector<double> > > DSS0_lib;
    DSS0_lib.resize(N_curves, vector<vector<double> > (250, vector<double>(200,0.)));

    for(j=0;j<N_curves;j++){
            if(curve_type[j]==0)
                load_simulations(support_name, element[REDOR_det_index[j][0]], DSS0_lib[j]);
    }

    vector< vector<double> > gaussians;
    gaussians.resize(401, vector<double>(101,0.));
    generate_gaussians(gaussians);

    vector< vector<double> > REDORs;
    REDORs.resize(10000, vector<double>(9,0.));
    generate_REDORs(REDORs);

    //writing the header to the CSV file
    fprintf(out,",");
    for(i=0;i<N_curves;i++){
        fprintf(out,",%s,,,",curve_filename[i]);
    }
    fprintf(out,"\ntime (s),,");
    for(i=0;i<N_curves;i++){
        fprintf(out,"best-fit,min,max,,");
    }
    fprintf(out,"\n");

    //writing out the dephasing values to the CSV file
    for(i=1;i<=250;i++){
        double time=i*0.0002;
        fprintf(out,"%lf,,",time);

        for(j=0;j<N_curves;j++){
            if(curve_type[j]==0){//surface curves
                DSS0 = DSS0_pred(time,distances[j][0]/10., stdevs[j][0]/10., gaussians, DSS0_lib[j],scaling_factor[j], Nspins[j]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_pred(time, distances[j][1]/10., stdevs[j][1]/10., gaussians, DSS0_lib[j],scaling_factor[j], Nspins[j]);
                fprintf(out,"%lf,",DSS0);
                DSS0 = DSS0_pred(time, distances[j][2]/10., stdevs[j][2]/10., gaussians, DSS0_lib[j],scaling_factor[j], Nspins[j]);
                fprintf(out,"%lf,,",DSS0);
            }
            else{
                double RDD1A=RDD_1A(element[REDOR_det_index[j][0]],element[REDOR_rec_index[j][0]]);
                int S=spin(element[REDOR_rec_index[j][0]]);
                DSS0 = REDOR_pred(time, distances[j][0]/10., stdevs[j][0]/10.,RDD1A, gaussians, REDORs,scaling_factor[j], Nspins[j],S);
                fprintf(out,"%lf,",DSS0);
                DSS0 = REDOR_pred(time, distances[j][1]/10., stdevs[j][1]/10., RDD1A, gaussians, REDORs,scaling_factor[j], Nspins[j],S);
                fprintf(out,"%lf,",DSS0);
                DSS0 = REDOR_pred(time, distances[j][2]/10., stdevs[j][2]/10., RDD1A, gaussians, REDORs,scaling_factor[j], Nspins[j],S);
                fprintf(out,"%lf,,",DSS0);
            }
        }
        fprintf(out,"\n");
    }

    fclose(out);

    return 0;
}
