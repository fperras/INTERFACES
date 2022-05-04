#include "statistics.hpp"

double get_average_distance(vector<vector<int>> &REDOR_atom_index, double (*xyz)[3], int N_curves, int curve_index){
    //This function returns the average vertical distance for the atoms that contribute to a given curve (curve_index)
    double NMR_z_coord_sum, NMR_z_avg;
    int j;

    NMR_z_coord_sum=0.;

    for(j=0; j<REDOR_atom_index[curve_index].size(); j++){
        NMR_z_coord_sum = NMR_z_coord_sum + xyz[REDOR_atom_index[curve_index][j]][2];
    }
    NMR_z_avg = NMR_z_coord_sum/REDOR_atom_index[curve_index].size();
    return NMR_z_avg;
}

double get_STDEV(vector<vector<int>> &REDOR_atom_index, double (*xyz)[3], int N_curves, int curve_index){
    //This function returns the standard deviation of the vertical distance for the atoms that contribute to a given curve (curve_index)
    double NMR_z_avg;
    double z_STDEV_numerator=0.;
    int j, k=0;

    NMR_z_avg = get_average_distance(REDOR_atom_index, xyz, N_curves, curve_index);

    for(j=0; j<REDOR_atom_index[curve_index].size(); j++){
        z_STDEV_numerator = z_STDEV_numerator + pow(xyz[REDOR_atom_index[curve_index][j]][2] - NMR_z_avg,2);
    }
    return sqrt(z_STDEV_numerator/REDOR_atom_index[curve_index].size());
}

int get_distance_index(vector<vector<int>> &REDOR_atom_index, double (*xyz)[3], int N_curves, int curve_index){
    //This function returns the average vertical distance for the atoms that contribute to a given curve (curve_index)
    //result is given as an index in the chi squared table (distance*10)
    double distance = get_average_distance(REDOR_atom_index, xyz, N_curves, curve_index);
    int d_index=round(distance*10.+0.1001);
    d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
    return d_index;
}

int get_STDEV_index(vector<vector<int>> &REDOR_atom_index, double (*xyz)[3], int N_curves, int curve_index){
    //This function returns the standard deviation of the vertical distance for the atoms that contribute to a given curve (curve_index)
    //result is given as an index in the chi squared table (STDEV*10)
    double d_std = get_STDEV(REDOR_atom_index, xyz, N_curves, curve_index);
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

	if(Nspins==2){
        distance=d_mean-d_std*d_std;
        d_index=round(distance*10.+0.1001);
        DSS0= DSS0 + 0.5*DSS0_lib[time_index][d_index];

        distance=d_mean+d_std*d_std;
        d_index=round(distance*10.+0.1001);
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
		else if(DSS0[i]<0.1){
            error_file=fopen(error_filename,"a");
            fprintf(error_file, "\nWARNING: Very low dephasing value of %lf in %s can lead to over fitting\n You may obtain a better fit by eliminating this data point.\n",DSS0[i], filename);
            fclose(error_file);
		}
	}

	fclose(fp);
	return 0;
}

int create_X2_table(const char *filename, const char *support_name, const char *element, vector< vector<double> > &X2, double scaling_factor, int Nspins){
	//This function is used to create the Chi^2 table for a single curve
	//X2 is ordered as X2[d_index][std_index]
	//X2 has dimensions of X2[200][101]

	int Npoints=find_Npoints(filename);

	vector<double> DSS0(Npoints,0.);
	vector<double> tmix(Npoints,0.);
	load_exp_curve(filename, DSS0, tmix);


	vector< vector<double> > DSS0_lib;
	DSS0_lib.resize(250, vector<double>(200,0.));
	//printf("finding support/n");
	load_simulations(support_name, element, DSS0_lib);

	vector< vector<double> > gaussians;
	gaussians.resize(401, vector<double>(101,0.));
	generate_gaussians(gaussians);


	printf("creating Chi2 table for %s\n",filename);
	int i, d_index;
	for(d_index=0;d_index<200;d_index++){
		for(i=0;i<101; i++){
			X2[d_index][i] = 0.;
		}
	}

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

	return 0;
}

int write_fits(int (*distances)[3], int (*stdevs)[3], char *base_filename, int N_curves, const char *support_name, vector<vector<int>> &REDOR_atom_index, char (*element)[3], char (*curve_filename)[120], double scaling_factor, int *Nspins){
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
            load_simulations(support_name, element[REDOR_atom_index[j][0]], DSS0_lib[j]);
    }

    vector< vector<double> > gaussians;
    gaussians.resize(401, vector<double>(101,0.));
    generate_gaussians(gaussians);

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
            DSS0 = DSS0_pred(time, distances[j][0]/10., stdevs[j][0]/10., gaussians, DSS0_lib[j],scaling_factor, Nspins[j]);
            fprintf(out,"%lf,",DSS0);
            DSS0 = DSS0_pred(time, distances[j][1]/10., stdevs[j][2]/10., gaussians, DSS0_lib[j],scaling_factor, Nspins[j]);
            fprintf(out,"%lf,",DSS0);
            DSS0 = DSS0_pred(time, distances[j][2]/10., stdevs[j][1]/10., gaussians, DSS0_lib[j],scaling_factor, Nspins[j]);
            fprintf(out,"%lf,,",DSS0);
        }
        fprintf(out,"\n");
    }

    fclose(out);

    return 0;
}
