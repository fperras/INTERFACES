#include <string.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <math.h>
#include <omp.h>
using namespace std;

double calc_chi2(vector<vector< vector<double> > > &X2, int *std_index, int *d_index, int N_curves) {
    //This function compiles the total Chi^2 of a given structure by summing the values from the individual curves
    double chi2 = 0;
    int i;

    for(i=0; i<N_curves; i++){
        chi2 = chi2+ X2[i][d_index[i]][std_index[i]];
    }
    return chi2;
}

double sign(const double x){
    //returns -1 for negative, 1 for positive, 0 for 0.
    return (0 < x) - (x < 0);
}

double inverf(const double x){
    //Function used to calculate the inverse error function.
    double r;
    double a[] = {0.886226899, -1.645349621, 0.914624893, -0.140543331};
    double b[] = {1, -2.118377725, 1.442710462, -0.329097515, 0.012229801};
    double c[] = {-1.970840454, -1.62490649, 3.429567803, 1.641345311};
    double d[] = {1, 3.543889200, 1.637067800};

    double z = sign(x) * x;

    if(z <= 0.7){
        double x2 = z * z;
        r = z * (((a[3] * x2 + a[2]) * x2 + a[1]) * x2 + a[0]);
        r /= (((b[4] * x2 + b[3]) * x2 + b[2]) * x2 + b[1])* x2 + b[0];
    } else {
     double y = sqrt( -log((1 - z)/2));
     r = (((c[3] * y + c[2]) * y + c[1]) * y + c[0]);
     r /= ((d[2] * y + d[1]) * y + d[0]);
   }

   r = r * sign(x);
   z = z * sign(x);

   r -= (erf(r) - z)/(2/sqrt(M_PI) *exp(-r * r));
   r -= (erf(r) - z)/(2/sqrt(M_PI) *exp(-r * r)); //Comment out if you want single refinement

   return r;
}

double max_Chi2(double Chi2_min, double confidence_interval){
    //function to calculate the max Chi^2 allowed for a given confidence interval
    //Note that the confidence interval is defined as a decimal, not a percentage
    return (2.*pow(inverf(confidence_interval),2.)+1.)*Chi2_min;
}

double max_Chi2_multi(double Chi2_best, double confidence_interval){
    //function to calculate the max Chi^2 allowed for a given confidence interval
    //Note that the confidence interval is defined as a decimal, not a percentage
    //This function differs in that it defines the Chi^2 level required to exclude
    //a given Chi^2 value.
    return Chi2_best/(2.*pow(inverf(confidence_interval),2.)+1.);
}

double return_CI(double Chi2_large, double Chi2_small){
    //This function uses the Chi^2 minimum and a given Chi^2 value to return the
    //corresponding confidence interval.
    return erf(sqrt((Chi2_large/Chi2_small-1.)/2.));
}

void generate_gaussians(vector< vector<double> > &gaussians){
    //The function pre-computes a series of gaussian distributions of distances
    //the vector is formatted with the first index for offset distance, in 0.1A increments from -20 to +20
    //the second index is the standard deviation from 0 to 10 A in 0.1A;
    //basically, double gaussians[401][101];
    int i,j;
    float d,dev;
    for(j=0;j<101;j++){
        for(i=0;i<401;i++){
            d=0.1*i;
            dev=0.1*j;
            gaussians[i][j] = 0.5*erf(0.707107*(20.-(d-0.05))/dev)-0.5*erf(0.707107*(20.-(d+0.05))/dev);
    }}
}

double d_prob(double distance, double d_mean, double d_std, vector< vector<double> > &gaussians){
	//This function uses the gaussians to calculate the probability of a given distance
	//note 0 is false 1 is true
	//signbit 0 for positive 1 for negative

	int d_index, std_index;

	std_index = round(d_std*10+0.001);
	std_index = (std_index>100)*100 + (std_index<=100)*std_index;

	distance = distance-d_mean;
	d_index= 200 + round(distance*10+0.001);
	d_index = (d_index>400)*400 + ((d_index<=400) && (d_index>0))*d_index;

	//printf("gaussians[%d][%d] = %lf\n",d_index,std_index,gaussians[d_index][std_index]);

	return gaussians[d_index][std_index];
}
