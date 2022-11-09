//#include "REDOR_data.hpp"

double DSS0_full_old(double time, double scaling_factor, double order_parameter, vector<int> &REDOR_det_index, vector<vector<double> > xyz, int curve_index, vector< vector<double> > DSS0_lib){
    //This function returns the dephasing form a group of detected spins as an average
    //for surface-to-atom REDOR

    int time_index = round((time/0.0002)*order_parameter - 1.);
	time_index = (time_index>249)*249 + (time_index<=249)*time_index;
	time_index = (time_index>-1)*time_index;
	double DSS0=0.;
	double distance;
	int i,d_index, ndet=REDOR_det_index.size();

    for(i=0; i<ndet; i++){
        distance = xyz[REDOR_det_index[i]][2];
        d_index=round(distance*10.+0.1001-1.);
        d_index = (d_index>199)*199 + ((d_index<=199) && (d_index>0))*d_index;
        DSS0 = DSS0 + DSS0_lib[time_index][d_index]/ndet;
    }

    return DSS0*scaling_factor;

}

double calculate_curve_Chi2_old(double scaling_factor, double order_parameter, vector<int> &REDOR_det_index, vector<vector<double> > xyz, int curve_index, int curve_type, vector< vector<double> > &REDORs, vector< vector<double> > DSS0_lib, vector<double> &DSS0,vector<double> tmix){
    //This function returns the Chi2 contribution from a single curve.


    int i, Npoints = DSS0.size();
    double X2=0, DSS0_calc;
   // double RDD1A=RDD_1A(element[REDOR_det_index[curve_index][0]],element[REDOR_rec_index[curve_index][0]]);
    //int S=spin(element[REDOR_rec_index[curve_index][0]]);


    if(curve_type==0){//surface-to-atom
        for(i=0; i<Npoints; i++){
            DSS0_calc=DSS0_full_old(tmix[i],scaling_factor,order_parameter,REDOR_det_index,xyz,curve_index,DSS0_lib);
            X2 = X2 + pow((DSS0[i]-DSS0_calc),2.) / pow((DSS0[i]),2.);
        }
        return X2;
    }

  //  for(i=0; i<Npoints; i++){
    //    DSS0_calc=REDOR_full(tmix[i],scaling_factor,order_parameter,RDD1A,S,REDOR_det_index,REDOR_rec_index,xyz,curve_index,REDORs);
      //  X2 = X2 + pow((DSS0[i]-DSS0_calc),2.) / pow((DSS0[i]),2.);
   // }
    return X2;
}


void precalculate_dephasing_old(vector<double> &base_DSS0, double scaling_factor, double order_parameter, vector<int> &REDOR_det_index, vector<vector<double> > xyz, int curve_index, int curve_type, vector< vector<double> > &REDORs, vector< vector<double> > DSS0_lib, vector<double> &DSS0,vector<double> tmix){
    //This function calculated the dephasing levels for a given set of xyz coordinates at the experimentally-specified recoupling times
    //Designed to be used with the calculate_curve_Chi2_multi() function.

    int i,j, Npoints = DSS0.size();
  //  double RDD1A=RDD_1A(element[REDOR_det_index[curve_index][0]],element[REDOR_rec_index[curve_index][0]]);
    //int S=spin(element[REDOR_rec_index[curve_index][0]]);

    if(curve_type==0){//surface-to-atom
        for(i=0; i<Npoints; i++){
            base_DSS0[i]=DSS0_full_old(tmix[i],scaling_factor,order_parameter,REDOR_det_index,xyz,curve_index,DSS0_lib);
        }
        return;
    }

//    for(i=0; i<Npoints; i++){
  //      base_DSS0[i]=REDOR_full(tmix[i],scaling_factor,order_parameter,RDD1A,S,REDOR_det_index,REDOR_rec_index,xyz,curve_index,REDORs);
    //}
    return;
}

double calculate_curve_Chi2_multi_old(double *curve_chi2, vector<double> base_DSS0, double scaling_factor, double order_parameter, vector<int> &REDOR_det_index, vector<vector<double> > xyz, int curve_index, vector< vector<double> > &REDORs, vector< vector<double> > DSS0_lib, vector<double> &DSS0,vector<double> tmix){
    //This function returns the Chi2 contributions from a single curve in a system containing multiple molecules
    //The function loops over site populations from 5% to 50 % in 5% increments and stores the curve_chi2 values
    //in an array of the same name. These values can then be compared to determine whether the site is statistically
    //significant or not. The dephasing levels determined from a prior structure determination are imported as
    //base_DSS0[Npoints].

    int i,j, Npoints = DSS0.size();
    double DSS0_calc, weight, minimum=100000000.;
    vector<double> new_DSS0;
    new_DSS0.resize(Npoints, 0.);
    precalculate_dephasing_old(new_DSS0,scaling_factor,order_parameter,REDOR_det_index,xyz,curve_index,0,REDORs,DSS0_lib,DSS0,tmix);
    int curve_type=0;
    if(curve_type==0){//surface-to-atom
        for(j=0;j<10;j++){
            weight = 0.05 + 0.05*j;
            curve_chi2[j]=0.0;
            for(i=0; i<Npoints; i++){
                DSS0_calc = weight*new_DSS0[i] + (1.-weight)*base_DSS0[i];
                curve_chi2[j] = curve_chi2[j] + pow((DSS0[i]-DSS0_calc),2.) / pow((DSS0[i]),2.);
            }
            minimum = (minimum < curve_chi2[j])*minimum + (minimum > curve_chi2[j])*curve_chi2[j];
        }
        return minimum;
    }

  /*  for(j=0;j<10;j++){
        weight = 0.05 + 0.05*j;
        curve_chi2[j]=0.0;
        for(i=0; i<Npoints; i++){
            DSS0_calc = weight*new_DSS0[i] + (1.-weight)*base_DSS0[i];
            curve_chi2[j] = curve_chi2[j] + pow((DSS0[i]-DSS0_calc),2.) / pow((DSS0[i]),2.);
        }
        minimum = (minimum < curve_chi2[j])*minimum + (minimum > curve_chi2[j])*curve_chi2[j];
    }*/
    return minimum;
}


