#include "stdio.h"
#include "math.h"
#include "stdlib.h"
#include <time.h>
#include <mpi.h>
float coswDt(float alpha, float beta, float gamma, float time, float RDD, float alpha2, float beta2, int spin);
void compile_curves();
static double pi = 3.1415926535897932384626433;

int main(){
//This program is used to generate a surface-to-atom REDOR or RESPDOR curve library between a given nuclide and a given support
//The structure is given in a file named "structure.txt" which lists the coordinates of the spins form the support in a 3-column: x y z format
//The surface place needs to appear at z=0 and the structure centered on x=y=0. The program calculates the RE(SP)DOR curves for a single spin
//situated at the coordinates (0,0,d), where d is the distance from the atom to the support interface.
//The program will finish by creating a file called RESPDOR_library.txt in the format expected by INTERFACES. This file should be renamed as: support_element.txt

// Initialize the MPI environment
    MPI_Init(NULL, NULL);
    int world_size, world_rank, name_len;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    MPI_Get_processor_name(processor_name, &name_len);

//These are the variables that need to be set by the user
    float NA=1.0; //natural abundance, or enrichment level
    int N_atoms = 104; //number of surface atoms in the structure.txt file
    int spin = 5; //spin quantum number of the surface atoms, multiplied by a factor of 2
    float RDD_1A = 4097.223; //dipolar coupling at a distance of 1A. (reduce further by a factor of 1.8 if it is a S-RE(SP)DOR experiment)

//variable declarations. Do not change time_steps, max_time, or init_time if you wish to maintain compatibility with INTERFACES
    int i,j,k,l,t, betasteps = 71, alphasteps = 51, time_steps = 250;
    float start_distance = (1 + world_rank)*0.1;
    float distance_increment = 0.1*world_size;
    float max_time = 0.050, init_time= 0.0002, distance, DSnorm[time_steps+1], DS, RDD,  betaD, alpha, beta, gamma;
    float x[N_atoms], y[N_atoms], z[N_atoms], xy[N_atoms], alphaD[N_atoms];
    FILE *fp, *struct_fp;
    char filename[30], buffer[256];

//Here the program reads the structure file, which is just a list of x y z coordinates in a 3-column text file
//It then converts these into a z distance, xy distance, and alpha angle.
        struct_fp=fopen("structure.txt","r");
        for(i=0;i<N_atoms;i++){
            fgets(buffer,48,struct_fp);
            sscanf(buffer,"%f  %f  %f\n",&x[i],&y[i],&z[i]);
            xy[i] = sqrt(x[i]*x[i]+y[i]*y[i]);
            z[i] = fabs(z[i]);
            if(x[i]>0.)
                alphaD[i]=atan(y[i]/x[i]);
            else if(x[i]<0.)
                alphaD[i]=atan(y[i]/x[i])+pi;
            else
                alphaD[i]=0.0;
        }
        fclose(struct_fp);

//Here the program starts to loop over the distances from the support
    for(distance=start_distance; distance<=20.0; distance = distance + distance_increment){
        printf("calculating distance of %.2f by processor %d\n", distance, world_rank);
        sprintf(filename, "RESPDOR_%.2f.txt", distance);

        for(i=0; i<time_steps+1; i++){
            DSnorm[i] = 0;
        }

//calculation of RESPDOR curve
        for(t=0; t<time_steps; t++){  //start of time loop
            for(i=0; i<betasteps; i++){  //start of beta loop (cosinusoidal to avoid powder average)
               beta = acos(1.-i/(0.5*betasteps));
               for(alpha=0.; alpha<2.*pi; alpha = alpha+ 2.*pi/alphasteps){  //start of alpha loop
                    for(gamma=0.; gamma<2.*pi; gamma = gamma+ 2.*pi/alphasteps){  //start of gamma loop
                    DS=1.;
                    for(j=0;j<N_atoms;j++){
                        RDD=RDD_1A/pow(xy[j]*xy[j]+pow(distance+z[j],2.),1.5);
                        betaD=atan(xy[j]/(distance+z[j]));
                        DS=DS*(1.-NA*(1.-    (-coswDt(alpha, beta, gamma, (max_time-init_time)*t/(time_steps-1)+init_time, RDD, alphaD[j], betaD, spin))));
                    }//end of pairs loop

                   //The dephasings are then multiplied for that particular orientation
                   DSnorm[t]=DSnorm[t] + (1. - DS)/(alphasteps*betasteps*alphasteps);
} //end of gamma loop
} //end of alpha loop
} //end of beta loop
} //end of time loop

//creating an output file
    fp=fopen(filename, "w");
    for(i=0; i<time_steps; i++){
        fprintf(fp, "%f  %f\n",(max_time-init_time)*i/(time_steps-1)+init_time, DSnorm[i]);
    }
    fclose(fp);
    }//end of distance loop

 //finalizing
    MPI_Barrier(MPI_COMM_WORLD);
    compile_curves();
    MPI_Finalize();
}

float coswDt(float alpha, float beta, float gamma, float time, float RDD, float alpha2, float beta2, int spin) {
    //This is a function to calculate the cosine of the dipolar frequency
    //equations taken from JMR 127, 147-154 (1997) and PCCP 12, 9395-9405 (2010)

       float x = sin(beta2)*cos(alpha2);
       float y = sin(beta2)*sin(alpha2);
       float z = cos(beta2);

       float niz = sin(beta)*cos(alpha)*x + sin(beta)*sin(alpha)*y + cos(beta)*z;
       float niy = (-cos(gamma)*sin(alpha) - cos(beta)*cos(alpha)*sin(gamma))*x + (cos(gamma)*cos(alpha) - cos(beta)*sin(alpha)*sin(gamma))*y + sin(gamma)*sin(beta)*z;
       float s2bca = 2*niz*niy;

       switch(spin){
			case 1 :
			return -cos(2.*sqrt(2.)*time*1*(RDD*s2bca));
			break;

			case 2 :
			return -1./3. -4./9.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca))-2./9.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca));
			break;

			case 3 :
			return -1./4. -3./8.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca))-1./4.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca))-1./8.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca));
			break;

			case 4 :
			return -1./5. -8./25.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca))-6./25.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca))-4./25.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca))-2./25.*cos(2.*sqrt(2.)*time*4*(RDD*s2bca));
			break;

			case 5 :
			return -1./6. - 5./18.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca)) - 2./9.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca)) - 1./6.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca)) - 1./9.*cos(2.*sqrt(2.)*time*4*(RDD*s2bca)) - 1./18.*cos(2.*sqrt(2.)*time*5*(RDD*s2bca));
			break;

			case 6 :
			return -1./7. - 12./49.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca)) - 10./49.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca)) - 8./49.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca)) - 6./49.*cos(2.*sqrt(2.)*time*4*(RDD*s2bca)) - 4./49.*cos(2.*sqrt(2.)*time*5*(RDD*s2bca)) - 2./49.*cos(2.*sqrt(2.)*time*6*(RDD*s2bca));
			break;

			case 7 :
			return -1./8. - 14./64.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca)) - 12./64.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca)) - 10./64.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca)) - 8./64.*cos(2.*sqrt(2.)*time*4*(RDD*s2bca)) - 6./64.*cos(2.*sqrt(2.)*time*5*(RDD*s2bca)) - 4./64.*cos(2.*sqrt(2.)*time*6*(RDD*s2bca)) - 2./64.*cos(2.*sqrt(2.)*time*7*(RDD*s2bca));
			break;

			case 8 :
			return -1./9. - 16./81.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca)) - 14./81.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca)) - 12./81.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca)) - 10./81.*cos(2.*sqrt(2.)*time*4*(RDD*s2bca)) - 8./81.*cos(2.*sqrt(2.)*time*5*(RDD*s2bca)) - 6./81.*cos(2.*sqrt(2.)*time*6*(RDD*s2bca)) - 4./81.*cos(2.*sqrt(2.)*time*7*(RDD*s2bca)) - 2./81.*cos(2.*sqrt(2.)*time*8*(RDD*s2bca));
			break;

			case 9 :
			return -1./10. - 9./50.*cos(2.*sqrt(2.)*time*1*(RDD*s2bca)) - 8./50.*cos(2.*sqrt(2.)*time*2*(RDD*s2bca)) - 7./50.*cos(2.*sqrt(2.)*time*3*(RDD*s2bca)) - 6./50.*cos(2.*sqrt(2.)*time*4*(RDD*s2bca)) - 5./50.*cos(2.*sqrt(2.)*time*5*(RDD*s2bca)) - 4./50.*cos(2.*sqrt(2.)*time*6*(RDD*s2bca)) - 3./50.*cos(2.*sqrt(2.)*time*7*(RDD*s2bca)) - 2./50.*cos(2.*sqrt(2.)*time*8*(RDD*s2bca)) - 1./50.*cos(2.*sqrt(2.)*time*9*(RDD*s2bca));
			break;
	   }
    return 0.0;
}

void compile_curves(){
    float time, distance, dephasing;
    int i, j;
    FILE *fp[200], *out;
    char buffer[64], filename[64];

    out=fopen("RESPDOR_library.txt","w");

    i=0;
    for(distance=0.1; distance<=20.01; distance=distance+0.1){
        sprintf(filename,"RESPDOR_%.2f.txt",distance);
        fp[i]=fopen(filename,"r");
        i++;
    }

    for(i=0;i<250;i++){
        fprintf(out,"%f  ", i*0.0002);
        for(j=0;j<200;j++){
            fgets(buffer,64,fp[j]);
            sscanf(buffer,"%f %f",&time,&dephasing);
            fprintf(out,"%f  ",dephasing);
        }
        fprintf(out,"\n");
    }

    i=0;
    for(distance=0.1; distance<=20.01; distance=distance+0.1){
        fclose(fp[i]);
        sprintf(filename,"RESPDOR_%.2f.txt",distance);
        remove(filename);
        i++;
    }
    fclose(out);
}
