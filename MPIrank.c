#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <mpi.h>
// main menu
/*Function to Compute 2-norm of a vector */
double vector_norm(double [], int );
int main (int argc, char *argv[])
 {
double vnorm;
int totalnodes=11;
double normalized=0;
double d=0.85;
 int k,p;

double M[11][11]={{0.09091,0.00000,0.00000,0.50000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,1.00000,0.50000,0.33333,0.50000,0.50000,0.50000,0.50000,0.00000,0.00000},
{0.09091,1.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.33333,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.00000,0.50000,0.50000,0.50000,0.50000,1.00000,1.00000},
{0.09091,0.00000,0.00000,0.00000,0.33333,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000},
{0.09091,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000,0.00000}};
double A[11][11],R[11],V[11],X[11],K[11],A1[4][11],A2[4][11],A3[3][11],V1[4],V2[4],V3[3],X1[4],X2[4],X3[3],K1[4],K2[4],K3[3];
int check,a,i,j,cnt;
i=0;
j=0;
cnt=1;

// Initialize the MPI environment
 // MPI_Init(NULL, NULL);
int ierr;
ierr=MPI_Init(&argc, &argv);
  // Find out rank, size
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
MPI_Status status;
  // We are assuming  4 processes for this task
  if (world_size!=4) {
    fprintf(stderr, "Please World size must be equalled to 4 for this task set \n");

    MPI_Abort(MPI_COMM_WORLD, 1);
  }
// Master distribute jobs from bucket row wise Matrix to workers
if (world_rank == 0) {
    for(k=0;k<11;k++){
			for(p=0;p<11;p++){
			   if((k>3)&&(k<8)){
                    A2[k-4][p]=d*M[k][p]+(1-d)/totalnodes;
			MPI_Send( &A2[k-4][p], 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);

			    }
			     if(k>=8){
                    A3[k-8][p]=d*M[k][p]+(1-d)/totalnodes;
			MPI_Send( &A3[k-8][p], 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);

			    }
			    if(k<=3){
                   A1[k][p]=d*M[k][p]+(1-d)/totalnodes;
		      MPI_Send( &A1[k][p], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);

			    }
			}
		}
    }
    // Synchronize
    MPI_Barrier(MPI_COMM_WORLD);
//workers reveive nodes from master
if(world_rank==1){
	for (i=0;i<11;i++){
		for(j=0;j<11;j++){
			if(i<=3){

				MPI_Recv(&A1[i][j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


			}
		}
	}
}
if(world_rank==2){

	for (i=0;i<11;i++){
		for(j=0;j<11;j++){
			if((i>3)&&(i<8)){
				MPI_Recv(&A2[i-4][j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


			}
		}
	}
}
if(world_rank==3){

	for (i=0;i<11;i++){
		for(j=0;j<11;j++){
			if(i>=8){
				MPI_Recv(&A3[i-8][j], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);


			}
		}
	}
}
// Synchronize
MPI_Barrier(MPI_COMM_WORLD);
//Master
//TODO- receive individual index,initialize send & receive buff for each worker - initialize global weights,send weights[index_i]

if(world_rank==0){

	V[0]=0.8147;
	V[1]=0.9058;
	V[2]=0.1270;
	V[3]=0.9134;
	V[4]=0.6324;
	V[5]=0.0975;
	V[6]=0.2785;
	V[7]=0.5469;
	V[8]=0.9575;
	V[9]=0.9649;
	V[10]=0.1576;

	vnorm = vector_norm(V, 11);
    for(i=0;i<11;i++){

       V[i] = V[i] / vnorm;
    }

    for(i=0;i<11;i++){


        if(i<4){
          V1[i]=V[i];
          K1[i]=V[i];
        }
        if((i>=4)&&(i<8)){
          V2[i-4]=V[i];
          K2[i-4]=V[i];
        }
        if(i>=8){
            V3[i-8]=V[i];
            K3[i-8]=V[i];
        }

	//Sending vector rank
    MPI_Send(&V[i], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
    MPI_Send(&V[i], 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
    MPI_Send(&V[i], 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);

    }

}
// Synchronize
MPI_Barrier(MPI_COMM_WORLD);
//Iteration
//TODO-Worker update local graph using receiving weight
int res=0;
k=1;
while(res==0){
MPI_Barrier(MPI_COMM_WORLD);

if(world_rank==1){
	for(i=0;i<11;i++){
    	MPI_Recv(&V[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}

    for(i=0;i<11;i++){
         if(i<4){
        X1[i]=0.0;
        }
        for(j=0;j<11;j++){
            if(i<4){

                X1[i] +=A1[i][j]*V[j];
            }

        }
         if(i<4){
            V1[i]=X1[i];

		MPI_Send(&V1[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

         }


    }

}


if(world_rank==2){
for(i=0;i<11;i++){
    	MPI_Recv(&V[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}

for (i=0;i<11;i++){
         if((i>=4)&&(i<8)){
        X2[i-4]=0.0;
        }
         for(j=0;j<11;j++){
    if((i>=4)&&(i<8)){
                X2[i-4] +=A2[i-4][j]*V[j];
            }

    }
    if((i>=4)&&(i<8)){
            V2[i-4]=X2[i-4];

		MPI_Send(&V2[i-4], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

         }

}

}
if(world_rank==3){
	for(i=0;i<11;i++){


    	MPI_Recv(&V[i], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	}


for (i=0;i<11;i++){
         if(i>=8){
        X3[i-8]=0.0;
        }
         for(j=0;j<11;j++){
            if(i>=8){
                X3[i-8] +=A3[i-8][j]*V[j];
            }

    }
    if(i>=8){
            V3[i-8]=X3[i-8];

		MPI_Send(&V3[i-8], 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);

         }

}

}
 // Synchronize
    MPI_Barrier(MPI_COMM_WORLD);
//TODO-send updated score back to master
//Master
//TODO-Gather individual updates from workers, update the global weight determined by index_i
if(world_rank==0){
cnt=1;
for(i=0;i<11;i++){
if(i<=3){
MPI_Recv(&V1[i], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
if((i>3)&&(i<8)){
MPI_Recv(&V2[i-4], 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
if(i>=8){
MPI_Recv(&V3[i-8], 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}
}
for (i=0;i<11;i++){
    if(i<4){
        V[i]=V1[i];

    }
    if((i>=4)&&(i<8)){
        V[i]=V2[i-4];

    }
    if(i>=8){
        V[i]=V3[i-8];

    }
    R[i]=V[i]-K[i];

}


if(vector_norm(R,11)<=0.001){
            cnt++;
        }




if(cnt==2){
	//send stop signals and house keeping

    normalized=K[0]+K[1]+K[2]+K[3]+K[4]+K[5]+K[6]+K[7]+K[8]+K[9]+K[10];
    printf(" pageRank of each nodes normalized  \n");

    for (i=0;i<11;i++){
        printf("node %d %f\n",i, K[i]/normalized);
    }


	res=1;
	MPI_Abort(MPI_COMM_WORLD, 1);
	MPI_Barrier(MPI_COMM_WORLD);
}else{
	vnorm = vector_norm(V, 11);
    for(i=0;i<11;i++){

       V[i] = V[i] / vnorm;

    }

	for (i=0;i<11;i++){
            if(i<4){
                K[i]=V[i];
		  //K[i]=K1[i];
            }
            if((i>=4)&&(i<8)){
                K[i]=V[i];
		  //K[i]=K2[i-4];
            }
            if(i>=8){
                K[i]=V[i];
		  //K[i]=K3[i-8];
            }

	//send global weights to workers
	MPI_Send(&V[i], 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
	MPI_Send(&V[i], 1, MPI_DOUBLE, 2, 0, MPI_COMM_WORLD);
	MPI_Send(&V[i], 1, MPI_DOUBLE, 3, 0, MPI_COMM_WORLD);

	}

	k++;



}

}
MPI_Barrier(MPI_COMM_WORLD);

}
MPI_Finalize();
//return 0;

}

/*Function to Compute 2-norm of a vector */
double vector_norm(double vector[], int n)
{
    int i;
    double sum = 0 ;

    for(i=0; i<n; i++)
    {
        sum+= vector[i] * vector[i];
    }

    return sqrt(sum);
}

