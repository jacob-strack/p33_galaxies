#include <iostream>
#include <omp.h>
#include <mpi.h>
#include <math.h>
#include <random>
#include <ctime>
using namespace std; 
int main(int argc, char *argv[])
{
    int globalrank,localrank, nproc; 
    MPI_Init(&argc, &argv);
    MPI_Comm nodecomm, mastercomm; //communicators for master and node level
    MPI_Comm_rank(MPI_COMM_WORLD, &globalrank); 
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalrank,MPI_INFO_NULL,&nodecomm); //create node communicator split from world communicator
    MPI_Comm_split(MPI_COMM_WORLD, localrank, globalrank, &mastercomm);//create master communicator split from world
    MPI_Comm_rank(mastercomm, &localrank); //get local rank
    MPI_Comm_free(&nodecomm);
    //need to make output happen on just one processor from each node with if statement w/ localrank
    if(localrank == 0){
        cout << "local rank is 0! " << " globalrank: " << globalrank << endl;
    }
    cout << "My local rank: " << localrank << " globalrank: " << globalrank << endl;
    omp_set_num_threads(1); //number of threads for each processor
    int count = 10000; //just keep adding zeros until int overflows huh  
    int count_in = 0;
    int final_count, final_ans; 
    int ans = 0; //global values to be used with MPI Reduce
    int total_count = 0; 
    srand (static_cast <unsigned> (time(0))); //seed random
    #pragma omp parallel for 
    for(int i=0; i<=count; i++)
    {
        srand (static_cast <unsigned> (time(0))); //I think the seeding should happen here
        float rand_x = -1.0 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2.0))); //random (x,y)
        float rand_y = -1.0 + static_cast <float> (rand()) / (static_cast <float> (RAND_MAX/(2.0)));
        if(rand_x*rand_x + rand_y*rand_y <= 1.0){
            count_in++; //increment if within circle
        }
    }
    //1) Reduce from nodecomm to mastercomm
    MPI_Reduce(&count_in, &ans, 1, MPI_INT, MPI_SUM, 0, mastercomm); 
    MPI_Reduce(&count, &total_count, 1, MPI_INT, MPI_SUM, 0, mastercomm); 
    //2) Reduce from mastercomm to world
    MPI_Reduce(&ans, &final_ans, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    MPI_Reduce(&total_count, &final_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD); 
    //calculate approximation once on rank 0
    if(globalrank == 0){
        float pi = 4.0 * static_cast<float>(final_ans) / static_cast<float>(final_count); //calculate answer on global rank zero only
        cout << "points in circle: " << ans << endl; 
        cout << "total points in square: " << total_count << endl; 
        cout << pi << endl;
    }
    MPI_Comm_free(&mastercomm);
    MPI_Finalize(); //get out of MPI process
}
