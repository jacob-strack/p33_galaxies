#include <iostream>
#include <omp.h> 
#include <mpi.h>
#include <string>
#include <fstream>
#include "simple_grid.h"
using namespace std; 

//code to parse hierarchy file
//copying enzo's homework pretty heavily
//ReadAllData.C
//need to read Left,Right,StartIndex,EndIndex,NextGrid

//all functions have been incorporated into simple_grid.h 

int main(int argc, char *argv[]){
	//put mpi comm splits here and go only on localrank = 0
	int globalrank, localrank;  
	MPI_Init(&argc, &argv);
	MPI_Comm nodecomm, mastercomm; //communicators for master and node level
	MPI_Comm_rank(MPI_COMM_WORLD, &globalrank); //get global rank 
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalrank,MPI_INFO_NULL,&nodecomm); //create node communicator split from world communicator
	MPI_Comm_split(MPI_COMM_WORLD, localrank, globalrank, &mastercomm);//create master communicator split from world
	MPI_Comm_rank(nodecomm, &localrank); //get local rank
	MPI_Comm_free(&nodecomm);
	if(localrank == 0){
		grid localgrid[5];
		int ids[5] = {1,2,3,4,5}; 
		//I only want to read in the hierarchy once for each node 
		int check;
		const char *filename = "time0000.hierarchy";	
		check = parse_hierarchy(filename, ids, localgrid, sizeof(ids)/sizeof(ids[0])); 
		cout << "check: " << check << endl;
	}
	MPI_Comm_free(&mastercomm); 
	MPI_Finalize();
}