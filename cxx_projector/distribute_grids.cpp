//code that creates grid objects on each node, but only actually read in 
//grids on certain nodes to balance the workload 
#include <iostream>
#include <fstream>
#include "simple_grid.h"
#include <mpi.h>  

//int Distribute_Grids(int num_nodes, int total_grids, int global_rank, const char* filename); 

int main(int argc, char *argv[]){
	int globalrank, localrank; 
    int num_nodes;
    ifstream file("time0000.hierarchy"); 
	int total_grids = 10; 
    grid localgrids[2000];
	MPI_Init(&argc, &argv);
	MPI_Comm nodecomm, mastercomm; //communicators for master and node level
	MPI_Comm_rank(MPI_COMM_WORLD, &globalrank); //get global rank 
	MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalrank,MPI_INFO_NULL,&nodecomm); //create node communicator split from world communicator
	MPI_Comm_size(nodecomm, &num_nodes);
	MPI_Comm_split(MPI_COMM_WORLD, localrank, globalrank, &mastercomm);//create master communicator split from world
	MPI_Comm_rank(nodecomm, &localrank); //get local rank
	MPI_Comm_free(&nodecomm);
    //loop to get the total number of grids from the file
    /*int TestGridID = -1;
    int CurrentGridID = -1; 
    char ch; 
    string line; 
    while(getline(file, line)){
            if(sscanf(line.c_str(), "\nGrid = %d\n", &TestGridID) == 1){
               if(TestGridID > CurrentGridID){
                    CurrentGridID = TestGridID;
                }
            }
    }*/
    //total_grids = CurrentGridID;
	if(localrank == 0){
		cout << "num nodes: " << num_nodes << endl;
		cout << "node number: " << globalrank / num_nodes << endl;
		Distribute_Grids(localgrids, num_nodes, globalrank, "time0000.hierarchy");  
        for(int i = 0; i < 2000; i++){
            localgrids[i].PrintAllData();
        }
	}
	MPI_Comm_free(&mastercomm);
	MPI_Finalize();
	return 0;

}

/*int Distribute_Grids(int num_nodes, int total_grids, int node_num, const char *filename){
	//determining the total number of grids 
	string line;
	ifstream file(filename);
	int TestGridID = -1;
	int CurrentGridID = -1;  
    while(getline(file, line)){
            if(sscanf(line.c_str(), "\nGrid = %d\n", &TestGridID) == 1){
               if(TestGridID > CurrentGridID){
                    CurrentGridID = TestGridID;
                }
            }
    }
	total_grids = CurrentGridID;
	cout << "total grids in function: " << total_grids << endl; 
	grid grids[total_grids]; //array of grids
	int lower_grid = node_num * (total_grids / num_nodes) + 1; //lower bound of indices of grids to read on each proc 
	int upper_grid = (node_num + 1) * (total_grids / num_nodes); //upper bound of indices of grids to read on each proc
    if(total_grids % 2 == 1 and node_num == num_nodes - 1){ //catch last grid if odd number of total grids
        upper_grid++; 
    }
    cout << "lower grid num: " << lower_grid << endl; 
    cout << "upper_grid num: " << upper_grid << endl;
	parse_hierarchy(filename, lower_grid, upper_grid, grids, total_grids, node_num);
	return 1; 
}*/
