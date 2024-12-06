//make grid object read fields from the baryon field file that each grid knows
//read only relevant quantities, even if more flags are part of grid members 
//magnetic field, density, position

#include<fstream>
#include "simple_grid.h"
#include <iostream>
#include <mpi.h>
#include <hdf5.h>

//int read_dataset(grid localgrids[], int size);

int main(int argc, char *argv[]){
    int globalrank, localrank, num_nodes; 
    grid localgrids[2000]; 
    MPI_Init(&argc, &argv); 
    MPI_Comm nodecomm, mastercomm; 
    MPI_Comm_rank(MPI_COMM_WORLD, &globalrank); 
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalrank, MPI_INFO_NULL, &nodecomm); 
    MPI_Comm_size(nodecomm, &num_nodes); 
    MPI_Comm_split(MPI_COMM_WORLD, localrank, globalrank, &mastercomm); 
    MPI_Comm_rank(nodecomm, &localrank); 
    MPI_Comm_free(&nodecomm);
    if(localrank==0){
        parse_hierarchy("time0000.hierarchy", localgrids, globalrank/num_nodes);
       Distribute_Grids(localgrids, num_nodes, globalrank, "time0000.hierarchy"); 
       //At this point localgrids is filled for the current node, ready for test case 
       read_dataset(localgrids, 2000, globalrank); 
    }
    MPI_Comm_free(&mastercomm); 
    MPI_Finalize();
    cout << "Finished" << endl;
}

/*int read_dataset(grid localgrids[], int num_grids){
    //get length of localgrids
    cout << "localgrids length: " << num_grids << endl; 
    //loop over localgrids
    for(int i = 0; i < num_grids; i++){
        char *filename = localgrids[i].GetBaryonFileName(); 
        int *FieldType = localgrids[i].GetFieldType(); 
        int *dims = localgrids[i].GetGridDimension();
        int rank = localgrids[i].GetGridRank();
        int GridID = localgrids[i].GetGridID();
        int densnum = 0; //field labels following enzo convention 
        int B1num = 49; 
        int B2num = 50; 
        int B3num = 51; 
        int labels[4] = {densnum, B1num, B2num, B3num};  
        if(GridID == 0){continue;}
        //checking that all the field labels that I expect exist and were read from .hierarchy
        for(int i = 0; i < 4; i++){
            for(int j = 0; j < 100; j++){
                if(FieldType[j] == labels[i]){break;}
                if(j == 99){
                    cout << "CANT FIND FIELD LABEL" << endl;
                    return 0;
                }
            }
        }
        const char *field_name = "Density";
        char group_name[100]; 
        sprintf(group_name,"Grid%08d",GridID); 
        cout << "Group name: " << group_name << endl;
        //Now field labels are sure to exist, start reading fields from .cpu files
        hid_t file_id, group_id, dset_id; 
        file_id = H5Fopen(filename, H5F_ACC_RDONLY,H5P_DEFAULT);
        group_id = H5Gopen2(file_id, group_name,H5P_DEFAULT);
        int size = 1; 
        for(int j = 0; j < rank; j++){size *= dims[j];}
        for(int field = 0; field < 4; field++){
            localgrids[i].BaryonField[field] = new float[size];
            for(int k = 0; k < size; k++){
                localgrids[i].BaryonField[field][k] = 0; 
            }
        }
        dset_id = H5Dopen(group_id, field_name); 
    }
    return 1;
}*/
