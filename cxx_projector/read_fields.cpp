//make grid object read fields from the baryon field file that each grid knows
//read only relevant quantities, even if more flags are part of grid members 
//magnetic field, density, position

#include<fstream>
#include "simple_grid.h"
#include <iostream>
#include <mpi.h>
#include <hdf5.h>

double my_derived_field(int index, double* in1, va_list args);

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
       for(int grid_num = 0; grid_num < 2000; grid_num++){
            build_isRefined(localgrids + grid_num);
        }
       for(int grid_num = 0; grid_num < 2000; grid_num++){
            set_not_refined(localgrids + grid_num);
        }
       for(int i = 0; i < 2000; i++){
            if(localgrids[i].GetNextGridThisLevelID() != 0 && localgrids[i].GetNextGridThisLevel() == NULL){cout << "Bad pointer setup " << localgrids[i].GetGridID() << endl;}
        }
       int num = GetNumNotRefinedGrids(localgrids, 2000); 
       cout << "num " << num << endl;
       int total_NR = GetTotalNotRefined(localgrids, 2000);
       cout << "total_NR " << total_NR << endl;
       cout << "About to make flat array object" << endl;
       flat_array FA_test(total_NR);
       flat_array dx(total_NR); 
       flat_array dy(total_NR); 
       flat_array dz(total_NR); 
       dx.Makedxyz(localgrids, 2000,0); 
       dy.Makedxyz(localgrids, 2000, 1); 
       dz.Makedxyz(localgrids, 2000, 2);
       flat_array FA_2(total_NR);
       FA_2.MakeCellVolume(&dx, &dy, &dz); 
       cout << "Setting Primative" << endl; 
       FA_test.SetPrimative(localgrids, 2000, "Density");
       flat_array mass = FA_2*FA_test;
       cout << "dx sum: " << dx.GetTotal() << endl;
       cout << "CV sum: " << FA_2.GetTotal() << endl;   
       cout << "density: " << FA_test.GetDataAtInd(10) << endl; 
       cout << "cell volume: " << FA_2.GetDataAtInd(10) << endl;
       cout << "mass arr test: " << mass.GetDataAtInd(10) << endl;
       mass.WriteData("mass_arr.h5");
       cout << "total mass: " << mass.GetTotal() << endl;
       flat_array x_arr(total_NR); 
       cout << "About to ccall makexyz FA" << endl;
       x_arr.Makexyz(localgrids, 2000, 0);
       x_arr.SetFieldName("x");
       flat_array y_arr(total_NR); 
       y_arr.Makexyz(localgrids, 2000, 1);
       x_arr.SetFieldName("y");
       flat_array xpy(total_NR);
       flat_array garbo(total_NR);
       xpy.SetDerivedFlatArray(&my_derived_field, x_arr.GetFieldData(), y_arr.GetFieldData());
       int res = plot_array(x_arr, y_arr, FA_test, 1000, 1000, 0.48, 0.52, 0.48, 0.52, 1);
       cout << FA_test.GetTotal() << endl;
    }
    cout << "done doing stuff" << endl;
    //localgrids[1000].PrintAllData();
    MPI_Comm_free(&mastercomm); 
    MPI_Finalize();
    cout << "done with everything" << endl;
}


double my_derived_field(int index, double* x_arr, va_list args_in){
    double* y_arr = va_arg(args_in, double*);
    double ans; 
    ans =  x_arr[index] + y_arr[index];
    return ans; 
}

