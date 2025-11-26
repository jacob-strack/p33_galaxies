#include "s1p.h"
#include <mpi.h>

using namespace std; 

int main(int argc, char *argv[]){ 
    
    int globalrank, localrank, num_nodes; 
    /*
    MPI_Init(&argc, &argv); 
    MPI_Comm nodecomm, mastercomm; 
    MPI_Comm_rank(MPI_COMM_WORLD, &globalrank); 
    MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalrank, MPI_INFO_NULL, &nodecomm); 
    MPI_Comm_size(nodecomm, &num_nodes); 
    MPI_Comm_split(MPI_COMM_WORLD, localrank, globalrank, &mastercomm); 
    MPI_Comm_rank(nodecomm, &localrank); 
    MPI_Comm_free(&nodecomm);
    */

    flat_array test_arr; 
    test_arr.SetPrimativeAmrex("plt35000", "gasDensity"); 
    
    flat_array x_arr, y_arr, z_arr;
    flat_array dx_arr, dy_arr, dz_arr; 
    
    x_arr.SetPrimativeAmrex("plt35000", "x",0);
    y_arr.SetPrimativeAmrex("plt35000", "y",0);
    z_arr.SetPrimativeAmrex("plt35000", "z",0);
    dx_arr.SetPrimativeAmrex("plt35000", "dx",0);
    dy_arr.SetPrimativeAmrex("plt35000", "dy",0);
    dz_arr.SetPrimativeAmrex("plt35000", "dz",0);

    vector<vector<double>> xyz(3); 
    vector<vector<double>> dxyz(3);
    vector<double> projax_arr(3); 

    xyz[0] = x_arr.GetFieldData();
    xyz[1] = y_arr.GetFieldData();
    xyz[2] = z_arr.GetFieldData();
    dxyz[0] = dx_arr.GetFieldData(); 
    dxyz[1] = dy_arr.GetFieldData(); 
    dxyz[2] = dz_arr.GetFieldData(); 
    
    float projax[3] = {1, 0, 0}; 
    float center[3] = {0.5, 0.5, 0.5}; //this code units? 

    projax_arr[0] = projax[0]; 
    projax_arr[1] = projax[1]; 
    projax_arr[2] = projax[2]; 
    
    cout << *min_element(xyz[0].begin(), xyz[0].end()) << " " << *max_element(xyz[0].begin(), xyz[0].end()) << endl;
    cout << *min_element(xyz[1].begin(), xyz[1].end()) << " " << *max_element(xyz[1].begin(), xyz[1].end()) << endl;
    cout << *min_element(xyz[2].begin(), xyz[2].end()) << " " << *max_element(xyz[2].begin(), xyz[2].end()) << endl;
    cout << xyz[0].size() << endl;   
    cout << xyz[1].size() << endl;   
    cout << xyz[2].size() << endl; 
    cout << x_arr.GetSize() << endl;
    vector<Healpix_Map<double>> res = project(test_arr.GetFieldData(), xyz, dxyz, center, projax, "density_proj.txt", 32,5); 
    //MPI_Comm_free(&mastercomm); 
    //MPI_Finalize();
    return 0; 
}
