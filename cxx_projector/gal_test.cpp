#include "s1p.h" 
#include <mpi.h>

using namespace std; 

vector<double> my_derived_field(vector<double> x_arr, va_list args_in);
vector<double> Q_int(vector<double> Bx, va_list args_in);
vector<double> U_int(vector<double> Bx, va_list args_in);

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
       parse_hierarchy("DD0070/DD0070.hierarchy", localgrids, globalrank/num_nodes);
       Distribute_Grids(localgrids, num_nodes, globalrank, "DD0070/DD0070.hierarchy"); 
       //At this point localgrids is filled for the current node, ready for test case
       for(int grid_num = 0; grid_num < 2000; grid_num++){
            build_isRefined(localgrids + grid_num);
        }
       for(int grid_num = 0; grid_num < 2000; grid_num++){
            set_not_refined(localgrids + grid_num);
        }
        int num = GetNumNotRefinedGrids(localgrids,2000); 
        int total_NR = GetTotalNotRefined(localgrids,2000); 
        flat_array Density(total_NR); 
        flat_array dx(total_NR); 
        flat_array dy(total_NR); 
        flat_array dz(total_NR); 
        flat_array derived_test(total_NR);
        dx.Makedxyz(localgrids, 2000,0); 
        dy.Makedxyz(localgrids, 2000, 1); 
        dz.Makedxyz(localgrids, 2000, 2);
        flat_array x_arr(total_NR); 
        flat_array y_arr(total_NR); 
        flat_array z_arr(total_NR); 
        flat_array Bx(total_NR);
        flat_array By(total_NR);
        flat_array Bz(total_NR);
        x_arr.Makexyz(localgrids,2000,0); 
        x_arr.SetFieldName("x");
        y_arr.Makexyz(localgrids,2000,1); 
        y_arr.SetFieldName("y");
        z_arr.Makexyz(localgrids,2000,2); 
        z_arr.SetFieldName("z");
        vector<vector<double>> dxyz(3); 
        vector<vector<double>> xyz(3); 
        vector<double> projax_arr(3);
        xyz[0] = x_arr.GetFieldData();
        xyz[1] = y_arr.GetFieldData();
        xyz[2] = z_arr.GetFieldData();
        dxyz[0] = dx.GetFieldData(); 
        dxyz[1] = dx.GetFieldData(); 
        dxyz[2] = dx.GetFieldData(); 
        Density.SetPrimative(localgrids, 2000, "Density");
        cout << "Setting Bx" << endl; 
        Bx.SetPrimative(localgrids, 2000, "Bx");
        cout << "Set Bx" << endl;
        cout.flush();
        By.SetPrimative(localgrids, 2000, "By");
        Bz.SetPrimative(localgrids, 2000, "Bz");
        float projax[3] = {0, 1, 0}; 
        projax_arr[0] = projax[0]; 
        projax_arr[1] = projax[1]; 
        projax_arr[2] = projax[2]; 
        float center[3] = {0.25,0.25,0.25};
        cout << "About to setderived" << endl;
        derived_test.SetDerivedFlatArray(&U_int, Bx.GetFieldData(), By.GetFieldData(), Bz.GetFieldData(), Density, projax_arr); 
        cout << "end of setderived" << endl; 
        cout.flush();
        //vector<Healpix_Map<double>> res = project(derived_test.GetFieldData(), xyz, dxyz, center, projax,"U_map.txt", 32, 5);
        vector<Healpix_Map<double>> res = project(Density.GetFieldData(), xyz, dxyz, center, projax,"Density_map.txt", 64, 1);
    }
        MPI_Comm_free(&mastercomm); 
        MPI_Finalize();

        return 0; 
}

vector<double> my_derived_field(vector<double> x_arr, va_list args_in){
    vector<double> y_arr = va_arg(args_in, vector<double>); 
    vector<double> ans;
    ans = x_arr + y_arr;
    return ans; 
}

vector<double> Q_int(vector<double> Bx, va_list args_in){
    float no_center[3] = {0.0, 0.0, 0.0};
    vector<double> By = va_arg(args_in,vector<double>);
    vector<double> Bz = va_arg(args_in,vector<double>);
    vector<double> density = va_arg(args_in,vector<double>);
    vector<double> projax = va_arg(args_in,vector<double>);
    for(int i = 0; i < 3; i++)
        projax[i] /= sqrt(projax[0]*projax[0] + projax[1]*projax[1] + projax[2]*projax[2]); 
    float p[3] = {(float)projax[0], (float)projax[1], (float)projax[2]}; //lol
    vector<vector<double>> B_arr(3); 
    B_arr[0] = Bx; 
    B_arr[1] = By; 
    B_arr[2] = Bz;
    vector<vector<double>> B_new = rotate(B_arr, p, no_center); 
    vector<double> B_new_sq = B_new[0]*B_new[0] + B_new[1]*B_new[1] + B_new[2]*B_new[2];
    return density*(B_new[0]*B_new[0] + B_new[1]*B_new[1]) / B_new_sq;
}

vector<double> U_int(vector<double> Bx, va_list args_in){
    float no_center[3] = {0.0, 0.0, 0.0};
    vector<double> By = va_arg(args_in,vector<double>);
    vector<double> Bz = va_arg(args_in,vector<double>);
    vector<double> density = va_arg(args_in,vector<double>);
    vector<double> projax = va_arg(args_in,vector<double>);
    for(int i = 0; i < 3; i++)
        projax[i] /= sqrt(projax[0]*projax[0] + projax[1]*projax[1] + projax[2]*projax[2]); 
    float p[3] = {(float)projax[0], (float)projax[1], (float)projax[2]}; //lol
    vector<vector<double>> B_arr(3); 
    B_arr[0] = Bx; 
    B_arr[1] = By; 
    B_arr[2] = Bz;
    vector<vector<double>> B_new = rotate(B_arr, p, no_center);
    vector<double> B_new_sq = B_new[0]*B_new[0] + B_new[1]*B_new[1] + B_new[2]*B_new[2];
    return 2.0*density*B_new[0]*B_new[1] / B_new_sq;
}
