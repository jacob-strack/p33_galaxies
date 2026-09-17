#include "s1p.h" 
#include <mpi.h>
#include<chrono>

using namespace std; 

vector<double> my_derived_field(vector<double> x_arr, va_list args_in);
vector<double> Q_int(vector<double> Bx, va_list args_in);
vector<double> U_int(vector<double> Bx, va_list args_in);

int main(int argc, char *argv[]){
    auto start = std::chrono::steady_clock::now(); 
    int globalrank, localrank, num_nodes; 
    int nside = 128; 
    grid localgrids[3000]; 
    MPI_Init(&argc, &argv); 
    MPI_Comm nodecomm, mastercomm; 
    MPI_Comm_rank(MPI_COMM_WORLD, &localrank); 
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes); 
    //MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, globalrank, MPI_INFO_NULL, &nodecomm); 
    //MPI_Comm_size(nodecomm, &num_nodes); 
    //MPI_Comm_split(MPI_COMM_WORLD, localrank, globalrank, &mastercomm); 
    //MPI_Comm_rank(nodecomm, &localrank); 
    std::cout << "num nodes " << localrank << " " << num_nodes << std::endl;
    int grids_this_proc;
    if(localrank==0){
       parse_hierarchy("DD0070/DD0070.hierarchy", localgrids, globalrank/num_nodes);
       Distribute_Grids(localgrids, num_nodes, globalrank, "DD0070/DD0070.hierarchy"); 
       //At this point localgrids is filled for the current node, ready for test case
       for(int grid_num = 0; grid_num < 3000; grid_num++){
            build_isRefined(localgrids + grid_num);
        }
       vector<int> unrefined_grids; 
       for(int grid_num = 0; grid_num < 3000; grid_num++){
            set_not_refined(localgrids + grid_num);
            if((localgrids+grid_num)->GetNumNotRefined() > 0){unrefined_grids.push_back(grid_num);}
        }
        int num = GetNumNotRefinedGrids(localgrids,3000); 
        int total_NR = GetTotalNotRefined(localgrids,3000);
        grids_this_proc = total_NR / num_nodes;
        for(int i = 1; i < num_nodes; i++) 
            MPI_Send(&grids_this_proc, 1, MPI_INT, i, 0, MPI_COMM_WORLD); 
    }
    std::cout << "done with init" << std::endl; 
    if(localrank > 0) 
        MPI_Recv(&grids_this_proc, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    std::cout << "array size " << grids_this_proc << std::endl;
    vector<double> rec_array(grids_this_proc);
    vector<double> rec_x_array(grids_this_proc);
    vector<double> rec_y_array(grids_this_proc);
    vector<double> rec_z_array(grids_this_proc);
    vector<double> rec_dx_array(grids_this_proc); 
    vector<double> rec_dy_array(grids_this_proc); 
    vector<double> rec_dz_array(grids_this_proc); 
    if(localrank==0){
        //Now i can split the work amongst local processors, each proc gets roughly 
        //Nzones / nproc 
        //it wont be exact since im not splitting down a grid 
        int total_NR = GetTotalNotRefined(localgrids,3000);
        flat_array Density(total_NR); 
        flat_array dx(total_NR); 
        flat_array dy(total_NR); 
        flat_array dz(total_NR); 
        flat_array derived_test(total_NR);
        dx.Makedxyz(localgrids, 3000,0); 
        dy.Makedxyz(localgrids, 3000, 1); 
        dz.Makedxyz(localgrids, 3000, 2);
        flat_array x_arr(total_NR); 
        flat_array y_arr(total_NR); 
        flat_array z_arr(total_NR); 
        flat_array Bx(total_NR);
        flat_array By(total_NR);
        flat_array Bz(total_NR);
        x_arr.Makexyz(localgrids,3000,0); 
        x_arr.SetFieldName("x");
        y_arr.Makexyz(localgrids,3000,1); 
        y_arr.SetFieldName("y");
        z_arr.Makexyz(localgrids,3000,2); 
        z_arr.SetFieldName("z");
        vector<vector<double>> dxyz(3); 
        vector<vector<double>> xyz(3); 
        vector<double> projax_arr(3);
        xyz[0] = x_arr.GetFieldData();
        xyz[1] = y_arr.GetFieldData();
        xyz[2] = z_arr.GetFieldData();
        dxyz[0] = dx.GetFieldData(); 
        dxyz[1] = dy.GetFieldData(); 
        dxyz[2] = dz.GetFieldData(); 
        Density.SetPrimative(localgrids, 1000, "Density");
        //cout << "Setting Bx" << endl; 
        //Bx.SetPrimative(localgrids, 1000, "Bx");
        //cout << "Set Bx" << endl;
        //cout.flush();
        //By.SetPrimative(localgrids, 1000, "By");
        //Bz.SetPrimative(localgrids, 1000, "Bz");
        float projax[3] = {0, 0, 1}; 
        projax_arr[0] = projax[0]; 
        projax_arr[1] = projax[1]; 
        projax_arr[2] = projax[2]; 
        vector<double> center = {0.25, 0.25, 0.25};
        vector<double> send_array = Density.GetFieldData(); 
        vector<double> send_x_array = x_arr.GetFieldData(); 
        vector<double> send_y_array = y_arr.GetFieldData(); 
        vector<double> send_z_array = z_arr.GetFieldData(); 
        vector<double> send_dx_array = dx.GetFieldData(); 
        vector<double> send_dy_array = dy.GetFieldData(); 
        vector<double> send_dz_array = dz.GetFieldData(); 
        for(int i = 0; i < grids_this_proc; i++){
            rec_array[i] = send_array[i]; 
            rec_x_array[i] = send_x_array[i]; 
            rec_y_array[i] = send_y_array[i]; 
            rec_z_array[i] = send_z_array[i]; 
            rec_dx_array[i] = send_dx_array[i]; 
            rec_dy_array[i] = send_dy_array[i]; 
            rec_dz_array[i] = send_dz_array[i]; 
        }
        //send array parts to other localranks 
        for(int send_rank = 1; send_rank < num_nodes; send_rank++){
            if(total_NR % 2 == 1 && send_rank == (num_nodes-1)){grids_this_proc++;} //capture last grid if odd number
            MPI_Send(&send_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
            MPI_Send(&send_x_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
            MPI_Send(&send_y_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
            MPI_Send(&send_z_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
            MPI_Send(&send_dx_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
            MPI_Send(&send_dy_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
            MPI_Send(&send_dz_array[send_rank*grids_this_proc], grids_this_proc, MPI_DOUBLE, send_rank, 0, MPI_COMM_WORLD); 
        }
    }
        std::cout << "test local rank " << localrank << std::endl;
        if(localrank > 0){
            MPI_Recv(&rec_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            MPI_Recv(&rec_x_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            MPI_Recv(&rec_y_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            MPI_Recv(&rec_z_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            MPI_Recv(&rec_dx_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            MPI_Recv(&rec_dy_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE); 
            MPI_Recv(&rec_dz_array[0], grids_this_proc, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,MPI_STATUS_IGNORE);
        }
        vector<double> projax_arr(3);
        vector<double> center = {0.25, 0.25, 0.25};
        float projax[3] = {0, 0, 1}; 
        projax_arr[0] = projax[0]; 
        projax_arr[1] = projax[1]; 
        projax_arr[2] = projax[2]; 
        Healpix_Map<double> global_res(nside, RING, SET_NSIDE);
        Healpix_Base base(nside, RING, SET_NSIDE); 
        global_res.fill(0); 
        vector<vector<double>> rec_xyz = {rec_x_array, rec_y_array, rec_z_array}; 
        vector<vector<double>> rec_dxyz = {rec_dx_array, rec_dy_array, rec_dz_array};
        std::cout << "start proj rank " << localrank << std::endl;
        vector<Healpix_Map<double>> res = project(rec_array, rec_xyz, rec_dxyz, center, projax_arr,"Density_map.txt", nside, 5, 0.2);
        std::cout << "finish proj rank " << localrank << std::endl;
        MPI_Reduce(&res[0][0], &global_res[0], nside2npix(nside), MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD); 
        if(localrank == 0){
            ofstream outfile("Density_map.txt"); 
            for(int i = 0; i < 12*nside*nside; i++){
                outfile << global_res[i] << endl;
            }
            outfile.close();
        }
        //MPI_Comm_free(&nodecomm);
        //MPI_Comm_free(&mastercomm); 
        auto end = std::chrono::steady_clock::now(); 
        std::chrono::duration<double> elapsed = end - start; 
        std::cout << "Total time rank " << localrank << " [s] = " << elapsed.count() << std::endl;
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
    vector<double> By = va_arg(args_in,vector<double>);
    vector<double> Bz = va_arg(args_in,vector<double>);
    vector<double> density = va_arg(args_in,vector<double>);
    vector<double> projax = va_arg(args_in,vector<double>);
    for(int i = 0; i < 3; i++)
        projax[i] /= sqrt(projax[0]*projax[0] + projax[1]*projax[1] + projax[2]*projax[2]); 
    vector<double> p = {projax[0], projax[1], projax[2]}; 
    vector<vector<double>> B_arr(3); 
    B_arr[0] = Bx; 
    B_arr[1] = By; 
    B_arr[2] = Bz;
    vector<vector<double>> B_new = rotate(B_arr, p); 
    vector<double> B_new_sq = B_new[0]*B_new[0] + B_new[1]*B_new[1] + B_new[2]*B_new[2];
    return density*(B_new[0]*B_new[0] + B_new[1]*B_new[1]) / B_new_sq;
}

vector<double> U_int(vector<double> Bx, va_list args_in){
    vector<double> By = va_arg(args_in,vector<double>);
    vector<double> Bz = va_arg(args_in,vector<double>);
    vector<double> density = va_arg(args_in,vector<double>);
    vector<double> projax = va_arg(args_in,vector<double>);
    for(int i = 0; i < 3; i++)
        projax[i] /= sqrt(projax[0]*projax[0] + projax[1]*projax[1] + projax[2]*projax[2]); 
    vector<double> p = {projax[0], projax[1], projax[2]}; 
    vector<vector<double>> B_arr(3); 
    B_arr[0] = Bx; 
    B_arr[1] = By; 
    B_arr[2] = Bz;
    vector<vector<double>> B_new = rotate(B_arr, p);
    vector<double> B_new_sq = B_new[0]*B_new[0] + B_new[1]*B_new[1] + B_new[2]*B_new[2];
    return 2.0*density*B_new[0]*B_new[1] / B_new_sq;
}
