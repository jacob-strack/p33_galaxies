#include<iostream> 
#include<hdf5.h>
#include<cmath>
#include<cstring>
#include"simple_grid.h"

using namespace std; 

vector<double> make_cube(int N);
vector<double> make_wire(int N); 
vector<vector<double>> make_xyz(int N);
vector<vector<double>> make_dxyz(int N);
vector<vector<double>> rotate(vector<vector<double>> xyz, float projax[3], float center[3]); 
vector<vector<double>> make_phi_theta(vector<vector<double>> xyz, float projax[3], float center[3]);

vector<double> make_cube(int N){
    double dx = 1/(double)N;
    vector<double> cube((int)pow(N,3));
    double c_1[3] = {0.2,0.8,0.4};
    double c_2[3] = {0.8, 0.1, 0.4};
    double r = 0.1;
    double x,y,z;
    //ans = new double[3][(int)pow(N,3)]; 
    for(int j = 0; j < (int)pow(N,3); j++){
        cube[j] = 0; 
        z = (0.5 + (double)(j%N)) * dx; 
        y = (0.5 + (double)((j/N)%N)) * dx; 
        x = (0.5 + (double)((j/(int)pow(N,2))%N)) * dx; 
        bool ok1 = pow(x - c_1[0],2) + pow(y - c_1[1], 2)  < pow(r,2);
        bool ok2 = pow(x - c_2[0],2) + pow(y - c_2[1], 2) < pow(r,2);
        if(ok1 && z < 0.5){
            cube[j] = 1.0; 
        }
        if(ok1 && z > 0.5){
            cube[j] = 2.0; 
        }
        if(ok2){
            cube[j] = 3.0;
        }
}
    return cube;
}

vector<double> make_wire(int N){
    double dx = 1/(double)N;
    double min = 1/(double(2*N)); 
    double max = 1 - min; 
    cout << dx << endl;
    cout << min << " " << max << endl;
    vector<double> cube((int)pow(N,3)); 
    //double ans[3][(int)pow(N,3)]; 
    double xyz[3];
    for(int j = 0; j < (int)pow(N,3); j++){
        cube[j] = 0;
        xyz[0] = (0.5 + (double)(j%N)) * dx; 
        xyz[1] = (0.5 + (double)((j/N)%N)) * dx; 
        xyz[2] = (0.5 + (double)((j/(int)pow(N,2))%N)) * dx; 
        for(int i = 0; i < 3; i++)
            for(int k = 0; k < 3; k++){ 
                if(xyz[k] == max && xyz[i] == max && i!=k)
                    cube[j] = 10.0; 
                if(xyz[k] == min && xyz[i] == min && i!=k)
                    cube[j] = 10.0;
                if(xyz[k] == max && xyz[i] == min)
                    cube[j] = 10.0;
            } 
                
    }
    return cube;

}

vector<vector<double>> make_xyz(int N){ 
    double dx = 1./N;
    vector<vector<double>> xyz(3); 
    for(int dim = 0; dim < 3; dim++){
        xyz[dim] = vector<double>((int)pow(N,3));
    }
    for(int j = 0; j < (int)pow(N,3); j++){
        xyz[2][j] = (0.5 + (double)(j%N)) * dx; 
        xyz[1][j] = (0.5 + (double)((j/N)%N)) * dx; 
        xyz[0][j] = (0.5 + (double)((j/(int)pow(N,2))%N)) * dx; 
    }
    return xyz; 
}

vector<vector<double>> make_dxyz(int N){ 
    double dx = 1./N;
    vector<vector<double>> dxyz(3);
    for(int dim = 0; dim < 3; dim++)
        dxyz[dim] = vector<double>((int)pow(N,3));
    for(int j = 0; j < (int)pow(N,3); j++){
        dxyz[2][j] = dx; 
        dxyz[1][j] = dx; 
        dxyz[0][j] = dx; 
    }
    return dxyz; 
}

vector<vector<double>> rotate(vector<vector<double>> xyz, float projax[3], float center[3]){
    vector<vector<double>> xyz_p(3); 
    for(int i = 0; i < 3; i++)
        xyz_p[i] = vector<double>(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size(); i++){
        xyz_p[0][i] = xyz[0][i] - center[0];  
        xyz_p[1][i] = xyz[1][i] - center[1];  
        xyz_p[2][i] = xyz[2][i] - center[2];  
    }
    double r_proj = sqrt(projax[0]*projax[0] + projax[1]*projax[1] + projax[2]*projax[2]); 
    for(int i = 0; i < 3; i++)
        projax[i] /= r_proj; 
    double theta = acos(projax[2]); 
    double phi = atan2(projax[1], projax[0]);
    double z_p[3] = {projax[0], projax[1], projax[2]}; 
    double y_p[3] = {-1*sin(phi), cos(phi), 0};
    double x_p[3] = {cos(theta)*cos(phi), cos(theta)*sin(phi), -1*sin(theta)};
    vector<vector<double>> xyz_new(3); 
    for(int i = 0; i < 3; i++)
        xyz_new[i] = vector<double>(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size(); i++){
        xyz_new[0][i] = xyz_p[0][i]*x_p[0] + xyz_p[1][i]*x_p[1] + xyz_p[2][i]*x_p[2]; 
        xyz_new[1][i] = xyz_p[0][i]*y_p[0] + xyz_p[1][i]*y_p[1] + xyz_p[2][i]*y_p[2]; 
        xyz_new[2][i] = xyz_p[0][i]*z_p[0] + xyz_p[1][i]*z_p[1] + xyz_p[2][i]*z_p[2]; 
    }
    return xyz_new;
}

vector<vector<double>> make_phi_theta(vector<vector<double>> xyz, float projax[3], float center[3]){
    //vector<vector<double>> xyz_n = rotate(xyz, projax, center); 
    vector<vector<double>> angles(2);
    for(int i = 0; i < 2; i++)
        angles[i] = vector<double>(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size(); i++){
        angles[1][i] = M_PI - acos(xyz[0][i] / sqrt(pow(xyz[0][i],2) + pow(xyz[1][i],2) + pow(xyz[2][i],2)));
        angles[0][i] = atan2(xyz[1][i], xyz[2][i]);
        if(angles[0][i] < 0) 
            angles[0][i] += 2*M_PI;
    }
    return angles;
}
