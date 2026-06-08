#include "s1p.h" 

using namespace std; 

int main(){
    int size = 8;
    vector<double> ahh = make_wire(size);
    cout << "made wire" << endl; 
    cout.flush();
    vector<double> a = make_cube(size); 
    cout << "made cube" << endl; 
    cout.flush();
    vector<vector<double>> b = make_xyz(size); 
    vector<vector<double>> c = make_dxyz(size); 
    float projax[3] = {1.,0.,0.}; 
    float center[3] = {0.0,0.0,0.5}; 
    cout << "about to project" << endl; 
    cout.flush();
    vector<Healpix_Map<double>> d = project(a, b, c, center, projax, "plot_array.txt", 32, 2);  
    
    return 0; 
}
