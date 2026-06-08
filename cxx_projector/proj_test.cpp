#include"proj.h"

using namespace std; 
int main(){
    vector<double> cb = make_cube(8); 
    vector<vector<double>>a = make_xyz(8);
    float projax[3] = {0,0,1}; 
    float center[3] = {0.5,0.5,0.5};
    vector<vector<double>> b = make_phi_theta(a, projax, center); 
    for(int i = 0; i < a[0].size(); i++){
        cout << a[0][i] << " " << a[1][i] << " " << a[2][i] << " " << b[0][i] << " " << b[1][i] << endl; 
    }
    
    for(int i = 0; i < a[0].size(); i++){
        if(cb[i] != 0){
            cout << cb[i] << " " << a[0][i] << " " << a[1][i] << " " << a[2][i] << endl;
        }
    }
    cout << "done" << endl;
    return 0;
}
