#include"proj.h"
#include<algorithm>
#include<healpix_base.h>
#include<healpix_map.h>
#include<healpix_map_fitsio.h>
#include<alm.h>
#include<pointing.h>
#include<Eigen/Dense>
#include<libqhullcpp/Qhull.h>
#include<libqhullcpp/QhullFacetList.h>
#include<libqhullcpp/QhullVertexSet.h>
#include<chealpix.h>
#include<set> 
#include<utility> 


vector<double> normalize_3d(vector<double> input);
vector<vector<double>> cube_face_planes(vector<vector<double>> corners, float tol = 1e-10);
vector<bool> points_in_cube(vector<vector<double>> pts, vector<vector<double>> face_normals, vector<double> face_c, float tol=1e-10);
bool origin_in_cube(vector<vector<double>> face_normals, vector<double> face_c, float tol = 1e-10);
vector<vector<double>> order_boundary_vecs_ccw(vector<vector<double>> ray_dirs, vector<double> center_vec);
vector<vector<double>> cone_side_normals(vector<vector<double>> ray_dirs);
vector<vector<vector<double>>> precompute_healpix_geometry(int nside);
bool cube_fully_inside_cone(vector<vector<double>> corners, vector<vector<double>> face_normals, float tol = 1e-10);
bool cube_fully_outside_cone(vector<vector<double>> corners, vector<vector<double>> face_normals, float tol = 1e-10);
vector<vector<bool>> classify_pixels_for_zone(vector<vector<double>> zone_corners, vector<vector<vector<double>>> cand_side_normals, float tol = 1e-10);
vector<bool> points_in_cone(vector<vector<double>> pts, vector<vector<double>> side_normals, float tol = 1e-10); 
pair<vector<vector<vector<double>>>, vector<vector<bool>>> segment_plane_intersections_batch(vector<vector<double>> p0, vector<vector<double>> p1, vector<vector<double>> plane_normals, float tol = 1e-10);
pair<vector<vector<vector<double>>>, vector<vector<bool>>> ray_plane_intersections_batch(vector<vector<double>> ray_dirs, vector<vector<double>> face_normals, vector<double> face_c, float tol = 1e-10);
vector<vector<double>> unique_points(vector<vector<double>> pts, float tol = 1e-10); 
vector<vector<double>> cube_cone_intersection_vertices_precomputed(vector<vector<double>> corners, vector<vector<double>> edge_p0, vector<vector<double>> edge_p1, vector<vector<double>> face_normals, vector<double> face_c, vector<vector<double>> ray_dirs, vector<vector<double>> side_normals, float tol = 1e-10);
double cube_cone_intersection_volume_precomputed(vector<vector<double>> corners, vector<vector<double>> edge_p0, vector<vector<double>> edge_p1, vector<vector<double>> face_normals, vector<double> face_c, vector<vector<double>> ray_dirs, vector<vector<double>> side_normals, float tol = 1e-10); 

vector<Healpix_Map<double>> project(vector<double> cube, vector<vector<double>> xyz, vector<vector<double>> dxyz, float center[3], float projax[3], const char* filename, int NSIDE, int exclude, float max_r = 1.0){
    vector<vector<int>> cube_edges = {{0,1},{1,2},{2,3},{3,0},
                                          {4,5},{5,6},{6,7},{7,4},
                                          {0,7},{1,6},{2,5},{3,4}};
    vector<vector<int>> cube_faces = {{0,1,6,7},
                                      {2,3,4,5},
                                      {0,1,2,3},
                                      {4,5,6,7},
                                      {0,3,4,7},
                                      {1,2,5,6}};
    //precompute healpix things 
    vector<vector<vector<double>>> healpix_cache = precompute_healpix_geometry(NSIDE); 
    //xyz from proj_center
    for(int dim = 0; dim < 3; dim++)
        for(int i = 0; i < xyz[0].size(); i++)
            xyz[dim][i] = xyz[dim][i] - center[dim]; 
    vector<float> r(xyz[0].size()); 
    vector<float> cell_scale(xyz[0].size()); 
    float min = 1e-30;
    for(int i = 0; i < xyz[0].size(); i++){
        cell_scale[i] = (dxyz[0][i] + dxyz[1][i] + dxyz[2][i]) / 3.0; 
        r[i] = xyz[0][i]*xyz[0][i] + xyz[1][i]*xyz[1][i] + xyz[2][i]*xyz[2][i];  
    }
    vector<float>nxyz(r.size()); 
    for(int i = 0; i < xyz[0].size(); i++){
        r[i] = sqrt(r[i]);
        nxyz[i] = r[i] / max(cell_scale[i], min); 
    }
    vector<bool> mask(r.size()); //mask out close points 
    int Nz = 0; 
    for(int i = 0; i < r.size(); i++){
        if(nxyz[i] > exclude){
            mask[i] = true;
            Nz++; 
        }
        else if(r[i] > max_r){
            mask[i] = false;
        }
        else
            mask[i] = false; 
    }
    std::cout << "Num zones: " << Nz << endl;
    Healpix_Map<double> counts(NSIDE, RING, SET_NSIDE); 
    Healpix_Map<double> output(NSIDE, RING, SET_NSIDE);
    Healpix_Base base(NSIDE, RING,SET_NSIDE);
    //fill maps with zeros initially
    counts.fill(0); 
    output.fill(0);
    //rotate axes 
    cout << "Rotation " << endl; 
    float no_center[3] = {0.0, 0.0, 0.0};//duct tape for now where i dont want the fucntion
                                         //to shift anymore bc i already did it earlier
    vector<vector<double>> xyz_p = rotate(xyz, projax, no_center); 

    //shifter object
    cout << "Shifter " << endl;
    vector<vector<vector<double>>> corners(8); 
    vector<vector<double>> center_vecs(3);
    float shifter[8][3] = {{-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5}, {0.5,-0.5,0.5},{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.5,0.5,0.5},{-0.5,0.5,0.5},{-0.5,0.5,-0.5}};
    for(int i = 0; i < 3; i++) 
        center_vecs[i] = vector<double>(xyz[0].size()); 
    for(int i = 0; i < 8; i++){
        corners[i] = vector<vector<double>>(3); 
        for(int j = 0; j < 3; j++){
            corners[i][j] = vector<double>(xyz[0].size());
            for(int k = 0; k < xyz[0].size(); k++){
                corners[i][j][k] = xyz[j][k] + dxyz[j][k]*shifter[i][j]; 
                center_vecs[j][k] = xyz_p[j][k];
            }
        }
    }
    cout << "rotate " << endl; 
    for(int i = 0; i < 8; i++)
        corners[i] = rotate(corners[i], projax, no_center);
    vector<vector<vector<double>>> corner_vecs = corners;    
    cout << "unit vectors " << endl; 
    //get the unit vectors for the centers and the corners using rotated frame 
    vector<vector<double>> dots(8);
    vector<double> circle_radius(xyz[0].size()); 
    for(int i = 0; i < 8; i++){
        dots[i] = vector<double>(xyz[0].size());
        }
    for(int i = 0; i < 8; i++)
        for(int k = 0; k < xyz[0].size(); k++){
            //make unit vectors for corners and center
            double corner_mag = sqrt(corner_vecs[i][0][k]*corner_vecs[i][0][k] + corner_vecs[i][1][k]*corner_vecs[i][1][k] + corner_vecs[i][2][k]*corner_vecs[i][2][k]); 
            double center_mag = sqrt(center_vecs[0][k]*center_vecs[0][k] + center_vecs[1][k]*center_vecs[1][k] + center_vecs[2][k]*center_vecs[2][k]); 
            for(int j = 0; j < 3; j++){
                corner_vecs[i][j][k] /= corner_mag; 
                center_vecs[j][k] /= center_mag;
            }
            dots[i][k] = corner_vecs[i][0][k]*center_vecs[0][k] + corner_vecs[i][1][k]*center_vecs[1][k] + corner_vecs[i][2][k]*center_vecs[2][k];
            if(dots[i][k] > 1.0) 
                dots[i][k] = 1.0; 
            if(dots[i][k] < -1.0) 
                dots[i][k] = -1.0;
            circle_radius[k] = acos(dots[0][k]);  
            for(int ind = 1; ind < 8; ind++)
                if(acos(dots[ind][k]) > circle_radius[k])
                    circle_radius[k] = acos(dots[ind][k]); 
        }
        

    //accumulation things 
    std::cout << "set up accumulation" << std::endl;
    vector<float> r_sq(xyz[0].size());
    vector<float> zone_volume(xyz[0].size()); 
    vector<float> zone_emission(xyz[0].size());
    for(int i = 0; i < xyz[0].size(); i++){
        r_sq[i] = xyz_p[0][i]*xyz_p[0][i] + xyz_p[1][i]*xyz_p[1][i] + xyz_p[2][i]*xyz_p[2][i]; 
        zone_volume[i] = dxyz[0][i]*dxyz[1][i]*dxyz[2][i]; 
        //zone_emission[i] = cube[i]/r_sq[i]*zone_volume[i];
        //change zone emission to cube * zone_volume and use unit vectors later 
        zone_emission[i] = cube[i] * zone_volume[i]; 
        if(mask[i] == true)
            cout << "cube " << i << " " << cube[i] << endl;
    }
    cout << "corner angles " << endl; 
    //get angles at edges of zone
    vector<vector<vector<double>>> corner_angles(8); 
    for(int i = 0; i < 8; i++){
        corner_angles[i] = vector<vector<double>>(2);
        vector<vector<double>> ans = make_phi_theta(corners[i], projax, no_center);
        for(int j = 0; j < 2; j++){
            corner_angles[i][j] = vector<double>(xyz[0].size());
            for(int k = 0; k < xyz[0].size(); k++){
                corner_angles[i][j][k] = ans[j][k];
                if(isnan(ans[j][k]) && mask[k] == true){
                    std::cout << "Found a nan " << j << " " << k << endl;
                    cout << corner_angles[i][0][k] << " " << corner_angles[i][1][k] << endl; 
                }
            }
        } 
    }
   std::cout << "About to start zone loop" << std::endl;    
   //ZONE LOOOOOOOOOOOOOOOOOOOOOP 
   int debug_ind = 0;
   cout << "Nz " << Nz << endl; 
   for(int izone = 0; izone < xyz[0].size(); izone++){
       cout << "izone " << izone << " " << xyz[0].size() << endl;
       if(mask[izone]==false){
           cout << "mask false iter" << endl; 
           continue;
        }
       double quantity = zone_emission[izone]; 
       vector<double> zone_center_vec(3); 
       vector<vector<double>> zone_corners(8); 
       double this_zone_volume = zone_volume[izone];  
       for(int i = 0; i < 3; i++) 
           zone_center_vec[i] = center_vecs[i][izone]; 
       for(int i = 0; i < 8; i++){
            zone_corners[i] = vector<double>(3);
            for(int j = 0; j < 3; j++) 
                zone_corners[i][j] = corners[i][j][izone]; 
        }
       //per-zone geometry precompute
       vector<vector<double>> cube_face_pl = cube_face_planes(zone_corners);  
       vector<vector<double>> edge_p0(12); 
       vector<vector<double>> edge_p1(12);
       for(int i = 0; i < 12; i++){ 
           edge_p0[i] = vector<double>(3); 
           edge_p1[i] = vector<double>(3); 
           for(int j = 0; j < 3; j++){
               edge_p0[i][j] = zone_corners[cube_edges[i][0]][j]; 
               edge_p1[i][j] = zone_corners[cube_edges[i][1]][j]; 
           }
        }
       //query disk call
       pointing this_pointing(vec3(zone_center_vec[0], zone_center_vec[1], zone_center_vec[2])); 
       vector<int> pix_nums;
       base.query_disc_inclusive(this_pointing, circle_radius[izone], pix_nums);
       if(pix_nums.size() == 0){
           cout << "pix size 0 iter" << endl; 
           continue;
        }
       vector<vector<vector<double>>> cand_side_normals(pix_nums.size()); 
       vector<vector<vector<double>>> ray_dirs(pix_nums.size());
       for(int i = 0; i < pix_nums.size(); i++){
           cout << "pix " << pix_nums[i] << endl;
           cand_side_normals[i] = vector<vector<double>>(4);
           ray_dirs[i] = vector<vector<double>>(4);
           for(int j = 0; j < 4; j++){
               cand_side_normals[i][j] = vector<double>(3);
               ray_dirs[i][j] = vector<double>(3);
               for(int k = 0; k < 3; k++){
                    cand_side_normals[i][j][k] = healpix_cache[pix_nums[i]][j+5][k];
                    ray_dirs[i][j][k] = healpix_cache[pix_nums[i]][j+1][k]; 
               }
           }
        }
       vector<vector<bool>> classify = classify_pixels_for_zone(zone_corners, cand_side_normals);
       for(int i = 0; i < classify[0].size(); i++){
       //add any full accept pixels 
            if(classify[0][i]){
                output[pix_nums[i]] += quantity;
                counts[pix_nums[i]] += this_zone_volume;
            }
       //exact pixels 
            if(classify[2][i]){
                double intersect_volume = 0.0; 
                if(cube_fully_inside_cone(zone_corners, cand_side_normals[i])){
                   intersect_volume = this_zone_volume;} 
                else if(cube_fully_outside_cone(zone_corners, cand_side_normals[i])){
                    intersect_volume = 0.0; }
                else{
                    vector<vector<double>> face_normals = {cube_face_pl[0], cube_face_pl[1], cube_face_pl[2]}; 
                    vector<double> face_c = cube_face_pl[3];
                    intersect_volume = cube_cone_intersection_volume_precomputed(zone_corners, edge_p0, edge_p1, face_normals, face_c, ray_dirs[i], cand_side_normals[i]); 
                }
                if(intersect_volume < 0.0){
                    cout << "negative volume iter" << endl; 
                    continue; //safeguard 
                }
                cout << "intersect volume " << izone << " " << pix_nums[i] << " "  << intersect_volume << endl; 
                double overlap_fraction = intersect_volume / max(this_zone_volume, 1.0e-30); 
                double net_light = quantity*overlap_fraction; 
                output[pix_nums[i]] += net_light; 
                counts[pix_nums[i]] += intersect_volume; 
            }
            if(!classify[0][i] && !classify[2][i])
                cout << "fully outside pixel " << pix_nums[i] << endl;
        }
    }
    ofstream outfile(filename); 
    for(int i = 0; i < 12*NSIDE*NSIDE; i++){
        outfile << output[i] << endl;
    }
    outfile.close();
    //fill vector of healpix maps for return
    vector<Healpix_Map<double>> out_maps(2); 
    out_maps.push_back(output); 
    out_maps.push_back(counts);
    cout << "done projector" << endl;
    return out_maps;
}

vector<double> normalize_3d(vector<double> input){
    double mag = sqrt(input[0]*input[0] + input[1]*input[1] + input[2]*input[2]); 
    vector<double> ans(3);  
    for(int i = 0; i < 3; i++)
        ans[i] = input[i] / mag; 
    return ans; 
}

vector<vector<double>> cube_face_planes(vector<vector<double>> corners, float tol){
    vector<vector<int>> cube_edges = {{0,1},{1,2},{2,3},{3,0},
                                          {4,5},{5,6},{6,7},{7,4},
                                          {0,7},{1,6},{2,5},{3,4}};
    vector<vector<int>> cube_faces = {{0,1,6,7},
                                      {2,3,4,5},
                                      {0,1,2,3},
                                      {4,5,6,7},
                                      {0,3,4,7},
                                      {1,2,5,6}};
    vector<vector<double>> res(4);
    for(int ind = 0; ind < 4; ind++)
        res[ind] = vector<double>(6); 
      
    //get center of cube 
    vector<double> center(3);
    for(int dim = 0; dim < 3; dim++){
         for(int corner = 0; corner < 8; corner++){
           center[dim] += corners[corner][dim];
        }
         center[dim] /= 8.0; 
    }
    vector<vector<double>> normals(6); 
    for(int ind = 0; ind < 6; ind++)
        normals[ind] = vector<double>(3); 
    vector<double> c(6); 
    for(int ind = 0; ind < 6; ind++)
        c[ind] = 0.0; 
    for(int ind = 0; ind < 6; ind++){
            vector<double> p0 = {corners[cube_faces[ind][0]][0], corners[cube_faces[ind][0]][1], corners[cube_faces[ind][0]][2]}; 
            vector<double> p1 = {corners[cube_faces[ind][1]][0], corners[cube_faces[ind][1]][1], corners[cube_faces[ind][1]][2]}; 
            vector<double> p2 = {corners[cube_faces[ind][2]][0], corners[cube_faces[ind][2]][1], corners[cube_faces[ind][2]][2]};
            vector<double> p1mp0(3); 
            vector<double> p2mp0(3);
            for(int k = 0; k < 3; k++){
                p1mp0[k] = p1[k] - p0[k]; 
                p2mp0[k] = p2[k] - p0[k]; 
            }
            vector<double> n(3); //result of cross product between p1mp0 and p2mp0
            //do the cross product between the vectors (there isn't a built-in cross product for vectors) 
            n[0] = p1mp0[1]*p2mp0[2] - p1mp0[2]*p2mp0[1]; 
            n[1] = p1mp0[2]*p2mp0[0] - p1mp0[0]*p2mp0[2]; 
            n[2] = p1mp0[0]*p2mp0[1] - p1mp0[1]*p2mp0[0]; 
            float nn = sqrt(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]); 
            //normalize result of cross product
            n[0] /= nn; 
            n[1] /= nn; 
            n[2] /= nn; 
            //cc = - n dot p0
            double cc = -1 * (n[0] * p0[0] + n[1] * p0[1] + n[2] * p0[2]);
            if(n[0]*center[0] + n[1]*center[1] + n[2]*center[2] + cc > 0.0){
                cc *= -1.0;
                n[0] *= -1.0; 
                n[1] *= -1.0; 
                n[2] *= -1.0; 
            }
            //fill res
            for(int ind_norm=0; ind_norm < 3; ind_norm++)
                res[ind_norm][ind] = n[ind_norm]; 
            res[3][ind] = cc; 
    }
   return res; 
}

vector<bool> points_in_cube(vector<vector<double>> pts, vector<vector<double>> face_normals, vector<double> face_c, float tol){
    vector<vector<bool>> vals(pts.size());
    vector<bool> collapsed_vals(pts.size()); 
    for(int i = 0; i < pts.size(); i++){
        vals[i] = vector<bool>(6);
        collapsed_vals[i] = true; 
    }  
    for(int i = 0; i < pts.size(); i++) 
        for(int j = 0; j < 6; j++){
            double res = 0.0;
            vals[i][j] = true; 
            for(int k = 0; k < 3; k++){
                res += pts[i][k] * face_normals[k][j];
            }
            res += face_c[j]; 
            if(res > tol)
                vals[i][j] = false;
            if(vals[i][j] == false)
                collapsed_vals[i] = false; //if one dimension outside of cube, then point is not in cube 
            }
        
    return collapsed_vals; 
}

bool origin_in_cube(vector<vector<double>> face_normals, vector<double> face_c, float tol){
    for(int i = 0; i < face_c.size(); i++)
        if(face_c[i] > tol)
            return false; 
    return true; 
}

vector<vector<double>> order_boundary_vecs_ccw(vector<vector<double>> ray_dirs, vector<double> center_vec){
    vector<vector<double>> res = ray_dirs; 
    double norm = sqrt(center_vec[0]*center_vec[0] + center_vec[1]*center_vec[1] + center_vec[2]*center_vec[2]); 
    if(norm < 1.e-30) 
        norm = 1.e-30; 
    for(int i = 0; i < 3; i++) 
        center_vec[i] /= norm; 
    vector<double> ref = {0.0, 0.0, 0.0}; 
    if(abs(center_vec[2]) < 0.9)
        ref = {0.0, 0.0, 1.0};
    else{ref = {1.0, 0.0, 0.0};}
    vector<double> e1 = {0.0, 0.0, 0.0}; 
    e1[0] = ref[1]*center_vec[2] - ref[2]*center_vec[1]; 
    e1[1] = ref[2]*center_vec[0] - ref[0]*center_vec[2]; 
    e1[2] = ref[0]*center_vec[1] - ref[1]*center_vec[0]; 
    double e1_norm = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]); 
    if(e1_norm <= 1e-30)
        e1_norm = 1.e-30; 
    for(int i = 0; i < 3; i++) 
        e1[i] /= e1_norm; 
    vector<double> e2 = {0.0, 0.0, 0.0};
    e2[0] = center_vec[1]*e1[2] - center_vec[2]*e1[1]; 
    e2[1] = center_vec[2]*e1[0] - center_vec[0]*e1[2]; 
    e2[2] = center_vec[0]*e1[1] - center_vec[1]*e1[0]; 
    double e2_norm = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]); 
    if(e2_norm <= 1.e-30)
        e2_norm = 1.e-30; 
    for(int i = 0; i < 3; i++) 
        e2[i] /= e2_norm;
    vector<double> x = {0.0, 0.0, 0.0, 0.0}; 
    vector<double> y = {0.0, 0.0, 0.0, 0.0}; 
    for(int i = 0; i < 4; i++) 
        for(int j = 0; j < 3; j++){
            x[i] += ray_dirs[i][j] * e1[j]; 
            y[i] += ray_dirs[i][j] * e2[j]; 
        }
    vector<double> ang = {0.0, 0.0, 0.0, 0.0}; 
    for(int i = 0; i < 4; i++) 
        ang[i] = atan2(y[i],x[i]); 
    vector<int> indices = {0,1,2,3}; 
    std::sort(indices.begin(), indices.end(), [&](int a, int b) {
        return ang[a] < ang[b]; //returns true if ang[a] < ang[b]? i think.  
    });
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 3; j++) 
            res[i][j] = ray_dirs[indices[i]][j]; 
    return res; 
}

vector<vector<double>> cone_side_normals(vector<vector<double>> ray_dirs){
    vector<vector<double>> res(8);
    for(int i = 0; i < 8; i++)
        res[i] = vector<double>(3); 
    for(int i = 0; i < 4; i++){ 
        double norm = sqrt(ray_dirs[i][0]*ray_dirs[i][0] + ray_dirs[i][1]*ray_dirs[i][1] + ray_dirs[i][2]*ray_dirs[i][2]);
        if(norm <= 1e-30)
            norm = 1.e-30; 
        for(int j = 0; j < 3; j++){ 
            ray_dirs[i][j] /= norm; 
            res[i+4][j] = ray_dirs[i][j]; 
        }
    }
    vector<double> interior = {0.0, 0.0, 0.0}; 
    for(int i = 0; i < 3; i++)
        interior[i] = (res[0][i] + res[1][i] + res[2][i] + res[3][i]) / 4.0; 
    double norm = sqrt(interior[0]*interior[0] + interior[1]*interior[1] + interior[2]*interior[2]); 
    if(norm <= 1.e-30)
        norm = 1.e-30; 
    for(int i = 0; i < 3; i++)
        interior[i] /= norm; 
    vector<vector<double>> u1 = ray_dirs; 
    for(int i = 0; i < 4; i++) 
        for(int j = 0; j < 3; j++){ 
            if(i == 3)
                u1[3][j] = ray_dirs[0][j];//cyclic permute back to 3  
            else{u1[i][j] = ray_dirs[i+1][j];} 
        }
    vector<vector<double>> normals = ray_dirs; 
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 3; j++){
            int index_1 = (j + 1) % 3; 
            int index_2 = (j + 2) % 3;
            normals[i][j] = ray_dirs[i][index_1]*u1[i][index_2] - ray_dirs[i][index_2]*u1[i][index_1];
        }
    }
    vector<double> nn = {0.0, 0.0, 0.0, 0.0}; 
    for(int i = 0; i < 4; i++){
        nn[i] = sqrt(normals[i][0]*normals[i][0] + normals[i][1]*normals[i][1] + normals[i][2]*normals[i][2]); 
    }
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 3; j++){
            normals[i][j] /= nn[i];
        }
    for(int i = 0; i < 4; i++){
        double cross_res = 0.0; 
        for(int j = 0; j < 3; j++){
            cross_res += normals[i][j] * interior[j]; 
        }
        if(cross_res < 0.0) 
           for(int j = 0; j < 3; j++)
               normals[i][j] *= -1.0; 
        for(int j = 0; j < 3; j++)
            res[i][j] = normals[i][j]; 
    }
    return res; 
}

vector<vector<vector<double>>> precompute_healpix_geometry(int nside){
    int npix = nside2npix(nside);
    Healpix_Base base(nside, RING, SET_NSIDE);  
    vector<vector<double>> pix_centers(npix); 
    vector<vector<vector<double>>> pix_rays(npix); 
    vector<vector<vector<double>>> pix_side_normals(npix); 
    vector<vector<vector<double>>> actual_res(npix); 
    for(int i = 0; i < npix; i++){
        pix_centers[i] = vector<double>(3); 
        pix_rays[i] = vector<vector<double>>(4); 
        pix_side_normals[i] = vector<vector<double>>(4); 
        actual_res[i] = vector<vector<double>>(9);
        for(int j = 0; j < 4; j++){
            pix_rays[i][j] = vector<double>(3); 
            pix_side_normals[i][j] = vector<double>(3);
        }
        for(int j = 0; j < 9; j++)
            actual_res[i][j] = vector<double>(3); 
        double center[3]; 
        pix2vec_ring(nside, i, center); 
        vector<double> center_vec = {center[0], center[1], center[2]};
        vector<vec3> rays;
        size_t step = 1; 
        base.boundaries(i, step, rays);
        vector<vector<double>> rays_input(4); 
        for(int j = 0; j < 4; j++){
            rays_input[j] = {rays[j].x, rays[j].y, rays[j].z}; 
        }
        rays_input = order_boundary_vecs_ccw(rays_input, center_vec);
        vector<vector<double>> res = cone_side_normals(rays_input); //first 4 elements are normals, last 4 are U 
        for(int j = 0; j < 9; j++){
            for(int k = 0; k < 3; k++){
                if(j == 0)
                    actual_res[i][j][k] = center[k]; 
                if(j >= 1 && j < 5){
                    actual_res[i][j][k] = rays_input[j-1][k];
                }
                if(j >= 5){ 
                    actual_res[i][j][k] = res[j-5][k]; 
                }
            }
        }
    }
    return actual_res; //format: 2nd index - 0 = centers, 1-4 = normals, 4-8 = U  
}

bool cube_fully_inside_cone(vector<vector<double>> corners, vector<vector<double>> face_normals, float tol){
    for(int i = 0; i < 8; i++)
        for(int j = 0; j < 4; j++){
            double res = 0.0; 
            for(int k = 0; k < 3; k++){
            res += corners[i][k]*face_normals[j][k]; 
            }
            if(res < -1.0 * tol)
                return false;
        }
    return true; 
}

bool cube_fully_outside_cone(vector<vector<double>> corners, vector<vector<double>> face_normals, float tol){
    vector<vector<double>> vals(8); 
    for(int i = 0; i < 8; i++)
        vals[i] = vector<double>(4); 
    for(int i = 0; i < 8; i++)
        for(int j = 0; j < 4; j++){
            double res = 0.0; 
            for(int k = 0; k < 3; k++){
            res += corners[i][k]*face_normals[j][k]; 
            }
            vals[i][j] = res; 
        }
    for(int j = 0; j < 4; j++) 
        for(int i = 0; i < 8; i++){ 
            if(vals[i][j] >= -1.0*tol) 
                break;
            if(i == 7) 
                return true; 
        }
    return false;
}

vector<vector<bool>> classify_pixels_for_zone(vector<vector<double>> zone_corners, vector<vector<vector<double>>> cand_side_normals, float tol){
    vector<vector<vector<double>>> vals(cand_side_normals.size()); 
    for(int i = 0; i < cand_side_normals.size(); i++){
        vals[i] = vector<vector<double>>(8); 
        for(int j = 0; j < 8; j++){ 
            vals[i][j] = vector<double>(4); 
            for(int k = 0; k < 4; k++){
                double res = 0.0; 
                for(int l = 0; l < 3; l++){
                    res += zone_corners[j][l] * cand_side_normals[i][k][l]; 
                }
                vals[i][j][k] = res; 
                }
        }
    }
    vector<vector<bool>> res(3); //0th axis: 0 = full_inside, 1 = full_outside, 2 = need_exact
    for(int i = 0; i < 3; i++)
        res[i] = vector<bool>(cand_side_normals.size()); 
    for(int i = 0; i < cand_side_normals.size(); i++){ 
        res[0][i] = true;//check for fully inside  
        res[1][i] = false; //check for fully outside
        for(int j = 0; j < 8; j++){
            for(int k = 0; k < 4; k++){
                if(vals[i][j][k] < -1.0*tol)
                    res[0][i] = false; 
            }
        }
        for(int k = 0; k < 4; k++){
            bool outside_check = true; 
            for(int j = 0; j < 8; j++){
                if(vals[i][j][k] >= -1.0*tol){
                    outside_check = false; 
                    break;
                }
            }
            if(outside_check){
                res[1][i] = true; 
                break; 
            }
        }
        //if neither fully inside or outside, needs exact intersection
        if(!(res[0][i] || res[1][i]))
            res[2][i] = true; 
        else{res[2][i] = false;}
    }
    return res; 
}

vector<bool> points_in_cone(vector<vector<double>> pts, vector<vector<double>> side_normals, float tol){
    vector<bool> ans(pts.size()); 
    vector<vector<double>> vals(pts.size()); 
    for(int i = 0; i < pts.size(); i++){ 
        ans[i] = true;
        vals[i] = vector<double>(4); 
        for(int j = 0; j < 4; j++){
            double res = 0.0; 
            for(int k = 0; k < 3; k++)
                res += pts[i][k] * side_normals[j][k];
            if(res < -1.0 * tol)
                ans[i] = false;

        }
    }
    return ans; 
}

pair<vector<vector<vector<double>>>, vector<vector<bool>>> segment_plane_intersections_batch(vector<vector<double>> p0, vector<vector<double>> p1, vector<vector<double>> plane_normals, float tol){
    vector<vector<double>> d(p1.size());
    vector<vector<double>> denom(d.size()); 
    vector<vector<double>> numer(d.size()); 
    vector<vector<vector<double>>> pts(d.size()); 
    for(int i = 0; i < d.size(); i++){
        denom[i] = vector<double>(plane_normals.size());
        numer[i] = vector<double>(plane_normals.size());
        d[i] = vector<double>(3); 
        pts[i] = vector<vector<double>>(plane_normals.size());
        for(int j = 0; j < 3; j++){
            d[i][j] = p1[i][j] - p0[i][j]; 
        }
        for(int j = 0; j < plane_normals.size(); j++)
            pts[i][j] = vector<double>(3); 
    }
    for(int i = 0; i < d.size(); i++){
        for(int j = 0; j < plane_normals.size(); j++){
            double denom_res = 0.0; 
            double numer_res = 0.0; 
            for(int k = 0; k < 3; k++){ 
                denom_res += d[i][k]*plane_normals[j][k]; 
                numer_res += -1.0 * p0[i][k] * plane_normals[j][k]; 
            }
            denom[i][j] = denom_res; 
            numer[i][j] = numer_res; 
        }
    }
    vector<vector<bool>> valid(denom.size());  
    vector<vector<double>> t(denom.size()); 
    for(int i = 0; i < denom.size(); i++){
        valid[i] = vector<bool>(denom[0].size());
        t[i] = vector<double>(denom[0].size());
        for(int j = 0; j < denom[0].size(); j++){
            t[i][j] = 0.0; 
            if(abs(denom[i][j]) > tol){
                valid[i][j] = true; 
                t[i][j] = numer[i][j] / denom[i][j];  
            }
            else{valid[i][j] = false;}
            if(t[i][j] < -1.0*tol || t[i][j] > (1.0 + tol))
                valid[i][j] = false; 
            //clip t_ij
            if(t[i][j] > 1.0)
                t[i][j] = 1.0; 
            if(t[i][j] < 0.0)
                t[i][j] = 0.0; 
            for(int k = 0; k < 3; k++){
                pts[i][j][k] = p0[i][k] + t[i][j]*d[i][k]; 
            }
        }
    }
   return {pts, valid};  
}

pair<vector<vector<vector<double>>>, vector<vector<bool>>> ray_plane_intersections_batch(vector<vector<double>> ray_dirs, vector<vector<double>> face_normals, vector<double> face_c, float tol){
    vector<vector<double>> denom(ray_dirs.size());
    vector<double> numer(face_normals[0].size());
    vector<vector<double>> t(ray_dirs.size());
    vector<vector<bool>> valid(ray_dirs.size()); 
    vector<vector<vector<double>>> pts(ray_dirs.size()); 
    for(int i = 0; i < ray_dirs.size(); i++){
        denom[i] = vector<double>(face_normals[0].size()); 
        t[i] = vector<double>(face_normals[0].size());
        valid[i] = vector<bool>(face_normals[0].size()); 
        pts[i] = vector<vector<double>>(face_normals[0].size()); 
        for(int j = 0; j < face_normals[0].size(); j++){
            pts[i][j] = vector<double>(3); 
            double res = 0.0; 
            for(int k = 0; k < 3; k++)
                res += ray_dirs[i][k]*face_normals[k][j]; 
            denom[i][j] = res; 
            numer[j] = -1.0 * face_c[j];
            t[i][j] = 0.0;
            if(abs(res) > tol)
                valid[i][j] = true;
            else{valid[i][j] = false;}
            if(valid[i][j])
                t[i][j] = numer[j]/denom[i][j];
            if(valid[i][j] && t[i][j] >= -1*tol)
                valid[i][j] = true; 
            else{valid[i][j] = false;}
            t[i][j] = max(t[i][j], 0.0); 
            for(int k = 0; k < 3; k++)
                pts[i][j][k] = ray_dirs[i][k] * t[i][j]; 
        }
    }
    return {pts, valid}; 
}

vector<vector<double>> unique_points(vector<vector<double>> pts, float tol){
    if(pts.empty()){
        return pts; 
    }
    if(pts.size() == 0){
        return pts; 
    }
    vector<vector<double>> unique; 
    vector<vector<long>> key(pts.size()); 
    for(int i = 0; i < pts.size(); i++){
        key[i] = vector<long>(pts[0].size()); 
        for(int j = 0; j < pts[0].size(); j++){
            key[i][j] = lround(pts[i][j]/tol); 
        }
    }
    vector<vector<long>> unique_rows; 
    vector<int> return_index; 

    map<vector<long>, int> seen_rows; 

    for(int i = 0; i < key.size(); i++){ 
        if(seen_rows.count(key[i]) == 0){
            seen_rows[key[i]] = i; 
            unique_rows.push_back(key[i]); 
            return_index.push_back(i); 
        }
    }


    for(int i = 0; i < return_index.size(); i++) 
           unique.push_back({pts[return_index[i]][0], pts[return_index[i]][1], pts[return_index[i]][2]});

    return unique; 
    
}

vector<vector<double>> cube_cone_intersection_vertices_precomputed(vector<vector<double>> corners, vector<vector<double>> edge_p0, vector<vector<double>> edge_p1, vector<vector<double>> face_normals, vector<double> face_c, vector<vector<double>> ray_dirs, vector<vector<double>> side_normals, float tol){
   vector<vector<double>> pts_list; 
   //1) cube corners inside cone 
   vector<bool> mask_corners = points_in_cone(corners, side_normals, tol = tol);
   for(int i = 0; i < mask_corners.size(); i++) 
       if(mask_corners[i] == true)
           pts_list.push_back(corners[i]); 
   //2) cone apex inside cube
   if(origin_in_cube(face_normals, face_c, tol=tol))
       pts_list.push_back({0.0, 0.0, 0.0});
   //3) cube edge/cone side plane intersections
   pair<vector<vector<vector<double>>>, vector<vector<bool>>> seg_res = segment_plane_intersections_batch(edge_p0, edge_p1, side_normals, tol=tol);  
   vector<vector<double>> cand; 
   for(int i = 0; i < seg_res.second.size(); i++) 
       for(int j = 0; j < seg_res.second[0].size(); j++){
           if(seg_res.second[i][j] == true){
               cand.push_back(seg_res.first[i][j]);
           }   
        }
    vector<bool> pts_cone = points_in_cone(cand, side_normals, tol = tol);
    vector<bool> pts_cube = points_in_cube(cand, face_normals, face_c, tol=tol);
    for(int k = 0; k < cand.size(); k++){
        if(pts_cone[k] && pts_cube[k]) 
            pts_list.push_back(cand[k]); 
    }
   
   
   //4) cone edge ray/cube face intersections
   pair<vector<vector<vector<double>>>, vector<vector<bool>>> ray_res = ray_plane_intersections_batch(ray_dirs, face_normals, face_c, tol=tol); 
   cand.clear();
   for(int i = 0; i < ray_res.second.size(); i++) 
       for(int j = 0; j < ray_res.second[0].size(); j++){ 
            if(ray_res.second[i][j] == true){
                cand.push_back(ray_res.first[i][j]); 
            }
        }
    pts_cone = points_in_cone(cand, side_normals, tol = tol); 
    pts_cube = points_in_cube(cand, face_normals, face_c, tol = tol); 
    for(int k = 0; k < cand.size(); k++){
        if(pts_cone[k] && pts_cube[k])
            pts_list.push_back(cand[k]); 
    }

   return unique_points(pts_list, tol = tol); 
}

double cube_cone_intersection_volume_precomputed(vector<vector<double>> corners, vector<vector<double>> edge_p0, vector<vector<double>> edge_p1, vector<vector<double>> face_normals, vector<double> face_c, vector<vector<double>> ray_dirs, vector<vector<double>> side_normals, float tol){
    vector<vector<double>> verts = cube_cone_intersection_vertices_precomputed(corners, edge_p0, edge_p1, face_normals, face_c, ray_dirs, side_normals, tol=1e-9);
    if(verts.empty())
        return 0.0;
    if(verts.size() < 4)
        return 0.0; 
    vector<vector<double>> centered = verts; 
    vector<double> mean(3); 
    for(int i = 0; i < 3; i++) 
        mean[i] = 0.0; 
    for(int i = 0; i < verts.size(); i++)
        for(int j = 0; j < 3; j++)
            mean[j] += verts[i][j];
    for(int i = 0; i < 3; i++) 
        mean[i] /= verts.size(); 
    for(int i = 0; i < verts.size(); i++) 
        for(int j = 0; j < 3; j++) 
            centered[i][j] -= mean[j]; 
    Eigen::MatrixXd mat(centered.size(), centered[0].size()); 
    for(int i = 0; i < centered.size(); i++) 
        for(int j = 0; j < centered[0].size(); j++)
            mat(i,j) = centered[i][j]; 
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinU | Eigen::ComputeThinV); 
    int rank = 0; 
    for(int i = 0; i < svd.singularValues().size(); i++)
        if(svd.singularValues()(2) > 1e-8 * svd.singularValues()(0))
            rank++;
    if(rank < 3)
        return 0.0;
    if(svd.singularValues()(2) < 1e-6 * svd.singularValues()(0)) 
        return 0.0; 
    cout << "rank " << rank << endl; 
    cout << "points" << endl;
    for(int i = 0; i < verts.size(); i++) 
        cout << centered[i][0] << " " << centered[i][1] << " " << centered[i][2] << endl;
    //qhull wants a flat vector
    vector<double> flat_points(3 * centered.size()); 
    for(int i = 0; i < centered.size(); i++){ 
        for(int j = 0; j < 3; j++){
            flat_points[3*i + j] = centered[i][j]; 
        }
    }
    //cout << "about to run qhull" << endl; 
    //make convex hull 
    orgQhull::Qhull qhull("", 3, flat_points.size()/3, flat_points.data(), "Qt QJ"); 
    return qhull.volume(); 
}

