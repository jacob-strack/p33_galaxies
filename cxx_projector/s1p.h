#include<geos/geom/LinearRing.h>
#include<geos/geom/GeometryFactory.h>
#include<geos/geom/Polygon.h>
#include<geos/algorithm/ConvexHull.h>
#include<gmp.h>
#include"proj.h"
#include<algorithm>
#include<healpix_base.h>
#include<healpix_map.h>
#include<healpix_map_fitsio.h>
#include<alm.h>
#include<deque>
#include<pointing.h>
//using namespace std; 

vector<Healpix_Map<double>> project(vector<double> cube, vector<vector<double>> xyz, vector<vector<double>> dxyz, float center[3], float projax[3], int NSIDE,int exclude = 1); 

vector<Healpix_Map<double>> project(vector<double> cube, vector<vector<double>> xyz, vector<vector<double>> dxyz, float center[3], float projax[3], int NSIDE, int exclude){
    geos::geom::GeometryFactory *factory;
    bool verbose = false;
    bool use_geos = true; 
    for(int dim = 0; dim < 3; dim++)
        for(int i = 0; i < xyz[0].size(); i++)
            xyz[dim][i] = xyz[dim][i] - center[dim]; 
    vector<float> r(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size(); i++)
        r[i] = xyz[0][i]*xyz[0][i] + xyz[1][i]*xyz[1][i] + xyz[2][i]*xyz[2][i];  
    vector<int>nxyz(r.size()); 
    for(int i = 0; i < xyz[0].size(); i++){
        r[i] = sqrt(r[i]);
        nxyz[i] = r[i] / dxyz[0][i]; 
    }
    vector<bool> mask(r.size()); //mask out close points 
    int Nz = 0; 
    for(int i = 0; i < r.size(); i++){
        if(nxyz[i] > exclude){
            mask[i] = true;
            Nz++; 
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
    float no_center[3] = {0.0, 0.0, 0.0};//duct tape for now where i dont want the fucntion
                                         //to shift anymore bc i already did it earlier
    vector<vector<double>> xyz_p = rotate(xyz, projax, no_center); 

    //shifter object
    vector<vector<vector<double>>> corners(8); 
    float shifter[8][3] = {{-0.5, -0.5, -0.5}, {-0.5, -0.5, 0.5}, {0.5,-0.5,0.5},{0.5,-0.5,-0.5},{0.5,0.5,-0.5},{0.5,0.5,0.5},{-0.5,0.5,0.5},{-0.5,0.5,-0.5}};
    for(int i = 0; i < 8; i++){
        corners[i] = vector<vector<double>>(3); 
        for(int j = 0; j < 3; j++){
            corners[i][j] = vector<double>(xyz[0].size());
            for(int k = 0; k < xyz[0].size(); k++)
                corners[i][j][k] = xyz[j][k] + dxyz[j][k]*shifter[i][j]; 
        }
    }
    for(int i = 0; i < 8; i++)
        corners[i] = rotate(corners[i], projax, no_center);
    //accumulation things 
    vector<float> r_sq(xyz[0].size());
    vector<float> zone_volume(xyz[0].size()); 
    vector<float> zone_emission(xyz[0].size());
    for(int i = 0; i < xyz[0].size(); i++){
        r_sq[i] = xyz_p[0][i]*xyz_p[0][i] + xyz_p[1][i]*xyz_p[1][i] + xyz_p[2][i]*xyz_p[2][i]; 
        zone_volume[i] = dxyz[0][i]*dxyz[1][i]*dxyz[2][i]; 
        zone_emission[i] = cube[i]/r_sq[i]*zone_volume[i];
    }

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
   // for(int i = 0; i < xyz[0].size(); i++){ 
     //   cout << "corners " << xyz_p[0][i] << " " << xyz_p[1][i] << " " << xyz_p[2][i] << i << endl; 
      //  for(int j = 0; j < 8; j++)
       //     cout << corner_angles[j][0][i] << " " << corner_angles[j][1][i] << endl;
//}
    //take out periodic wrap
    vector<float> phi_min(xyz[0].size()); 
    vector<float> phi_max(xyz[0].size());
    vector<float> theta_min(xyz[0].size()); 
    vector<float> theta_max(xyz[0].size()); 
    for(int i = 0; i < 8; i++)
        for(int j = 0; j < xyz[0].size(); j++){
           if(i==0){
                phi_min[j] = corner_angles[0][0][j]; 
                phi_max[j] = corner_angles[0][0][j]; 
                theta_min[j] = corner_angles[0][1][j]; 
                theta_max[j] = corner_angles[0][1][j];
           }
           else{
                if(corner_angles[i][0][j] > phi_max[j])
                    phi_max[j] = corner_angles[i][0][j];
                if(corner_angles[i][0][j] < phi_min[j])
                    phi_min[j] = corner_angles[i][0][j];
                if(corner_angles[i][1][j] > theta_max[j])
                    theta_max[j] = corner_angles[i][1][j];
                if(corner_angles[i][1][j] < theta_min[j])
                    theta_min[j] = corner_angles[i][1][j]; 
            }
        }
    vector<float> dphi(xyz[0].size()); 
    vector<float> seam_shift(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size();i++){
        dphi[i] = phi_max[i] - phi_min[i];
        if(dphi[i] > M_PI/2)
            seam_shift[i] = 2*M_PI; 
        else
            seam_shift[i] = 0; 
    }
    //need to add shift only if the corner phi is < pi/2
    for(int i = 0; i < 8; i++)
        for(int j = 0; j < xyz[0].size(); j++){
                if(corner_angles[i][0][j]<M_PI/2){
                    corner_angles[i][0][j] += seam_shift[j]; 
                }
            }

    //center of each zone in theta and phi 
    vector<double>phi_cen(xyz[0].size()); 
    vector<double>theta_cen(xyz[0].size()); 

    for(int i = 0; i < xyz[0].size(); i++){
        float mean_phi = 0;
        float mean_theta = 0;
        for(int j = 0; j < 8; j++){
            mean_phi += corner_angles[j][0][i]; 
            mean_theta += corner_angles[j][1][i]; 
        }
        phi_cen[i] = mean_phi/8.0; 
        theta_cen[i] = mean_theta/8.0;
    }

    //angular size and cut out circle with radius of max angular distance 
    vector<vector<float>> distance(8);
    for(int i = 0; i < 8; i++){
            distance[i] = vector<float>(xyz[0].size()); 
            for(int k = 0; k < xyz[0].size(); k++){
                distance[i][k] = sqrt((corner_angles[i][0][k] - phi_cen[k])*(corner_angles[i][0][k] - phi_cen[k]) + (corner_angles[i][1][k] - theta_cen[k])*(corner_angles[i][1][k] - theta_cen[k]));  
            }
    }

    vector<float> circle_radius(xyz[0].size()); 

    for(int i = 0; i < xyz[0].size(); i++){
        float max_distance = 0;
        for(int j = 0; j < 8; j++){
            if(j==0)
                max_distance = distance[j][i]; 
            else{
                if(distance[j][i] > max_distance)
                    max_distance = distance[j][i]; 
            }
        }
        circle_radius[i] = max_distance; 
    }
    
    //get pole zones
    vector<vector<float>> max_pts(8); 
    vector<vector<float>> min_pts(8);
    vector<bool> theta_pole_check(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size(); i++)
        theta_pole_check[i] = false;
    for(int i = 0; i < 8; i++){
        max_pts[i] = vector<float>(xyz[0].size()); 
        min_pts[i] = vector<float>(xyz[0].size()); 
    }
    for(int i = 0; i < xyz[0].size(); i++)
        for(int j = 0; j < 3; j++) 
           for(int k = 0; k < 8; k++){
                if(k==0){
                    max_pts[j][i] = corners[0][j][i]; 
                    min_pts[j][i] = corners[0][j][i]; 
                }
                if(corners[k][j][i] > max_pts[j][i])
                    max_pts[j][i] = corners[k][j][i];
                if(corners[k][j][i] < min_pts[j][i])
                    min_pts[j][i] = corners[k][j][i];
                if(abs(cos(corner_angles[k][1][i])) > 2./3.){
                    theta_pole_check[i] = true;
                }
            }
    //create array of map length to keep track of pole zones
    vector<bool> is_pole(xyz[0].size());
    vector<bool> south_pole(xyz[0].size()); 
    for(int i = 0; i < xyz[0].size(); i++){
        if(max_pts[0][i] < 0)
            is_pole[i] = false;
        else{is_pole[i] = false;}
        if(theta_pole_check[i] == true)
           is_pole[i] = true; 
        if(min_pts[0][i] >=  0)
            south_pole[i] = true; 
        if(max_pts[0][i] <= 0)
            south_pole[i] = false;
    }


   //ZONE LOOOOOOOOOOOOOOOOOOOOOP 
   int debug_ind = 0;
   for(int izone = 0; izone < xyz[0].size(); izone++){
       cout << izone << endl;
       if(mask[izone]==false)
            continue; 
        vector<float> all_phi(8);
        vector<float> all_theta(8); 
        for(int i = 0; i < 8; i++){
            all_phi[i] = corner_angles[i][0][izone]; 
            all_theta[i] = corner_angles[i][1][izone]; 
        }
        
        //switch coordinate systems at the poles 
        vector<float> all_horz(8); 
        vector<float> all_vert(8);
        if(is_pole[izone]){
            vector<float> this_r(8); 
            if(south_pole[izone]){
                for(int i = 0; i < 8; i++)
                    this_r[i] = M_PI - all_theta[i];
            }
            else{
                for(int i = 0; i < 8; i++)
                    this_r[i] = all_theta[i]; 
            }
            for(int i = 0; i < 8; i++){
                all_horz[i] = this_r[i]*sin(all_phi[i]); 
                all_vert[i] = this_r[i]*cos(all_phi[i]);
            }
        }
        else{
            for(int i = 0; i < 8; i++){
                all_horz[i] = all_phi[i]; 
                all_vert[i] = all_theta[i];
            }   
        }
        vector<double> in_pts; 
        for(int i = 0; i < 8; i++){
            in_pts.push_back(all_horz[i]); 
            in_pts.push_back(all_vert[i]);
        }
        //make convex hull of zone corners, get edges
        vector<double> phi_edge; 
        vector<double> theta_edge;
        if(1){
        geos::geom::Coordinate hull_coords; 
        vector<unique_ptr<geos::geom::Geometry>> hull_geo;
        for(int i = 0; i < 8; i++){
            hull_coords = geos::geom::Coordinate(all_horz[i], all_vert[i]);
            hull_geo.push_back(factory->createPoint(hull_coords));
        }
        unique_ptr<geos::geom::GeometryCollection> collection = factory->createGeometryCollection(move(hull_geo));
        geos::algorithm::ConvexHull chull(collection.get()); 
        unique_ptr<geos::geom::Geometry> res = chull.getConvexHull();
        unique_ptr<geos::geom::CoordinateSequence> res_pts = res->getCoordinates();
        for(int i = 0; i < res_pts->getSize(); i++){
            phi_edge.push_back(res_pts->getAt(i).x); 
            theta_edge.push_back(res_pts->getAt(i).y);
        }
   }
        unique_ptr<geos::geom::Polygon> geos_hull_poly;
        if(use_geos){
            unique_ptr<geos::geom::CoordinateSequence> crds = make_unique<geos::geom::CoordinateSequence>(); 
            //sort points clockwise 
            double phi_center = 0; 
            double theta_center = 0; 
            for(int i = 0; i < phi_edge.size(); i++){
                phi_center += phi_edge[i]; 
                theta_center += theta_edge[i]; 
            }
            theta_center/=phi_edge.size(); 
            phi_center/=phi_edge.size(); 
            vector<float>sort_angle; 
            for(int i = 0; i < phi_edge.size(); i++){
                sort_angle.push_back(atan2(theta_edge[i] - theta_center, phi_edge[i] - phi_center));
            }
            vector<float>sorted_angle = sort_angle;
            sort(sorted_angle.begin(), sorted_angle.end());
            int first_ind;
            for(int i = 0; i < sorted_angle.size(); i++){
                int index; 
                for(int j = 0; j < sort_angle.size(); j++){
                   if(sort_angle[j] == sorted_angle[i])
                       index = j;
                }
                if(i==0){
                    crds->add(phi_edge[index], theta_edge[index]);
                    first_ind = index;
                }
                else{
                    if(sorted_angle[i] == sorted_angle[i - 1])
                        continue; 
                    crds->add(phi_edge[index], theta_edge[index]); 
                }
            }
            crds->add(phi_edge[first_ind], theta_edge[first_ind]); 
            unique_ptr<geos::geom::LinearRing> ring = factory->createLinearRing(move(crds));
            geos_hull_poly = factory->createPolygon(move(ring));
        }
        double zone_min_theta;
        double zone_max_theta;
        for(int i = 0; i < 8; i++){
            if(i == 0){
                zone_min_theta = all_theta[i]; 
                zone_max_theta = all_theta[i];
            }
            else{
                if(all_theta[i] < zone_min_theta)
                    zone_min_theta = all_theta[i]; 
                if(all_theta[i] > zone_max_theta)
                    zone_max_theta = all_theta[i];
            }
        }
        //make convex hull exterior polygon
        pointing this_pointing(theta_cen[izone], phi_cen[izone]);
        vector<int> pix_nums; 
        base.query_disc_inclusive(this_pointing, circle_radius[izone], pix_nums);
        //loop over found pixels in hp map within radius
        for(int ipix = 0; ipix < pix_nums.size(); ipix++){
           vector<vec3> xyz_pointing; 
           base.boundaries(pix_nums[ipix],1,xyz_pointing);
           //things are returned in cartesian, convert to a pointing 
           vector<pointing> pix_angs(4);
           for(int i = 0; i < 4; i++){
                pix_angs[i] = pointing(xyz_pointing[i]);
            }
            float min_theta; 
            float max_theta;
            vector<double> phi(4); 
            vector<double> theta(4);
            for(int i = 0; i < 4; i++){
                phi[i] = pix_angs[i].phi; 
                theta[i] =  pix_angs[i].theta;
                if(i == 0){
                    min_theta = theta[i]; 
                    max_theta = theta[i];
                }
                else{
                    if(theta[i] < min_theta)
                        min_theta = theta[i]; 
                    if(theta[i] > max_theta)
                        max_theta = theta[i];
                }
            }
    
            if(min_theta > zone_max_theta)
                continue;
            if(max_theta < zone_min_theta)
                continue;
            //find max/min phi for the zone
            float zone_min_phi; 
            float zone_max_phi;
            for(int z = 0; z < all_phi.size(); z++){
                if(z==0){
                    zone_min_phi = all_phi[0]; 
                    zone_max_phi = all_phi[0]; 
                }
                else{
                    if(all_phi[z] > zone_max_phi)
                        zone_max_phi = all_phi[z]; 
                    if(all_phi[z] < zone_min_phi)
                        zone_min_phi = all_phi[z];
                }
            }
           //check for concavity
           //step around zone, checking area
           float zxprod[4]; 
           for(int i = 0; i < 4; i++){
                int ind1 = i; 
                int ind2 = i + 1; 
                int ind3 = i + 2; 
                if(ind2 > 3)
                    ind2 -= 4;
                if(ind3 > 3)
                    ind3 -= 4;
                vector<float> x(4);
                vector<float> y(4);
                for(int i = 0; i < 4; i++){
                    x[i] = theta[i]; 
                    y[i] = phi[i];
                }
                float dx1 = x[ind2] - x[ind1]; 
                float dy1 = y[ind2] - y[ind1]; 
                float dx2 = x[ind3] - x[ind2]; 
                float dy2 = y[ind3] - y[ind2]; 
                zxprod[i] = dx1*dy2 - dx2*dy1;
            }
           float abs_zxprod_sum = 0;
           float zxprod_abs_sum = 0; 
           int nneg = 0; 
           for(int i = 0; i < 4; i++){
               zxprod_abs_sum += abs(zxprod[i]); 
               abs_zxprod_sum += zxprod[i];
               if(zxprod[i] < 0)
                   nneg++; 
            }
           abs_zxprod_sum = abs(abs_zxprod_sum); 
           bool shifted = false;
           if(zxprod_abs_sum - abs_zxprod_sum > 0){
               //shift appropriate points if the coordinates are not in fake cartesian
               //I think there's something missing in the logic with nneg=1
               if(nneg == 3){
                    int ind_to_shift; 
                    for(int i = 0; i < 4; i++)
                        if(zxprod[i] > 0)
                            ind_to_shift = i - 1; 
                        if(ind_to_shift < 0)
                            ind_to_shift = 3;
                            cout << "shifting nneg 3 " << ind_to_shift << endl;
                            shifted = true;
                        if(pix_angs[ind_to_shift].phi < M_PI/2)
                            pix_angs[ind_to_shift].phi += 2*M_PI;
                        else{pix_angs[ind_to_shift].phi -= 2*M_PI;} 
                }
               if(nneg == 1){
                    int ind_to_shift; 
                    for(int i = 0; i < 4; i++)
                        if(zxprod[i] < 0)
                            ind_to_shift = i; 
                        if(ind_to_shift < 0)
                            ind_to_shift += 4;
                        shifted = true;
                        cout << "shifting nneg 1 " << ind_to_shift << endl;
                        cout << "before shift: " << endl; 
                        for(int i = 0; i < 4; i++) 
                            cout << pix_angs[i].phi << " " << pix_angs[i].theta << endl;
                        if(pix_angs[ind_to_shift].phi < M_PI/2)
                            pix_angs[ind_to_shift].phi += 2*M_PI;
                        else{pix_angs[ind_to_shift].phi -= 2*M_PI;}
                        cout << "after shift: " << endl; 
                        for(int i = 0; i < 4; i++) 
                            cout << pix_angs[i].phi << " " << pix_angs[i].theta << endl;
                }
           }

            //get min/max of pixels post correction

            float max_phi, min_phi; 
           
            for(int i = 0; i < 4; i++){
                phi[i] = pix_angs[i].phi; 
                theta[i] =  pix_angs[i].theta;
                if(i == 0){
                    min_phi = phi[i]; 
                    max_phi = phi[i];
                }
                else{
                    if(phi[i] < min_phi)
                        min_phi = phi[i]; 
                    if(phi[i] > max_phi)
                        max_phi = phi[i]; 
                }
            }
            //undo the periodic wrap for healpix pixels
            if(zone_min_phi > max_phi){
                for(int i = 0; i < 4; i++)
                    phi[i] += 2*M_PI; 
            }
            if(zone_max_phi < min_phi){
                for(int i = 0; i < 4; i++)
                    phi[i] -= 2*M_PI;
            } 
            //change to fake cartesian if we are in a pole zone
            for(int i = 0; i < 4; i++){
                if(is_pole[izone]){
                    float this_r = 0; 
                    if(south_pole[izone]){
                        this_r = M_PI - theta[i]; 
                    }
                    else{this_r = theta[i]; 
                    }
                theta[i] = this_r*cos(phi[i]); 
                phi[i] = this_r*sin(phi[i]);
                //
            }
            }
            double phi_center = 0; 
            double theta_center = 0; 
            unique_ptr<geos::geom::Polygon> geos_pix_poly;
            if(use_geos){
                unique_ptr<geos::geom::CoordinateSequence> pix_crds = make_unique<geos::geom::CoordinateSequence>(); 
                //sort points clockwise 
                for(int i = 0; i < phi.size(); i++){
                    phi_center += phi[i]; 
                    theta_center += theta[i]; 
                }
                theta_center/=phi.size(); 
                phi_center/=phi.size(); 
                vector<float>sort_angle; 
                for(int i = 0; i < phi.size(); i++){
                    sort_angle.push_back(atan2(theta[i] - theta_center, phi[i] - phi_center));
                }
                vector<float>sorted_angle = sort_angle;
                sort(sorted_angle.begin(), sorted_angle.end());
                int first_ind;
                for(int i = 0; i < sorted_angle.size(); i++){
                    int index; 
                    for(int j = 0; j < sort_angle.size(); j++){
                       if(sort_angle[j] == sorted_angle[i])
                           index = j;
                    }
                    if(i==0){
                        pix_crds->add(phi[index], theta[index]);
                        first_ind = index;
                    }
                    else{
                        if(sorted_angle[i] == sorted_angle[i - 1])
                            continue; 
                        pix_crds->add(phi[index], theta[index]); 
                    }
                }
                pix_crds->add(phi[first_ind], theta[first_ind]); 
                unique_ptr<geos::geom::LinearRing> pix_ring = factory->createLinearRing(move(pix_crds));
                geos_pix_poly = factory->createPolygon(move(pix_ring));
            }
            //intersect the polygons and store in intersection
            if(use_geos){
                unique_ptr<geos::geom::Geometry> intersection = unique_ptr<geos::geom::Geometry>(geos_hull_poly->intersection(geos_pix_poly.get())); 
                double myarea = intersection->getArea();
                if(1){
                    counts[pix_nums[ipix]] += 1;
                    output[pix_nums[ipix]] += myarea*cube[izone]/(r_sq[izone]*(geos_hull_poly->getArea())); 
                }
                if(verbose){
                    cout << "pole " << is_pole[izone] << endl;
                    cout << "south pole " << south_pole[izone] << endl;
                    cout << "pix num " << pix_nums[ipix] << endl;
                    cout << "zone cen" << endl; 
                    cout << phi_cen[izone] << " " << theta_cen[izone] << endl;
                    cout << "healpix angles" << endl; 
                    for(int i = 0; i < phi.size(); i++)
                        cout << phi[i] << " " << theta[i] << endl; 
                    cout << "zone corners" << endl;
                    for(int i = 0; i < phi_edge.size(); i++){
                        cout << phi_edge[i] << " " << theta_edge[i] << " " << min_pts[0][izone] << " " << max_pts[0][izone] << endl; 
                    }
                    cout << "all zone points" << endl; 
                    for(int i = 0; i < 8; i++)
                        cout << corner_angles[i][0][izone] << " " << corner_angles[i][1][izone] << endl; 
                    cout << "area: " << myarea << endl;
                    cout << output[pix_nums[ipix]] << endl;
                }
            }
        }
   }
    ofstream outfile("plot_array.txt"); 
    for(int i = 0; i < 12*NSIDE*NSIDE; i++){
        if(output[i] > 1.4)
            cout << i << " " << output[i] << endl;
        outfile << output[i] << endl;
    }
    outfile.close();
    //fill vector of healpix maps for return
    vector<Healpix_Map<double>> out_maps(2); 
    out_maps.push_back(output); 
    out_maps.push_back(counts);
    return out_maps;
}
