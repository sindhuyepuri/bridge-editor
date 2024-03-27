#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include <cassert>
#include "gurobi_c++.h"
#include "gurobi_c.h"

using namespace std;

pair<int, int> make_edge(int x, int y) {
    if (x < y) {
        return pair<int, int>(x, y);
    } else {
        return pair<int, int>(y, x);
    }
}

void add_neighbor(map<int, set<int>>& m, int key, int neighbor) {
    if (m.find(key) == m.end()) {
        set<int> s;
        s.insert(neighbor);
        m[key] = s;
    } else {
        m[key].insert(neighbor);
    }
}

float get_random() {
    static std::default_random_engine e;
    static std::uniform_real_distribution<> dist(1, 11);
    return dist(e);
}

int main(int argc, char **argv) {
    polyscope::init();
    vector<glm::vec3> points;

    /*
        File name for triangle edges and (x, y) points of vertices
    */
    string file_name;
    getline(cin, file_name);
    
    ifstream file("../data/" + file_name + "-xy.txt");
    int num_nodes;
    file >> num_nodes;
    int ignore;
    file >> ignore;
    file >> ignore;
    file >> ignore;

    /*
        Mapping of vertex index to (x, y) coordinate
    */
    map<int, float> x_i;
    map<int, float> y_i;

    for (int i = 0; i < num_nodes; i++) {
        int p_idx;
        file >> p_idx;
        float x, y;
        file >> x;
        file >> y;
        file >> ignore;
        x_i.insert(pair<int, float>(p_idx, x));
        y_i.insert(pair<int, float>(p_idx, y));
    }

    file.close();

    file.open("../data/" + file_name + "-triangles.txt");
    
    cout << file_name << endl;

    int num_triangles;
    file >> num_triangles;
    file >> ignore;
    file >> ignore;
    
    /*
        Edges contains only unique edges, for instance (1, 3) == (3, 1)
        So only contains information about vertex connectivity
    */
    set<pair<int, int>> edges;
    /*
        Maps vertex index to set of neighbor vertex indices
    */
    map<int, set<int>> neighbors; 
    for (int i = 0; i < num_triangles; i++) {
        file >> ignore;
        int p1, p2, p3;

        file >> p1;
        file >> p2;
        file >> p3;

        add_neighbor(neighbors, p1, p2);
        add_neighbor(neighbors, p2, p1);
        add_neighbor(neighbors, p2, p3);
        add_neighbor(neighbors, p3, p2);
        add_neighbor(neighbors, p1, p3);
        add_neighbor(neighbors, p3, p1);

        edges.insert(make_edge(p1, p2));
        edges.insert(make_edge(p1, p3));
        edges.insert(make_edge(p3, p2));
    }

    /*
        Pinned vertex indices represents vertices that are rooted in ground
        -> should handle infinite load
        -> used for visualization purposes, since no added constraint
    */
    set<int> pinned_indices;
    for (int i = 0; i < num_nodes; i++) {
        if (x_i[i + 1] == 210 || x_i[i + 1] == -210 || y_i[i + 1] == 210 || y_i[i + 1] == -210) {
            pinned_indices.insert(i + 1);
        }
    }

    try {
        double total_load = 1000.0;
        double load = total_load / num_nodes;
        
        set<pair<int, int>>::iterator itr;
        map<pair<int, int>, int> scalar_idx;
        int i = 0;
        for (itr = edges.begin(); itr != edges.end(); itr++) {
            auto edge = *itr;
            scalar_idx[edge] = i++;
        }

        // Start with random edge scalars
        float subs_scalars[edges.size()];
        for (int i = 0; i < edges.size(); i++) {
            subs_scalars[i] = get_random();
            cout << "old scalar " << subs_scalars[i] << endl;
        }
        float subs_z[num_nodes];

        int iters = 0;
        bool toggling = true; // while true situation
        bool toggle_solve_z = true;
        bool toggle_solve_scalar = false;

        float prev_z_abs_diff = 1e8;
        float prev_scalar_abs_diff = 1e8;
        while (toggling) {
            iters++;
            float z_abs_diff = 0;
            float scalar_abs_diff = 0;
            if (toggle_solve_z) {
                cout << "solving for z" << endl;
                GRBEnv env = GRBEnv();
                GRBModel model = GRBModel(env);

                GRBVar z_i[num_nodes + 1];
                for (int i = 0; i < num_nodes; i++) {
                    string z = "z_" + to_string(i + 1);
                    z_i[i + 1] = model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, z);
                }
                for (int i = 0; i < num_nodes; i++) {
                    set<int> neighbors_i = neighbors[i + 1];
                    set<int>::iterator itr;
                    GRBLinExpr x_comp;
                    GRBLinExpr y_comp;
                    GRBLinExpr z_comp;

                    float e_ix = x_i[i + 1];
                    float e_iy = y_i[i + 1];
                    GRBVar e_iz = z_i[i + 1];

                    // pinned vertices (replace with better method)
                    auto pos = pinned_indices.find(i + 1);
                    if (pos == pinned_indices.end()) {
                        for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                            int neighbor = *itr;
                            pair<int, int> edge = make_edge(i + 1, neighbor);
                            float e_jx = x_i[neighbor];
                            float e_jy = y_i[neighbor];
                            GRBVar e_jz = z_i[neighbor];
                            float s = subs_scalars[scalar_idx[edge]]; // use substituded scalar instead of var

                            x_comp += (e_ix - e_jx) * s;
                            y_comp += (e_iy - e_jy) * s;
                            z_comp += (e_iz - e_jz) * s;
                        }
                        string x_constr = to_string(i + 1)+"_xcmpnt";
                        string y_constr = to_string(i + 1)+"_ycmpnt";
                        string z_constr = to_string(i + 1)+"_zcmpnt";
                        model.addConstr(x_comp == 0.0, x_constr);
                        model.addConstr(y_comp == 0.0, y_constr);
                        model.addConstr((z_comp - load) == 0.0, z_constr);
                    } else {
                        model.addConstr(z_i[i + 1] == 0.0);
                        // pinned_indices.push_back(i + 1); // visualization purposes
                    }
                }
                model.update();
                model.feasRelax(0, false, false, true);
                model.write("debug.lp");
                model.optimize();

                float solved_z[num_nodes + 1]; // unneeded, but just for understanding
                float sum_difference = 0;
                for(int i = 0; i < num_nodes; i++) {
                    solved_z[i + 1] = z_i[i + 1].get(GRB_DoubleAttr_X);
                    if (iters > 1) { // no valid comparison atp
                        sum_difference += std::abs(solved_z[i + 1] - subs_z[i + 1]);
                    }
                    subs_z[i + 1] = solved_z[i + 1];
                }
                z_abs_diff = sum_difference;
                toggle_solve_scalar = true;
                toggle_solve_z = false;
            }
            else if (toggle_solve_scalar) {
                cout << "solving for scalar" << endl;
                GRBEnv env = GRBEnv();
                GRBModel model = GRBModel(env);

                GRBVar scalars[edges.size()];
                set<pair<int, int>>::iterator itr;
                int i = 0;
                for (itr = edges.begin(); itr != edges.end(); itr++) {
                    auto edge = *itr;
                    string s_i = "s_" + to_string(edge.first) + "_" + to_string(edge.second);
                    scalars[i++] = model.addVar(1.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, s_i);
                }

                for (int i = 0; i < num_nodes; i++) {
                    set<int> neighbors_i = neighbors[i + 1];
                    set<int>::iterator itr;
                    GRBLinExpr x_comp;
                    GRBLinExpr y_comp;
                    GRBLinExpr z_comp;

                    float e_ix = x_i[i + 1];
                    float e_iy = y_i[i + 1];
                    float e_iz = subs_z[i + 1]; // use solved z 

                    auto pos = pinned_indices.find(i + 1);
                    if (pos == pinned_indices.end()) {
                        for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                            int neighbor = *itr;
                            pair<int, int> edge = make_edge(i + 1, neighbor);
                            float e_jx = x_i[neighbor];
                            float e_jy = y_i[neighbor];
                            float e_jz = subs_z[neighbor]; // use solved z
                            GRBVar s = scalars[scalar_idx[edge]];
                            x_comp += (e_ix - e_jx) * s;
                            y_comp += (e_iy - e_jy) * s;
                            z_comp += (e_iz - e_jz) * s;
                        }
                        string x_constr = to_string(i + 1)+"_xcmpnt";
                        string y_constr = to_string(i + 1)+"_ycmpnt";
                        string z_constr = to_string(i + 1)+"_zcmpnt";
                        model.addConstr(x_comp == 0.0, x_constr);
                        model.addConstr(y_comp == 0.0, y_constr);
                        model.addConstr((z_comp - load) == 0.0, z_constr);
                    } else {
                        // model.addConstr(z_i[i + 1] == 0.0); <- not needed
                        // pinned_indices.push_back(i + 1); // visualization purposes
                    }
                }
                model.update();
                model.feasRelax(0, false, false, true);
                model.write("debug.lp");
                model.optimize();

                float solved_scalar[edges.size()]; // unneeded, but just for understanding
                float sum_difference = 0;
                for(int i = 0; i < edges.size(); i++) {
                    solved_scalar[i] = scalars[i].get(GRB_DoubleAttr_X);
                    sum_difference += std::abs(solved_scalar[i] - subs_scalars[i]);
                    subs_scalars[i] = solved_scalar[i];
                }
                scalar_abs_diff = sum_difference;
                toggle_solve_scalar = false;
                toggle_solve_z = true;
            } else {
                cout << "ARRIVED AT INVALID STATE" << endl;
            } 
            if (iters > 100) {
                cout << "Hit 100 iters, terminating now" << endl;
                toggling = false;
            }

            // Checking for convergence
            float z_update_diff = std::abs(prev_z_abs_diff - z_abs_diff);
            float scalar_update_diff = std::abs(prev_scalar_abs_diff - scalar_abs_diff);
            if (z_update_diff < 1e-8 && scalar_update_diff < 1e-8) {
                cout << "Converged, terminating now" << endl;
                cout << "Num iters: " << iters << endl;
                cout << z_update_diff << endl;
                cout << scalar_update_diff << endl;
                toggling = false;
            } else {
                cout << "Didn't converge" << endl;
                cout << z_update_diff + scalar_update_diff << endl;
            }
            prev_z_abs_diff = z_abs_diff;
            prev_scalar_abs_diff = scalar_abs_diff;
        }

        for (int i = 0; i < num_nodes; i++) {
            float x = x_i[i + 1];
            float z = subs_z[i + 1];
            cout << "node " << subs_z[i + 1] << endl;
            float y = y_i[i + 1];
            points.push_back(
                glm::vec3{x, z, y}
            );
        }

        
        for (int i = 0; i < num_nodes; i++) {
            cout << "node " << i << " ";
            set<int> neighbors_i = neighbors[i + 1];
            set<int>::iterator itr;
            for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) { 
                int neighbor = *itr;
                pair<int, int> edge = make_edge(i + 1, neighbor);
                float s = subs_scalars[scalar_idx[edge]];
                cout << s << " ";
            }
            cout << endl;
        }

        for (int i = 0; i < edges.size(); i++) {
            cout << "scalar " << subs_scalars[i] << endl;
        }
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    vector<glm::vec3> pinned;
    for (auto& index : pinned_indices) {
        pinned.push_back(glm::vec3{x_i[index], 0, y_i[index]});
    }

    polyscope::PointCloud* pinnedCloud = polyscope::registerPointCloud("pinned vertices", pinned);
    pinnedCloud->setPointRenderMode(polyscope::PointRenderMode::Quad);

    vector<array<int, 2>> vec_edges;
    for (auto &p : edges) {
        array<int, 2> arr = {p.first - 1, p.second - 1};
        vec_edges.push_back(arr);
    }

    if (points.size() == 0) {
        for (int i = 0; i < num_nodes; i++) {
            int x = x_i[i + 1];
            int y = y_i[i + 1];
            points.push_back(glm::vec3{x, 0, y});
        }
    }

    polyscope::registerCurveNetwork("my network", points, vec_edges);

    polyscope::show();

    return 0;
}