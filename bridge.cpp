#include <unordered_map>
#include <iostream>
#include <fstream>

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "polyscope/pick.h"

#include "gurobi_c++.h"
#include "gurobi_c.h"

using namespace std;

// // // // // // // // // // // // // // // // 
// No existing C++ default pair<T, T> hash   //
// // // // // // // // // // // // // // // //

struct PairHash {
    template <typename T1, typename T2>
    size_t operator()(const pair<T1, T2>& p) const {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2; // Combining the hash values
    }
};

// // // // // // // // // // // // // // // // // //
//  Piping required to use the Triangle library.   //
// // // // // // // // // // // // // // // // // //

#define REAL double
#define VOID int

extern "C" void triangulate(char *, struct triangulateio *, struct triangulateio *, struct triangulateio *);

extern "C" struct triangulateio {
  REAL *pointlist;                                               /* In / out */
  REAL *pointattributelist;                                      /* In / out */
  int *pointmarkerlist;                                          /* In / out */
  int numberofpoints;                                            /* In / out */
  int numberofpointattributes;                                   /* In / out */

  int *trianglelist;                                             /* In / out */
  REAL *triangleattributelist;                                   /* In / out */
  REAL *trianglearealist;                                         /* In only */
  int *neighborlist;                                             /* Out only */
  int numberoftriangles;                                         /* In / out */
  int numberofcorners;                                           /* In / out */
  int numberoftriangleattributes;                                /* In / out */

  int *segmentlist;                                              /* In / out */
  int *segmentmarkerlist;                                        /* In / out */
  int numberofsegments;                                          /* In / out */

  REAL *holelist;                        /* In / pointer to array copied out */
  int numberofholes;                                      /* In / copied out */

  REAL *regionlist;                      /* In / pointer to array copied out */
  int numberofregions;                                    /* In / copied out */

  int *edgelist;                                                 /* Out only */
  int *edgemarkerlist;            /* Not used with Voronoi diagram; out only */
  REAL *normlist;                /* Used only with Voronoi diagram; out only */
  int numberofedges;                                             /* Out only */
};

// // // // // // // // // //
//   End Triangle pipes.   //
// // // // // // // // // //

/********************************************************************************************************
    We now want to define 2D closed boundaries in the x-y plane.
    We construct this boundary using a set of Bezier curves that form a closed loop, defined as follows:
    ∀ curves i with endpoints (a, b)->(c, d): ∃ curves x, y with endpoints (a, b) and (c, d)
*********************************************************************************************************/
vector<array<pair<double, double>, 3>> triangles;
unordered_map<pair<double, double>, set<int>, PairHash> point_tri_map;
unordered_map<int, pair<double, double>> tri_centroid;
double increment = 0.1;
bool pinSelectMode = false;

struct Bezier {
  // User-defined parameters
  vector<pair<double, double>> controlPoints;
  bool pinnedCurve;
  // Internal points plotted along curve 
  vector<pair<double, double>> points;
};

/*
    Algorithm below implements calculating a bezier curve given a set of control points.
    Input: Bezier curve defined with a vector of control points.
    Output: void, sets `points` element in `Bezier` struct.
*/
pair<double, double> deCastelJau(const vector<pair<double, double>>& controlPoints, double t) {
  if (controlPoints.size() == 1) {
    return controlPoints[0];
  }
  vector<pair<double, double>> intermediatePoints;
  for (int i = 0; i < controlPoints.size() - 1; i++) {
    double x = (1 - t) * controlPoints[i].first + t * controlPoints[i + 1].first;
    double y = (1 - t) * controlPoints[i].second + t * controlPoints[i + 1].second;
    intermediatePoints.push_back({x, y});
  }
  return deCastelJau(intermediatePoints, t);
}

void bezierCurve(Bezier* curve) {
  vector<pair<double, double>> p;
  for (double t = 0; t <= 1.0; t += increment) {
    p.push_back(deCastelJau(curve->controlPoints, t));
  }
  curve->points = p;
}

/*
    We define our boundary by a set of Bezier curves.
    This method will calculate the points along each Bezier curves.
*/
void defineBezierCurves(vector<Bezier*> curves) {
    for(int i = 0; i < curves.size(); i++) {
        bezierCurve(curves[i]);

        // Assumptions: first and last control points are end points.
        // Goal: Amend lost precision in bezier curve fitting.
        vector<pair<double, double>> controlPoints = curves[i]->controlPoints;
        vector<pair<double, double>> points = curves[i]->points;
        points[0] = controlPoints[0];
        points[points.size() - 1] = controlPoints[controlPoints.size() - 1];
        curves[i]->points = points;
    }
}

/*
    Given a set of defined Bezier curves, we want to stitch them together and construct the closed shape.
    - Note: The present definition of "closed" does not ensure that curves do not cross.
            If the defined curves cross, this may result in error.
    Input: vector of Bezier curves
    Output: array of points, array of point indices where every 2 elements represent a segment
*/
void constructBoundary(vector<Bezier*> curves, vector<pair<double, double>>& boundaryPoints, vector<pair<int, int>>& boundarySegments, int& boundarySize) {
    defineBezierCurves(curves);

    // Stitch curves together by finding the neighbors of each point.
    // Each point has exactly 2 neighbors. -> Used to ensure "closed" boundary property.
    unordered_map<pair<double, double>, vector<pair<double, double>>, PairHash> pointNeighbors;
    for (int i = 0; i < curves.size(); i++) {
        vector<pair<double, double>> points = curves[i]->points;
        pair<double, double> prev_point = points[0];
        for (int j = 1; j < points.size(); j++) {
            pair<double, double> cur_point = points[j];
            pointNeighbors[cur_point].push_back(prev_point);
            pointNeighbors[prev_point].push_back(cur_point);
            prev_point = cur_point;
        }
    }

    // Looping over our pointNeighbors maps allows us to do the following:
    // (1) Construct `boundaryPoints` array.
    // (2) Check "closed" loop condition.
    // (3) Map points to indices s.t. can construct `boundarySegments` array.
    vector<pair<double, double>> bPoints;
    bool exactlyTwoNeighbors = true;
    unordered_map<pair<double, double>, int, PairHash> point_idx_map;
    int p_index = 0;
    for (const auto& [point, neighbors] : pointNeighbors) {
        bPoints.push_back(point);
        exactlyTwoNeighbors &= (neighbors.size() == 2);
        point_idx_map[point] = p_index++;
    }
    if (!exactlyTwoNeighbors) 
        cout << "BROKE!" << endl; // TODO: replace with exception handling

    //
    // We will now traverse our points like a circular linked list, to construct our segments.
    //
    vector<pair<int, int>> bSegments;
    int start_idx = 0;
    pair<double, double> first {bPoints[start_idx].first, bPoints[start_idx].second};
    pair<double, double> cur {bPoints[start_idx].first, bPoints[start_idx].second};
    pair<double, double> prev {bPoints[start_idx].first, bPoints[start_idx].second};
    int seg_idx = 0;
    do {
        // Note: exactly 2 neighbors per point
        pair<double, double> n_1 = pointNeighbors[cur][0];
        pair<double, double> n_2 = pointNeighbors[cur][1]; 

        // We arrived cur from either neighbor 1 or neighbor 2, do not want to revisit points
        if (n_1 != prev) {
            bSegments.push_back({point_idx_map[cur], point_idx_map[n_1]});
            prev = cur;
            cur = n_1;
        } else {
            bSegments.push_back({point_idx_map[cur], point_idx_map[n_2]});
            prev = cur;
            cur = n_2;
        }
        seg_idx++;
    } while (cur != first);

    boundaryPoints = bPoints;
    boundarySegments = bSegments;
    boundarySize = (int)pointNeighbors.size();
}

pair<double, double> triCentroid(pair<double, double> p1, pair<double, double> p2, pair<double, double> p3) {
    pair<double, double> p = pair<double, double>{(p1.first + p2.first + p3.first)/3.0, (p1.second + p2.second + p3.second)/3.0};
    return p;
}

void triangulateBoundary(vector<pair<double, double>> points, vector<pair<int, int>> segments, int size, bool q, bool D, double maxArea, bool Y, int maxS,
                         vector<pair<double, double>>& points_out, vector<pair<int, int>>& edges_out) {
    double* points_in = new double[2 * points.size()];
    int* segs_in = new int[2 * segments.size()];
    for (int i = 0; i < points.size(); i++) {
        points_in[2 * i] = points[i].first;
        points_in[2 * i + 1] = points[i].second;
    }
    for (int i = 0; i < segments.size(); i++) {
        segs_in[2 * i] = segments[i].first;
        segs_in[2 * i + 1] = segments[i].second;
    }

    struct triangulateio* trio_in = new triangulateio;
    trio_in->pointlist = points_in;
    trio_in->numberofpoints = size;
    trio_in->numberofpointattributes = 0;
    trio_in->pointmarkerlist = nullptr;
    trio_in->segmentlist = segs_in;
    trio_in->numberofsegments = size;
    trio_in->segmentmarkerlist = nullptr;
    trio_in->numberofholes = 0;
    trio_in->holelist = nullptr;
    trio_in->numberofregions = 0;
    trio_in->regionlist = nullptr;

    struct triangulateio* trio_out = new triangulateio;
    trio_out->pointlist = nullptr;
    trio_out->trianglelist = nullptr;
    trio_out->segmentlist = nullptr;
    trio_out->segmentmarkerlist = nullptr;
    trio_out->pointmarkerlist = nullptr;

    string triswitches = "pz";
    if (q) triswitches += "q";
    if (D) triswitches += "D";
    if (maxArea != 0.0) triswitches += "a" + to_string(maxArea);
    if (Y) triswitches += "Y";
    if (maxS > -1) triswitches += "S" + to_string(maxS);
    triangulate(const_cast<char*>(triswitches.c_str()), trio_in, trio_out, nullptr);
    vector<pair<double, double>> p_out;
    for (int i = 0; i < trio_out->numberofpoints; i++) {
        p_out.push_back(pair<double, double>{trio_out->pointlist[2 * i], trio_out->pointlist[2 * i + 1]});
    }
    vector<pair<int, int>> e_out;
    for (int i = 0; i < trio_out->numberoftriangles; i++) {
        e_out.push_back(pair<int, int>{trio_out->trianglelist[3 * i] , trio_out->trianglelist[3 * i + 1]});
        e_out.push_back(pair<int, int>{trio_out->trianglelist[3 * i], trio_out->trianglelist[3 * i + 2]});
        e_out.push_back(pair<int, int>{trio_out->trianglelist[3 * i + 1], trio_out->trianglelist[3 * i + 2]});
        int v1 = trio_out->trianglelist[3 * i];
        int v2 = trio_out->trianglelist[3 * i + 1];
        int v3 = trio_out->trianglelist[3 * i + 2];
        pair<double, double> t_1 = p_out[v1];
        pair<double, double> t_2 = p_out[v2];
        pair<double, double> t_3 = p_out[v3];
        triangles.push_back(array<pair<double, double>, 3>{t_1, t_2, t_3});
        
        int tri_idx = triangles.size() - 1;
        
        point_tri_map[t_1].insert(tri_idx);
        point_tri_map[t_2].insert(tri_idx);
        point_tri_map[t_3].insert(tri_idx);

        tri_centroid[tri_idx] = triCentroid(t_1, t_2, t_3);
    }
    points_out = p_out;
    edges_out = e_out;
}

// vector<Bezier*> curves = {new Bezier{{{0, 0}, {6, -3}, {12, 0}}}, 
//   new Bezier{{{12, 0}, {15, 9}}},
//   new Bezier{{{15, 9}, {10.5, 12}, {6, 15}}},
//   new Bezier{{{6, 15}, {-3, 9}}},
//   new Bezier{{{-3, 9}, {0, 0}}}};

// vector<Bezier*> curves = {new Bezier{{{0, 0}, {0, 10}}}, 
//   new Bezier{{{0, 10}, {10, 10}}},
//   new Bezier{{{10, 10}, {10, 0}}},
//   new Bezier{{{10, 0}, {0, 0}}}};

// vector<Bezier*> curves = {new Bezier{{{210, 210}, {-210, 210}}},
//     new Bezier{{{-210, 210}, {-210, -210}}},
//     new Bezier{{{-210, -210}, {210, -210}}},
//     new Bezier{{{210, -210}, {210, 210}}}};

vector<Bezier*> curves = {new Bezier{{{-210, 123}, {0, 56}, {210, 123}}},
    new Bezier{{{210, 123}, {210, -67}}},
    new Bezier{{{210, -67}, {0, 0}, {-210, -67}}},
    new Bezier{{{-210, -67}, {-210, 123}}}};
    
bool curvesReady;
vector<pair<double, double>> points;
vector<pair<int, int>> segments;
int in_size;
bool q = false;
bool D = false;
double maxArea = 0.0;
bool Y = false;
int maxS = -1;
vector<pair<double, double>> p_out;
vector<pair<int, int>> e_out;
vector<pair<double, double>> pinned;
bool triReady = false;
float bridge_load = 2000;

void visualizeBoundary() {
    vector<glm::vec3> bPoints;
    for (int i = 0; i < points.size(); i++) {
        glm::vec3 pnt{points[i].first, 0, points[i].second};
        bPoints.push_back(pnt);
    }
    polyscope::registerPointCloud("Boundary", bPoints);
}

void visualizePinned() {
    vector<glm::vec3> pinnedPoints;
    for (int i = 0; i < pinned.size(); i++) {
        glm::vec3 pnt{pinned[i].first, 0, pinned[i].second};
        pinnedPoints.push_back(pnt);
    }
    for (int i = 0; i < curves.size(); i++) {
        if (curves[i]->pinnedCurve) {
            vector<pair<double, double>> p = curves[i]->points;
            for (int j = 0; j < p.size(); j++) {
                glm::vec3 pnt{p[j].first, 0, p[j].second};
                pinnedPoints.push_back(pnt);
                pinned.push_back(p[j]);
            }
        }
    }
    polyscope::registerPointCloud("Pinned Points", pinnedPoints);
}

// void visualizeBoundary() {
//     vector<glm::vec3> bPoints;
//     vector<glm::vec3> pinnedPoints;
//     for (int i = 0; i < curves.size(); i++) {
//         vector<pair<double, double>> p = curves[i]->points;
//         for (int j = 0; j < p.size(); j++) {
//             glm::vec3 pnt{p[j].first, 0, p[j].second};
//             bPoints.push_back(pnt);
//             if (curves[i]->pinnedCurve) {
//                 pinnedPoints.push_back(pnt);
//                 pinned.push_back(p[j]);
//             }
//         }
//     }
//     polyscope::registerPointCloud("Boundary", bPoints);
//     polyscope::registerPointCloud("Pinned Points", pinnedPoints);
// }

void updateBoundary() {
    vector<glm::vec3> pinnedPoints;
    for (auto &p : pinned) {
        glm::vec3 pnt{p.first, 0, p.second};
        pinnedPoints.push_back(pnt);
    }
    polyscope::registerPointCloud("Pinned Points", pinnedPoints);
}

void visualizeTriangulation() {
    vector<glm::vec3> tPoints;
    vector<array<int, 2>> tEdges;
    for (int i = 0; i < p_out.size(); i++) {
        tPoints.push_back(glm::vec3{p_out[i].first, 0, p_out[i].second});
    }
    for (int i = 0; i < e_out.size(); i++) {
        tEdges.push_back({e_out[i].first, e_out[i].second});
    }
    polyscope::registerCurveNetwork("Triangulated Boundary", tPoints, tEdges);

    vector<array<int, 2>> dual_edges;
    vector<glm::vec3> dual_points;
    for (int i = 0; i < triangles.size(); i++) {
        dual_points.push_back(glm::vec3{1.05 * (tri_centroid[i].first), 0, 1.05 * (tri_centroid[i].second) - 300});
    }
    for (int i = 0; i < p_out.size(); i++) {
        set<int> triNeighbors = point_tri_map[p_out[i]];
        // cout << triNeighbors.size() << endl;
        set<int>::iterator itr;
        multiset<pair<double, int>> mset;
        for (itr = triNeighbors.begin(); itr != triNeighbors.end(); itr++) {
            int tri_idx = *itr;
            pair<double, double> centroid = tri_centroid[tri_idx];
            double angle = atan2((centroid.second - p_out[i].second), (centroid.first - p_out[i].first));
            mset.insert(pair<double, int>{angle, tri_idx});
        }

        vector<int> in_order;
        for (const auto& p: mset) {
            in_order.push_back(p.second);
        }
        for (int j = 0; j < in_order.size() - 1; j++) {
            dual_edges.push_back({in_order[j], in_order[j + 1]});
        }
        dual_edges.push_back({in_order[in_order.size() - 1], in_order[0]});
    }
    // polyscope::registerCurveNetwork("Dual figure", dual_points, dual_edges);
}

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
    return 1;
    static std::default_random_engine e;
    static std::uniform_real_distribution<> dist(1, 11);
    return dist(e);
}

void debugConstraints(string debugFileName, map<int, set<int>> neighbors, set<int> pinned, map<pair<int, int>, int> scalar_idx, 
                        vector<double> subs_scalars, vector<double> subs_z, unordered_map<int, double> x_i, unordered_map<int, double> y_i, int iter) {
    
    std::fstream debug(debugFileName, std::ios::out | std::ios::app);
    double load = bridge_load / x_i.size();
    double check_x = 0;
    double check_y = 0;
    double check_z = 0;
    for (int i = 0; i < x_i.size(); i++) {
        set<int> neighbors_i = neighbors[i];
        set<int>::iterator itr;
        double e_ix = x_i[i];
        double e_iy = y_i[i];
        double e_iz = subs_z[i];

        double x_comp = 0;
        double y_comp = 0;
        double z_comp = 0;
        if (pinned.find(i) == pinned.end()) {
            for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                int neighbor = *itr;
                pair<int, int> edge = make_edge(i, neighbor);
                double e_jx = x_i[neighbor];
                double e_jy = y_i[neighbor];
                double e_jz = subs_z[neighbor];
                double s = subs_scalars[scalar_idx[edge]];

                x_comp += (e_ix - e_jx) * s;
                y_comp += (e_iy - e_jy) * s;
                z_comp += (e_iz - e_jz) * s;
            }
            check_x += abs(x_comp);
            check_y += abs(y_comp);
            check_z += abs(z_comp - load);
        }
    }
    debug << "ITERATION: " << iter << endl;
    debug << "constraint error in x: " << check_x << endl;
    debug << "constraint error in y: " << check_y << endl;
    debug << "constraint error in z: " << check_z << endl;
}

vector<float> visualizeForce(map<int, set<int>> neighbors, set<int> pinned, map<pair<int, int>, int> scalar_idx, 
                        vector<double> subs_scalars, vector<double> subs_z, vector<double> x_i, vector<double> y_i) {
    
    double load = bridge_load / x_i.size();
    double check_x = 0;
    double check_y = 0;
    double check_z = 0;
    vector<float> forces;
    for (int i = 0; i < x_i.size(); i++) {
        set<int> neighbors_i = neighbors[i];
        set<int>::iterator itr;
        double e_ix = x_i[i];
        double e_iy = y_i[i];
        double e_iz = subs_z[i];

        double x_comp = 0;
        double y_comp = 0;
        double z_comp = 0;
        if (pinned.find(i) == pinned.end()) {
            for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                int neighbor = *itr;
                pair<int, int> edge = make_edge(i, neighbor);
                double e_jx = x_i[neighbor];
                double e_jy = y_i[neighbor];
                double e_jz = subs_z[neighbor];
                double s = subs_scalars[scalar_idx[edge]];

                x_comp += (e_ix - e_jx) * s;
                y_comp += (e_iy - e_jy) * s;
                z_comp += (e_iz - e_jz) * s;
            }
            float forceResidual = abs(x_comp) + abs(y_comp) + abs(z_comp - load);
            forces.push_back(forceResidual);
        } else {
            forces.push_back(0);
        }
    }
    return forces;
}

void constructBridge() {
    vector<double> x_i;
    vector<double> y_i;
    unordered_map<pair<double, double>, int, PairHash> pnt_idx_map;
    for (int i = 0; i < p_out.size(); i++) {
        x_i.push_back(p_out[i].first);
        y_i.push_back(p_out[i].second);
        pnt_idx_map[p_out[i]] = i;
    }

    map<int, set<int>> neighbors;
    set<pair<int, int>> edges;
    for (int i = 0; i < e_out.size(); i++) {
        edges.insert(make_edge(e_out[i].first, e_out[i].second));
        edges.insert(make_edge(e_out[i].second, e_out[i].first));
        add_neighbor(neighbors, e_out[i].first, e_out[i].second);
        add_neighbor(neighbors, e_out[i].second, e_out[i].first);
    }

    set<int> pinned_idx;
    cout << "PINNED" << endl;
    cout << pinned.size() << endl;
    for (int i = 0; i < pinned.size(); i++) {
        cout << pinned[i].first << " " << pinned[i].second << ": " << pnt_idx_map.count(pinned[i]) << endl;
        pinned_idx.insert(pnt_idx_map[pinned[i]]);
        // cout << "pinned " << pnt_idx_map[pinned[i]] << endl;
    }
    cout << pinned_idx.size() << endl;

    double load = bridge_load / (x_i.size());

    set<pair<int, int>>::iterator itr;
    map<pair<int, int>, int> scalar_idx;
    int i = 0;
    for (itr = edges.begin(); itr != edges.end(); itr++) {
        auto edge = *itr;
        scalar_idx[edge] = i++;
    }

    vector<double> subs_scalars;
    for (int i = 0; i < edges.size(); i++) {
        subs_scalars.push_back(1.0);
    }
    vector<double> subs_z(x_i.size(), 0.0);

    vector<float> forceResiduals;
    try {
        int iters = 0;
        float prev_z_abs_diff = 1e8;
        float prev_scalar_abs_diff = 1e8;

        while (true) {
            float z_abs_diff = 0;
            float scalar_abs_diff = 0;

            cout << "Z STEP" << endl;

            GRBEnv z_env = GRBEnv();
            z_env.set(GRB_DoubleParam_TimeLimit, 15);
            GRBModel z_model = GRBModel(z_env);

            GRBVar z_i[x_i.size()];
            for (int i = 0; i < x_i.size(); i++) {
                string z = "z_" + to_string(i);
                z_i[i] = z_model.addVar(-GRB_INFINITY, GRB_INFINITY, 0.0, GRB_CONTINUOUS, z);
            }
            GRBQuadExpr forceResidualsZ;
            for (int i = 0; i < x_i.size(); i++) {
                set<int> neighbors_i = neighbors[i];
                set<int>::iterator itr;
                GRBLinExpr x_comp;
                GRBLinExpr y_comp;
                GRBLinExpr z_comp;

                float e_ix = x_i[i];
                float e_iy = y_i[i];
                GRBVar e_iz = z_i[i];
                if (pinned_idx.find(i) == pinned_idx.end()) {
                    for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                        int neighbor = *itr;
                        pair<int, int> edge = make_edge(i, neighbor);
                        float s = subs_scalars[scalar_idx[edge]];

                        float e_jx = x_i[neighbor];
                        float e_jy = y_i[neighbor];
                        GRBVar e_jz = z_i[neighbor];

                        x_comp += (e_ix - e_jx) * s;
                        y_comp += (e_iy - e_jy) * s;
                        z_comp += (e_iz - e_jz) * s;
                    }
                    forceResidualsZ += (x_comp * x_comp) + (y_comp * y_comp) + ((z_comp - load) * (z_comp - load));
                } else {
                    z_model.addConstr(z_i[i] == 0.0);
                }
            }
            z_model.setObjective(forceResidualsZ);
            z_model.update();
            z_model.write("debug.lp");
            z_model.optimize();
            if (z_model.get(GRB_IntAttr_Status) == 9) {
                forceResiduals = visualizeForce(neighbors, pinned_idx, scalar_idx, subs_scalars, subs_z, x_i, y_i);
                break;
            }

            for(int i = 0; i < x_i.size(); i++) {
                subs_z[i] = z_i[i].get(GRB_DoubleAttr_X);
            }
            
            cout << "SCALAR STEP" << endl;

            GRBEnv s_env = GRBEnv();
            s_env.set(GRB_DoubleParam_TimeLimit, 15);
            GRBModel s_model = GRBModel(s_env);

            GRBVar scalars[edges.size()];
            set<pair<int, int>>::iterator itr;
            int i = 0;
            for (itr = edges.begin(); itr != edges.end(); itr++) {
                auto edge = *itr;
                string s_i = "s_" + to_string(edge.first) + "_" + to_string(edge.second);
                scalars[i] = s_model.addVar(0.0, GRB_INFINITY, 0.0, GRB_CONTINUOUS, s_i);
                scalars[i].set(GRB_DoubleAttr_Start, subs_scalars[i]);
                i++;
            }

            GRBQuadExpr forceResidualsS;
            for (int i = 0; i < x_i.size(); i++) {
                set<int> neighbors_i = neighbors[i];
                set<int>::iterator itr;
                GRBLinExpr x_comp;
                GRBLinExpr y_comp;
                GRBLinExpr z_comp;
                GRBLinExpr scalar_const;

                float e_ix = x_i[i];
                float e_iy = y_i[i];
                float e_iz = subs_z[i];

                if (pinned_idx.find(i) == pinned_idx.end()) {
                    for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                        int neighbor = *itr;
                        pair<int, int> edge = make_edge(i, neighbor);
                        GRBVar s = scalars[scalar_idx[edge]];

                        float e_jx = x_i[neighbor];
                        float e_jy = y_i[neighbor];
                        float e_jz = subs_z[neighbor];
                        
                        x_comp += (e_ix - e_jx) * s;
                        y_comp += (e_iy - e_jy) * s;
                        z_comp += (e_iz - e_jz) * s;
                        scalar_const += s;
                        forceResidualsS += 1e-9 * s * s;
                    }
                    string scalar_constr = to_string(i)+"_scalar";
                    s_model.addConstr(scalar_const >= .1, scalar_constr);
                    forceResidualsS += (x_comp * x_comp) + (y_comp * y_comp) + ((z_comp - load) * (z_comp - load));
                } 
            }
            s_model.setObjective(forceResidualsS);
            s_model.update();
            s_model.write("debug.lp");
            s_model.set(GRB_IntParam_Presolve, 0);
            s_model.optimize();
            if (s_model.get(GRB_IntAttr_Status) == 9) {
                forceResiduals = visualizeForce(neighbors, pinned_idx, scalar_idx, subs_scalars, subs_z, x_i, y_i);
                break;
            }

            float solved_scalar[edges.size()]; // unneeded, but just for understanding

            i = 0;
            for(itr = edges.begin(); itr != edges.end(); itr++) {
                subs_scalars[i] = scalars[i].get(GRB_DoubleAttr_X);
                i++;
            }

            if (s_model.getObjective().getValue() < 1e-3) {
                cout << "Objective is < .001" << endl;
                forceResiduals = visualizeForce(neighbors, pinned_idx, scalar_idx, subs_scalars, subs_z, x_i, y_i);
                break;
            }
            if (iters > 1000) {
                forceResiduals = visualizeForce(neighbors, pinned_idx, scalar_idx, subs_scalars, subs_z, x_i, y_i);
                break;
            }   
            prev_z_abs_diff = z_abs_diff;
            prev_scalar_abs_diff = scalar_abs_diff;
            iters++;
        }
    } catch(GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        if (e.getErrorCode() == 9) {
            forceResiduals = visualizeForce(neighbors, pinned_idx, scalar_idx, subs_scalars, subs_z, x_i, y_i);
        }
        cout << e.getMessage() << endl;
    } catch(...) {
        cout << "Exception during optimization" << endl;
    }

    vector<glm::vec3> bridge_points;
    for (int i = 0; i < x_i.size(); i++) {
        float x = x_i[i];
        float z = subs_z[i];
        float y = y_i[i];
        bridge_points.push_back(
            glm::vec3{x, z, y}
        );
    }

    vector<array<int, 2>> vec_edges;
    for (auto &p : edges) {
        array<int, 2> arr = {p.first, p.second};
        vec_edges.push_back(arr);
    }
    std::vector<std::array<double, 3>> edgeColor(edges.size());
    i = 0;
    for (itr = edges.begin(); itr != edges.end(); itr++) {
        edgeColor[i] = {subs_scalars[i], 0, 0};
        i++;
    }
    std::vector<std::array<double, 3>> forceColor(x_i.size());
    for (size_t i = 0; i < x_i.size(); i++) {
        forceColor[i] = {0, forceResiduals[i], 0};
    }
    polyscope::registerPointCloud("network points", bridge_points);
    polyscope::getPointCloud("network points")->addColorQuantity("forces", forceColor);

    polyscope::registerCurveNetwork("my network", bridge_points, vec_edges);
    polyscope::getCurveNetwork("my network")->addEdgeColorQuantity("scalars", edgeColor);
}

void drawImGui() {
    ImGuiIO& io = ImGui::GetIO();

    if (ImGui::Checkbox("Pin Select Mode", &pinSelectMode)) {
        
    }
    if (triReady && pinSelectMode) {
        vector<glm::vec3> triPoints;
        for (int i = 0; i < p_out.size(); i++) {
            triPoints.push_back({p_out[i].first, 0, p_out[i].second});
        }
        polyscope::registerPointCloud("Triangle Points", triPoints);
        if (io.MouseClicked[0]) { // if the left mouse button was clicked
            // gather values
            glm::vec2 screenCoords{ io.MousePos.x * io.DisplayFramebufferScale.x, io.MousePos.y * io.DisplayFramebufferScale.y}; 
            glm::vec3 worldRay = polyscope::view::screenCoordsToWorldRay(screenCoords);
            glm::vec3 worldPos = polyscope::view::screenCoordsToWorldPosition(screenCoords);
            std::pair<polyscope::Structure*, size_t> pickPair =
                polyscope::pick::evaluatePickQuery(screenCoords.x, screenCoords.y);
            if (pickPair.first != nullptr) {
                // std::cout << "    structure: " << "none" << std::endl;
                std::cout << "    structure: " << pickPair.first << " element id: " << pickPair.second << std::endl;
                if (pickPair.first->getName() == "Triangle Points") {
                    cout << "here?" << endl;
                    pinned.push_back(p_out[pickPair.second]);
                }
                // cout << pickPair.first->typeName() << endl;
                updateBoundary();
            }
        }
    }

    if (ImGui::Button("Add Curve")) {
        curves.push_back(new Bezier());
    }
    for (int i = 0; i < curves.size(); i++) {
        if (ImGui::CollapsingHeader(("Curve " + std::to_string(i)).c_str())) {
            if (ImGui::Checkbox(("Pin edge##_" + to_string(i)).c_str(), &curves[i]->pinnedCurve)) {

            }
            if (ImGui::Button(("Delete##_" + to_string(i)).c_str())) {
                curves.erase(curves.begin() + i);
                continue;
            }
            if (ImGui::Button(("Add Control Point##_" + to_string(i)).c_str())) {
                curves[i]->controlPoints.push_back({0.0, 0.0});
            }
            for (int j = 0; j < curves[i]->controlPoints.size(); j++) {
                ImGui::Text("Control Point %d", j);
                ImGui::InputDouble((std::to_string(i) + "_" + std::to_string(j)+ "_x").c_str(), &curves[i]->controlPoints[j].first);
                ImGui::InputDouble((std::to_string(i) + "_" + std::to_string(j) + "_y").c_str(), &curves[i]->controlPoints[j].second);
            }
        }
    }

    ImGui::Separator();
    ImGui::InputDouble("Bezier increment granularity", &increment);
    if (ImGui::Button("Construct Boundary")) {
        constructBoundary(curves, points, segments, in_size);
        curvesReady = true;
        visualizeBoundary();
        visualizePinned();
    }

    ImGui::Separator();
    
    if(ImGui::CollapsingHeader("Triangulation Parameters")) {
        if (ImGui::Checkbox("No angles smaller than 20 degrees", &q)) {
        }
        if (ImGui::Checkbox("Conforming Delaunay", &D)) {
        }
        ImGui::InputDouble("Maximum triangle area", &maxArea);
        if (ImGui::Checkbox("No added Steiner points along boundary", &Y)) {
        }
        ImGui::InputInt("Max number of added Steiner points", &maxS);
    }
    if (ImGui::Button("Triangulate")) {
        if (curvesReady) {
            triangulateBoundary(points, segments, in_size, q, D, maxArea, Y, maxS, p_out, e_out);
            visualizeTriangulation();
            triReady = true;
        }
    }
    ImGui::Separator();
    ImGui::InputFloat("Total vertical load on bridge", &bridge_load);
    if (ImGui::Button("Construct Bridge")) {
        if (triReady) {
            constructBridge();
        }
    }
}

int main() {
  polyscope::init();

  polyscope::state::userCallback = drawImGui;

  polyscope::show();
}