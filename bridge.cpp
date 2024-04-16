#include <unordered_map>
#include <iostream>
#include <fstream>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

#include "gurobi_c++.h"
#include "gurobi_c.h"

using std::cout;
using std::endl;
using std::vector;
using std::array;
using std::pair;
using std::unordered_map;
using std::set;
using std::map;
using std::string;
using std::ofstream;
using std::to_string;
using std::multiset;
using std::hash;


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
  if (increment <= 0 || increment >= 1) increment = 0.1;
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
    // ofstream myfile;
    // myfile.open ("test-xy.txt");
    // myfile << trio_out->numberofpoints << " " << 2 << " " << 0 << " " << 0 << endl;
    for (int i = 0; i < trio_out->numberofpoints; i++) {
        p_out.push_back(pair<double, double>{trio_out->pointlist[2 * i], trio_out->pointlist[2 * i + 1]});
        // myfile << (i + 1) << " " << trio_out->pointlist[2 * i] << " " << trio_out->pointlist[2 * i + 1] << " " << 0 << endl;
    }
    vector<pair<int, int>> e_out;
    // myfile.close();
    // myfile.open("test-triangles.txt");
    // myfile << trio_out->numberoftriangles << " " << 3 << " " << 0 << endl;
    for (int i = 0; i < trio_out->numberoftriangles; i++) {
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
        
        e_out.push_back(pair<int, int>{trio_out->trianglelist[3 * i] , trio_out->trianglelist[3 * i + 1]});
        e_out.push_back(pair<int, int>{trio_out->trianglelist[3 * i], trio_out->trianglelist[3 * i + 2]});
        e_out.push_back(pair<int, int>{trio_out->trianglelist[3 * i + 1], trio_out->trianglelist[3 * i + 2]});
        if (trio_out->trianglelist[3 * i] == -1 || trio_out->trianglelist[3 * i + 1] == -1 || trio_out->trianglelist[3 * i + 2] == -1) cout << "what is happening" << endl;
        // myfile << (i + 1) << " " << trio_out->trianglelist[3 * i] + 1 << " " << trio_out->trianglelist[3 * i + 1] + 1 << " " << trio_out->trianglelist[3 * i + 2] + 1 << endl;
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

vector<Bezier*> curves = {new Bezier{{{210, 210}, {-210, 210}}},
    new Bezier{{{-210, 210}, {-210, -210}}},
    new Bezier{{{-210, -210}, {210, -210}}},
    new Bezier{{{210, -210}, {210, 210}}}};
    
// vector<Bezier*> curves = {new Bezier{{{-210, 123}, {0, 56}, {210, 123}}},
//     new Bezier{{{210, 123}, {210, -67}}},
//     new Bezier{{{210, -67}, {0, 0}, {-210, -67}}},
//     new Bezier{{{-210, -67}, {-210, 123}}}};

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
double bridge_load = 2000.0;

void visualizeBoundary() {
    vector<glm::vec3> bPoints;
    vector<glm::vec3> pinnedPoints;
    for (int i = 0; i < curves.size(); i++) {
        vector<pair<double, double>> p = curves[i]->points;
        for (int j = 0; j < p.size(); j++) {
            glm::vec3 pnt{p[j].first, 0, p[j].second};
            bPoints.push_back(pnt);
            if (curves[i]->pinnedCurve) {
                pinnedPoints.push_back(pnt);
                pinned.push_back(p[j]);
            }
        }
    }
    polyscope::registerPointCloud("Boundary", bPoints);
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

double get_random() {
    static std::default_random_engine e;
    static std::uniform_real_distribution<> dist(1, 100);
    return dist(e);
}

void debugConstraints(string debugFileName, map<int, set<int>> neighbors, set<int> pinned, map<pair<int, int>, int> scalar_idx, 
                        vector<double> subs_scalars, vector<double> subs_z, vector<double> x_i, vector<double> y_i, int iter) {
    
    std::fstream debug(debugFileName, std::ios::out | std::ios::app);
    // debug.open(debugFileName);
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
    for (int i = 0; i < pinned.size(); i++) {
        pinned_idx.insert(pnt_idx_map[pinned[i]]);
    }
    cout << "Number of pinned vertices: " << pinned_idx.size() << endl;
    vector<glm::vec3> bridge_points;
    double load = bridge_load / (x_i.size());

    set<pair<int, int>>::iterator itr;
    map<pair<int, int>, int> scalar_idx;
    int i = 0;
    for (itr = edges.begin(); itr != edges.end(); itr++) {
        auto edge = *itr;
        scalar_idx[edge] = i++;
    }

    // Start with random edge scalars
    vector<double> subs_z(x_i.size(), 0.0);
    for (int i = 0; i < x_i.size(); i++) {
        subs_z[i] = get_random();
    }
    vector<double> subs_scalars(edges.size(), 0.0);
    for (int i = 0; i < edges.size(); i++) {
        subs_scalars[i] = get_random();
    }

    int num_nodes = (int)x_i.size();
    int num_edges = (int)edges.size();

    Eigen::MatrixXd scalarConstraints(3 * num_nodes, num_edges);
    Eigen::VectorXd bsolveScalar(3 * num_nodes);
    for (int i = 0; i < num_nodes; i++) {
        if (pinned_idx.find(i) == pinned_idx.end()) {
            set<int> neighbors_i = neighbors[i];
            set<int>::iterator itr;
            for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                int neighbor = *itr;
                pair<int, int> edge = make_edge(i, neighbor);
                int s_idx = scalar_idx[edge];
                scalarConstraints(i, s_idx) = x_i[i] - x_i[neighbor];
                scalarConstraints(i + 1, s_idx) = y_i[i] - y_i[neighbor];
                scalarConstraints(i + 2, s_idx) = subs_z[i] - subs_z[neighbor];
            }
            bsolveScalar(i) = 0;
            bsolveScalar(i + 1) = 0;
            bsolveScalar(i + 2) = load;
        }
    }
    auto scalars = scalarConstraints.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bsolveScalar);
    cout << scalarConstraints << endl;
    cout << bsolveScalar << endl;
    cout << scalars << endl;
    


    Eigen::MatrixXd nodeConstraints(num_nodes, num_nodes);
    for (int i = 0; i < x_i.size(); i++) {
        if (pinned_idx.find(i) == pinned_idx.end()) {
            set<int> neighbors_i = neighbors[i];
            set<int>::iterator itr;
            double z_i_scalar = 0.0;
            for (itr = neighbors_i.begin(); itr != neighbors_i.end(); itr++) {
                int neighbor = *itr;
                pair<int, int> edge = make_edge(i, neighbor);
                double s = subs_scalars[scalar_idx[edge]];
                nodeConstraints(i, neighbor) = -s;
                z_i_scalar += s;
            }
            nodeConstraints(i, i) = z_i_scalar;
        }
    }
    Eigen::VectorXd bsolveZ(num_nodes);
    for (int i =0 ; i < x_i.size(); i++) {
        if (pinned_idx.find(i) == pinned_idx.end()) {
            bsolveZ(i) = load;
        } else {
            bsolveZ(i) = 0;
        }
    }
    bsolveZ.setConstant(load);
    cout << nodeConstraints << endl;
    cout << bsolveZ << endl;
    auto z_values = nodeConstraints.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(bsolveZ);
    for (int i = 0; i < num_nodes; i++) {
        if (pinned_idx.find(i) != pinned_idx.end()) {
            subs_z[i] = 0.0;
        } else {
            subs_z[i] = z_values[i];
        }
    }
    cout << z_values << endl;

    debugConstraints("test.txt", neighbors, pinned_idx, scalar_idx, subs_scalars, subs_z, x_i, y_i, 0);
    for (int i = 0; i < x_i.size(); i++) {
        double x = x_i[i];
        double z = subs_z[i];
        double y = y_i[i];
        bridge_points.push_back(
            glm::vec3{x, z, y}
        );
    }

    vector<array<int, 2>> vec_edges;
    for (auto &p : edges) {
        array<int, 2> arr = {p.first, p.second};
        vec_edges.push_back(arr);
    }

    polyscope::registerCurveNetwork("my network", bridge_points, vec_edges);

    vector<glm::vec3> bridge_pnts_shadow;
    vector<array<int, 2>> shadow_edges;
    for (auto pnt: bridge_points) {
        bridge_pnts_shadow.push_back(pnt);
    }
    for (auto pnt: bridge_points) {
        bridge_pnts_shadow.push_back(glm::vec3{pnt[0], -20, pnt[2]});
    }
    for (int i = 0; i < bridge_points.size(); i++) {
        shadow_edges.push_back({i, i + (int)bridge_points.size()});
    }
    polyscope::registerCurveNetwork("shadow", bridge_pnts_shadow, shadow_edges);

}

void drawImGui() {
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
    ImGui::InputDouble("Total load on bridge", &bridge_load);
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