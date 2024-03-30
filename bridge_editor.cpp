#include <unordered_map>
#include <iostream>

#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"

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
  for (double t = 0; t <= 1.0; t += 0.1) {
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
    }
    points_out = p_out;
    edges_out = e_out;
}

vector<Bezier*> curves = {new Bezier{{{0, 0}, {6, -3}, {12, 0}}}, 
  new Bezier{{{12, 0}, {15, 9}}},
  new Bezier{{{15, 9}, {10.5, 12}, {6, 15}}},
  new Bezier{{{6, 15}, {-3, 9}}},
  new Bezier{{{-3, 9}, {0, 0}}}};
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

void visualizeBoundary() {
    vector<glm::vec3> bPoints;
    vector<glm::vec3> pinnedPoints;
    for (int i = 0; i < curves.size(); i++) {
        vector<pair<double, double>> p = curves[i]->points;
        for (int j = 0; j < p.size(); j++) {
            glm::vec3 pnt{p[j].first, 0, p[j].second};
            bPoints.push_back(pnt);
            if (curves[i]->pinnedCurve) pinnedPoints.push_back(pnt);
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
    
    if (ImGui::Button("Construct Boundary")) {
        constructBoundary(curves, points, segments, in_size);
        curvesReady = true;
        visualizeBoundary();
    }
    
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
        }
    }
}

int main() {
  polyscope::init();

  polyscope::state::userCallback = drawImGui;

  polyscope::show();
}