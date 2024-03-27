#include "polyscope/polyscope.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "polyscope/pick.h"
#include <unordered_map>

using namespace std;

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

struct Bezier {
  vector<array<double, 2>> controlPoints; // TODO: maybe change to pair<double, double>
  bool pinnedCurve;
  vector<pair<double, double>> points;
};

// No existing C++ default pair<T, T> hash
struct PairHash {
    template <typename T1, typename T2>
    size_t operator()(const pair<T1, T2>& p) const {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2; // Combining the hash values
    }
};

// Default values, 10 x 10 square with (0, 0) as leftmost corner
// vector<Bezier> curves = {Bezier{{{0, 0}, {0, 10}}}, 
//                          Bezier{{{0, 10}, {10, 10}}},
//                          Bezier{{{10, 10}, {10, 0}}},
//                          Bezier{{{10, 0}, {0, 0}}}};

// Pentagon
vector<Bezier> curves = {Bezier{{{0, 0}, {6, -3}, {12, 0}}}, 
  Bezier{{{12, 0}, {15, 9}}},
  Bezier{{{15, 9}, {10.5, 12}, {6, 15}}},
  Bezier{{{6, 15}, {-3, 9}}},
  Bezier{{{-3, 9}, {0, 0}}}};

unordered_map<pair<double, double>, vector<pair<double, double>>, PairHash> neighbors;

unordered_map<pair<double, double>, int, PairHash> point_idx_map; // Maps each point to its index in points (used to denote segments)

bool q = false;
bool D = false;
double maxArea = 0.0;

array<double, 2> deCastelJau(const vector<array<double, 2>>& controlPoints, double t) {
  if (controlPoints.size() == 1) {
    return controlPoints[0];
  }

  vector<array<double, 2>> intermediatePoints;
  for (int i = 0; i < controlPoints.size() - 1; i++) {
    double x = (1 - t) * controlPoints[i][0] + t * controlPoints[i + 1][0];
    double y = (1 - t) * controlPoints[i][1] + t * controlPoints[i + 1][1];
    intermediatePoints.push_back({x, y});
  }

  return deCastelJau(intermediatePoints, t);
}

vector<array<double, 2>> bezierCurve(int curveIdx) {
  vector<array<double, 2>> points;
  for (double t = 0; t <= 1.0; t += 0.1) {
    points.push_back(deCastelJau(curves[curveIdx].controlPoints, t));
  }
  return points;
}

void drawCurves() {
  vector<glm::vec3> vis_points;
  for (int i = 0; i < curves.size(); i++) {
    vector<array<double, 2>> bezier_points = bezierCurve(i);
    // TODO: Change to set the endpoints more generally, currently assumes EP are first and last control points
    bezier_points[0] = curves[i].controlPoints[0];
    bezier_points[bezier_points.size() - 1] = curves[i].controlPoints[curves[i].controlPoints.size() - 1];

    // Assign neighbors along each bezier curve s.t. all the curves can be stitched together
    pair<double, double> prev_point;
    for (int j = 0; j < bezier_points.size(); j++) {
      vis_points.push_back(glm::vec3{bezier_points[j][0], 0, bezier_points[j][1]});
      pair<double, double> p{bezier_points[j][0], bezier_points[j][1]};
      curves[i].points.push_back(p);
      if (j != 0) {
        neighbors[p].push_back(prev_point);
        neighbors[prev_point].push_back(p);
      }
      prev_point = p;
    }  
  }
  int p_index = 0;
  for (const auto& [v, n] : neighbors) {
    point_idx_map[v] = p_index;
    p_index++;
  }

  polyscope::registerPointCloud("curve", vis_points);
}

void callTriangulate() {
  // Want to construct consecutive segments of the bezier curves that are stitched together:
  // Choose an arbitrary start point, find an unvisited neighbor and make current point, repeat until we revisit the start node

  double points[neighbors.size() * 2]; // Input to `triangulate()`: flattened double array representing outer boundary points
                                       // neighbors.size() -> returns number of keys in hashmap, corresponds to number of unique points in boundary
  int p_index = 0;
  for (const auto& [v, n] : neighbors) {
    points[2 * p_index] = v.first;
    points[2 * p_index + 1] = v.second;
    p_index++;
  }
  
  vector<pair<double, double>> points_in_order; // Collects points in order of traversal

  int start_index = 0; // Arbitrarily starts at first point in hashmap
  pair<double, double> first_point{points[2 * start_index], points[2 * start_index + 1]}; // Used for termination condition
  pair<double, double> cur_point{points[2 * start_index], points[2 * start_index + 1]};
  pair<double, double> prev_point{points[2 * start_index], points[2 * start_index + 1]};
  
  while (true) {
    vector<pair<double, double>> p_neighbors = neighbors[cur_point]; // neighbors of current point (should be 2 neighbors)
    pair<double, double> n_1 = p_neighbors[0];
    pair<double, double> n_2 = p_neighbors[1];

    // We either got to the current point by n_1 or n_2, don't want to revisit that previous point
    if (n_1 != prev_point) {
      // The new segment is cur_point -> n_1
      points_in_order.push_back(cur_point);
      points_in_order.push_back(n_1);

      if (n_1 == first_point) break; // Looped back to starting point
      prev_point = cur_point;
      cur_point = n_1;
    } else if (n_2 != prev_point) {
      // The new segment is cur_point -> n_2
      points_in_order.push_back(cur_point);
      points_in_order.push_back(n_2);
      
      if (n_2 == first_point) break; // Looped back to starting point
      prev_point = cur_point;
      cur_point = n_2;
    }
  }

  vector<glm::vec3> boundary_points; // Used for visualizing boundary
  for (int i = 0; i < neighbors.size(); i++) {
    boundary_points.push_back(glm::vec3{points[2 * i], 0, points[2 * i + 1]});
  }

  vector<array<int, 2>> boundary_edges; // Used for visualizing boundary
  int segs[points_in_order.size()]; // Input to `triangulate()`: flattened int array
                                    // Every two values represents the indices of the endpoints of all segments
  for (int i = 0; i < points_in_order.size(); i += 2) {

    // TODO: seeing if zero-indexing causes different behavior
    segs[i] = point_idx_map[points_in_order[i]] + 1;
    segs[i + 1] = point_idx_map[points_in_order[i + 1]] + 1; 

    boundary_edges.push_back({segs[i] - 1, segs[i + 1] - 1});
  }

  auto boundary = polyscope::registerCurveNetwork("mesh boundary", boundary_points, boundary_edges);
  polyscope::registerPointCloud("boundary points", boundary_points);
  
  struct triangulateio* trio_in = new triangulateio;
  trio_in->pointlist = points;
  trio_in->numberofpoints = neighbors.size();
  trio_in->numberofpointattributes = 0;
  trio_in->pointmarkerlist = nullptr;
  trio_in->segmentlist = segs;
  trio_in->numberofsegments = points_in_order.size() / 2;
  // for (int i = 0; i < points_in_order.size(); i += 2) {
  //   cout << "seg " << trio_in->segmentlist[i] << " -> " << trio_in->pointlist[i + 1] << endl;
  // }
  trio_in->segmentmarkerlist = nullptr;
  trio_in->numberofholes = 0;
  trio_in->holelist = nullptr;
  trio_in->numberofregions = 0;
  trio_in->regionlist = nullptr;

  double* points_out;
  int* triangles;
  int* segments_out;
  int* segmentmarker_out;
  int* pointmarkerlist_out;

  struct triangulateio* trio_out = new triangulateio;
  trio_out->pointlist = nullptr;
  trio_out->trianglelist = nullptr;
  trio_out->segmentlist = nullptr;
  trio_out->segmentmarkerlist = nullptr;
  trio_out->pointmarkerlist = nullptr;

  // -p Triangulates a Planar Straight Line Graph
  // -z Numbers all items starting from zero (rather than one)
  // -a Imposes a maximum triangle area constraint
  string switches = "p";
  if (q) switches += "q";
  if (D) switches += "D";
  if (maxArea > 0.0) {
    switches += "a";
    switches += to_string(maxArea);
  }
  char* triswitches = const_cast<char*>(switches.c_str());
  
  triangulate(triswitches, trio_in, trio_out, nullptr);

  cout << "number of boundary edges " << trio_out->numberofsegments << endl;

  vector<glm::vec3> points_2d;
  vector<array<int, 2>> edges_2d;
  for (int i = 0; i < trio_out->numberofpoints; i++) {
    double x = trio_out->pointlist[2 * i];
    double y = trio_out->pointlist[2 * i + 1];
    points_2d.push_back(glm::vec3{x, 0, y});
  }
  for (int i = 0; i < trio_out->numberoftriangles; i++) {
    array<int, 2> seg1 = {trio_out->trianglelist[3 * i] - 1, trio_out->trianglelist[3 * i + 1] - 1};
    array<int, 2> seg2 = {trio_out->trianglelist[3 * i] - 1, trio_out->trianglelist[3 * i + 2] - 1};
    array<int, 2> seg3 = {trio_out->trianglelist[3 * i + 1] - 1, trio_out->trianglelist[3 * i + 2] - 1};
    edges_2d.push_back(seg1);
    edges_2d.push_back(seg2);
    edges_2d.push_back(seg3);
  }


  set<int> pinned_indices;
  vector<glm::vec3> pinned_points;
  for (int i = 0; i < curves.size(); i++) {
    if (curves[i].pinnedCurve) {
      for (int j = 0; j < curves[i].points.size(); j++) {
        pinned_indices.insert(point_idx_map[curves[i].points[j]]);
        pinned_points.push_back(glm::vec3{curves[i].points[j].first, 0, curves[i].points[j].second});
      }
    }
  }

  polyscope::registerPointCloud("pinned points", pinned_points);
  polyscope::registerCurveNetwork("triangulated mesh", points_2d, edges_2d);
}

void drawImGui() {

  // TODO: use mouse clicks to pin/un-pin vertices
  ImGui::ShowDemoWindow();
  if (ImGui::Button("Add Curve")) {
    curves.push_back(Bezier());
  }

  for (int i = 0; i < curves.size(); i++) {
    // Begin a collapsible header for each curve
    if (ImGui::CollapsingHeader(("Curve " + std::to_string(i)).c_str())) {
        if (ImGui::Checkbox(("Pin edge##_" + to_string(i)).c_str(), &curves[i].pinnedCurve)) {
            cout << "pinned " << i << endl;
        }

        if (ImGui::Button(("Delete##_" + to_string(i)).c_str())) {
            curves.erase(curves.begin() + i);
            continue;
        }

        if (ImGui::Button(("Add Control Point##_" + to_string(i)).c_str())) {
            curves[i].controlPoints.push_back({0.0, 0.0});
        }

        for (int j = 0; j < curves[i].controlPoints.size(); j++) {
            ImGui::Text("Control Point %d", j);
            ImGui::InputDouble((std::to_string(i) + "_" + std::to_string(j)+ "_x").c_str(), &curves[i].controlPoints[j][0]);
            ImGui::InputDouble((std::to_string(i) + "_" + std::to_string(j) + "_y").c_str(), &curves[i].controlPoints[j][1]);
        }
    }
  }

  if (ImGui::Button("Draw Curves")) {
    drawCurves();
  }

  if (ImGui::CollapsingHeader("Triangulation Parameters")) {
    if (ImGui::Checkbox("No angles smaller than 20 degrees", &q)) {
      // switches += "q";
    }
    if (ImGui::Checkbox("Conforming Delaunay", &D)) {
      // switches += "D";
    }
    // if (ImGui::Checkbox("Constrain max triangle area", &enableMax)) {
    ImGui::InputDouble("Maximum triangle area", &maxArea);
    // }
  }

  if (ImGui::Button("Triangulate")) {
    callTriangulate();
  }
}

int main() {
  polyscope::init();

  polyscope::state::userCallback = drawImGui;

  polyscope::show();
}