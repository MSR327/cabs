#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <limits>
#include <cmath>
#include <string>
#include <sstream>
#include <algorithm>
#include "json.hpp" // Download from https://github.com/nlohmann/json

void printGrid(const vector<Node>& nodes, int grid_size, int userNodeId, const vector<Cab>& cabs) {
    vector<string> symbols(nodes.size(), ".");
    
    symbols[userNodeId] = "U";
    for (const auto& cab : cabs) {
        if (symbols[cab.node_id] == "U")
            symbols[cab.node_id] = "U/C"; // User and cab at same node
        else if (symbols[cab.node_id] == ".")
            symbols[cab.node_id] = "C";
        else if (symbols[cab.node_id] == "C")
            symbols[cab.node_id] = "C"; // Already marked as cab
    }

    cout << "\n===== ASCII Map of Grid =====\n";
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            int idx = i * grid_size + j;
            cout.width(4);
            cout << symbols[idx];
        }
        cout << endl;
    }
    cout << "=============================\n\n";
}

using json = nlohmann::json;
using namespace std;

struct Node {
    double lat, lon;
    vector<pair<int, double>> neighbors; // (neighbor index, distance)
};

struct Cab {
    int user_id;
    string name;
    double lat, lon;
    int node_id;
    double route_distance_dijkstra;
    double route_distance_astar;
};

// Haversine distance in meters
inline double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double R = 6371000.0;
    double dLat = (lat2 - lat1) * M_PI / 180.0;
    double dLon = (lon2 - lon1) * M_PI / 180.0;
    double a = pow(sin(dLat/2),2) + cos(lat1*M_PI/180.0)*cos(lat2*M_PI/180.0)*pow(sin(dLon/2),2);
    double c = 2*atan2(sqrt(a), sqrt(1-a));
    return R * c;
}

int mapToNode(const vector<Node>& nodes, double lat, double lon) {
    int best = 0;
    double best_dist = numeric_limits<double>::max();
    for (int i = 0; i < nodes.size(); ++i) {
        double d = haversine(lat, lon, nodes[i].lat, nodes[i].lon);
        if (d < best_dist) {
            best_dist = d;
            best = i;
        }
    }
    return best;
}

// Dijkstra's algorithm
vector<double> dijkstra(const vector<Node>& nodes, int start) {
    int n = nodes.size();
    vector<double> dist(n, numeric_limits<double>::infinity());
    dist[start] = 0;
    using P = pair<double, int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    pq.push({0, start});
    while (!pq.empty()) {
        auto [d, u] = pq.top(); pq.pop();
        if (d > dist[u]) continue;
        for (const auto& [v, w] : nodes[u].neighbors) {
            if (dist[u] + w < dist[v]) {
                dist[v] = dist[u] + w;
                pq.push({dist[v], v});
            }
        }
    }
    return dist;
}

// A* algorithm
struct AStarNode {
    int id;
    double g, h, f;
    int parent;
    bool operator>(const AStarNode& other) const {
        return f > other.f;
    }
};

double heuristic(const Node& a, const Node& b) {
    // Euclidean distance in meters (Haversine)
    return haversine(a.lat, a.lon, b.lat, b.lon);
}

vector<double> astar(const vector<Node>& nodes, int start, int goal) {
    int n = nodes.size();
    vector<double> dist(n, numeric_limits<double>::infinity());
    vector<int> parent(n, -1);
    vector<bool> closed(n, false);
    dist[start] = 0;

    priority_queue<AStarNode, vector<AStarNode>, greater<AStarNode>> open;
    open.push({start, 0, heuristic(nodes[start], nodes[goal]), heuristic(nodes[start], nodes[goal]), -1});

    while (!open.empty()) {
        AStarNode current = open.top();
        open.pop();

        if (closed[current.id]) continue;
        closed[current.id] = true;

        if (current.id == goal) {
            // reconstruct distances
            vector<double> result(n, numeric_limits<double>::infinity());
            int cur = goal;
            while (cur != -1) {
                result[cur] = dist[cur];
                cur = parent[cur];
            }
            return result;
        }

        for (const auto& [neighbor, weight] : nodes[current.id].neighbors) {
            if (closed[neighbor]) continue;
            double tentative_g = dist[current.id] + weight;
            if (tentative_g < dist[neighbor]) {
                dist[neighbor] = tentative_g;
                parent[neighbor] = current.id;
                double h = heuristic(nodes[neighbor], nodes[goal]);
                open.push({neighbor, tentative_g, h, tentative_g + h, current.id});
            }
        }
    }
    return vector<double>(n, numeric_limits<double>::infinity()); // no path
}

int main() {
    // Build a simple 5x5 grid network
    vector<Node> nodes;
    double base_lat = 12.9716, base_lon = 77.5946;
    int grid_size = 5;
    double step = 0.009;
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            nodes.push_back({base_lat + i * step, base_lon + j * step, {}});
        }
    }
    for (int i = 0; i < grid_size; ++i) {
        for (int j = 0; j < grid_size; ++j) {
            int idx = i * grid_size + j;
            if (i > 0) {
                int up = (i-1) * grid_size + j;
                double d = haversine(nodes[idx].lat, nodes[idx].lon, nodes[up].lat, nodes[up].lon);
                nodes[idx].neighbors.push_back({up, d});
            }
            if (i < grid_size-1) {
                int down = (i+1) * grid_size + j;
                double d = haversine(nodes[idx].lat, nodes[idx].lon, nodes[down].lat, nodes[down].lon);
                nodes[idx].neighbors.push_back({down, d});
            }
            if (j > 0) {
                int left = i * grid_size + (j-1);
                double d = haversine(nodes[idx].lat, nodes[idx].lon, nodes[left].lat, nodes[left].lon);
                nodes[idx].neighbors.push_back({left, d});
            }
            if (j < grid_size-1) {
                int right = i * grid_size + (j+1);
                double d = haversine(nodes[idx].lat, nodes[idx].lon, nodes[right].lat, nodes[right].lon);
                nodes[idx].neighbors.push_back({right, d});
            }
        }
    }

    // Read JSON input
    cout << "Enter path to JSON file (e.g., cabs.json): ";
    string json_file;
    cin >> json_file;
    ifstream file(json_file);
    if (!file) {
        cout << "File not found!" << endl;
        return 1;
    }
    json data;
    file >> data;

    double user_lat = stod(data[0]["latitude"].get<string>());
    double user_lon = stod(data[0]["longitude"].get<string>());
    int userNodeId = mapToNode(nodes, user_lat, user_lon);

    vector<Cab> cabs;
    for (size_t i = 1; i < data.size(); ++i) {
        Cab cab;
        cab.user_id = data[i]["user_id"];
        cab.name = data[i]["name"];
        cab.lat = stod(data[i]["latitude"].get<string>());
        cab.lon = stod(data[i]["longitude"].get<string>());
        cab.node_id = mapToNode(nodes, cab.lat, cab.lon);
        cabs.push_back(cab);
    }

    // Run Dijkstra
    auto dist_dijkstra = dijkstra(nodes, userNodeId);

    // For each cab, calculate both Dijkstra and A* distances
    for (auto& cab : cabs) {
        cab.route_distance_dijkstra = dist_dijkstra[cab.node_id];
        auto dist_astar = astar(nodes, userNodeId, cab.node_id);
        cab.route_distance_astar = dist_astar[cab.node_id];
    }

    // Filter cabs within 50 km for both algorithms
    cout << "Nearby cabs within 50 km route distance (Dijkstra and A*):\n";
    int count = 0;
    for (const auto& cab : cabs) {
        if (cab.route_distance_dijkstra <= 50000 && cab.route_distance_astar <= 50000) {
            cout << "Cab Name: " << cab.name
                 << ", User ID: " << cab.user_id
                 << ", Dijkstra Distance: " << cab.route_distance_dijkstra / 1000 << " km"
                 << ", A* Distance: " << cab.route_distance_astar / 1000 << " km" << endl;
            count++;
        }
    }
    cout << "Total nearby cabs found: " << count << endl;

    return 0;
}
