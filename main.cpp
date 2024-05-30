#include <iostream>
#include <vector>
#include <queue>
#include <algorithm>
#include <cmath>
#include <limits>
#include <chrono>
#include <fstream>
#include <random>
#include <string>

using namespace std;

struct City {
    int id;
    string name;
    double x, y;
    vector<pair<int, double>> neighbors;
};

struct Node {
    vector<int> path;
    double cost;
    double max_distance;
};

// Custom comparator for priority queue
struct NodeComparator {
    bool operator()(const Node& a, const Node& b) const {
        return a.max_distance > b.max_distance;
    }
};

bool AllCitiesVisited(const vector<int>& path, size_t num_cities) {
    return path.size() == num_cities;
}

double CalculateRouteCost(const vector<int>& path, const vector<City>& cities) {
    double cost = 0;
    double max_distance = 0;
    for (size_t i = 0; i < path.size() - 1; ++i) {
        int from = path[i];
        int to = path[i + 1];
        double distance = 0;
        for (const auto& neighbor : cities[from].neighbors) {
            if (neighbor.first == to) {
                distance = neighbor.second;
                break;
            }
        }
        cost += distance;
        max_distance = max(max_distance, distance);
    }
    return max_distance;
}

void DFS(const vector<City>& cities, int start, vector<int>& path, double& best_cost, vector<int>& best_path, double& max_distance) {
    if (AllCitiesVisited(path, cities.size())) {
        path.push_back(start);  // Return to the start city
        double current_max_distance = CalculateRouteCost(path, cities);
        if (current_max_distance < max_distance) {
            max_distance = current_max_distance;
            best_path = path;
        }
        path.pop_back();
        return;
    }

    for (const auto& neighbor : cities[path.back()].neighbors) {
        if (find(path.begin(), path.end(), neighbor.first) == path.end()) {
            path.push_back(neighbor.first);
            DFS(cities, start, path, best_cost, best_path, max_distance);
            path.pop_back();
        }
    }
}

void LeastCostSearch(const vector<City>& cities, int start) {
    priority_queue<Node, vector<Node>, NodeComparator> pq;
    pq.push({ {start}, 0, 0 });

    double best_cost = numeric_limits<double>::infinity();
    vector<int> best_path;
    double max_distance = numeric_limits<double>::infinity();

    while (!pq.empty()) {
        Node current = pq.top();
        pq.pop();

        if (AllCitiesVisited(current.path, cities.size())) {
            current.path.push_back(start);  // Return to the start city
            double current_max_distance = CalculateRouteCost(current.path, cities);
            if (current_max_distance < max_distance) {
                max_distance = current_max_distance;
                best_path = current.path;
            }
            continue;
        }

        for (const auto& neighbor : cities[current.path.back()].neighbors) {
            if (find(current.path.begin(), current.path.end(), neighbor.first) == current.path.end()) {
                Node next = current;
                next.path.push_back(neighbor.first);
                next.cost += neighbor.second;
                next.max_distance = max(next.max_distance, neighbor.second);
                pq.push(next);
            }
        }
    }

    cout << "\tLeast-Cost Search: \n\tMax distance = " << max_distance << "\n";
    cout << "\tPath: ";
    for (int city : best_path) {
        cout << city << " ";
    }
    cout << "\n";
}

double Heuristic(const Node& node, const std::vector<City>& cities, int start) {
    // Calculate the Euclidean distance to the start city
    int current_city = node.path.back();
    double dx = cities[current_city].x - cities[start].x;
    double dy = cities[current_city].y - cities[start].y;
    return std::sqrt(dx * dx + dy * dy);
}

void AStarSearch(const vector<City>& cities, int start) {
    priority_queue<Node, vector<Node>, NodeComparator> pq;
    pq.push({ {start}, 0, 0 });

    double best_cost = numeric_limits<double>::infinity();
    vector<int> best_path;
    double max_distance = numeric_limits<double>::infinity();

    while (!pq.empty()) {
        Node current = pq.top();
        pq.pop();

        if (AllCitiesVisited(current.path, cities.size())) {
            current.path.push_back(start);  // Return to the start city
            double current_max_distance = CalculateRouteCost(current.path, cities);
            if (current_max_distance < max_distance) {
                max_distance = current_max_distance;
                best_path = current.path;
            }
            continue;
        }

        for (const auto& neighbor : cities[current.path.back()].neighbors) {
            if (find(current.path.begin(), current.path.end(), neighbor.first) == current.path.end()) {
                Node next = current;
                next.path.push_back(neighbor.first);
                next.cost += neighbor.second;
                next.max_distance = max(next.max_distance, neighbor.second);
                double priority = next.cost + Heuristic(next, cities, start);
                pq.push({ next.path, next.cost, next.max_distance });
            }
        }
    }

    cout << "\tA* Search: \n\tMax distance = " << max_distance << "\n";
    cout << "\tPath: ";
    for (int city : best_path) {
        cout << city << " ";
    }
    cout << "\n";
}

vector<City> generateRandomCities(int num_cities) {
    vector<City> cities;

    // Random number generator
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> dis(0.0, 100.0); // coordinates range from 0 to 100

    // Generate cities
    cout << "Generated Cities:\n";
    for (int i = 0; i < num_cities; ++i) {
        string name = "City_" + to_string(i);
        double x = dis(gen);
        double y = dis(gen);
        cities.push_back({ i, name, x, y, {} });
        cout << name << ": (" << x << ", " << y << ")\n";
    }

    // Calculate distances between cities
    for (int i = 0; i < num_cities; ++i) {
        for (int j = i + 1; j < num_cities; ++j) {
            double dx = cities[i].x - cities[j].x;
            double dy = cities[i].y - cities[j].y;
            double distance = sqrt(dx * dx + dy * dy);
            cities[i].neighbors.push_back({ j, distance });
            cities[j].neighbors.push_back({ i, distance });
        }
    }

    return cities;
}



int main() {



    vector<City> cities = generateRandomCities(10); // Generate 10 random cities

    int start = 0;
    cout << "-------------------------------------------------------\n";
    cout << "\t(This may take a while...)\n";
    cout << "-------------------------------------------------------\n";

    // DFS
    auto start_time = chrono::high_resolution_clock::now();
    vector<int> path = { start };
    double best_cost = numeric_limits<double>::infinity();
    vector<int> best_path;
    double max_distance = numeric_limits<double>::infinity();
    DFS(cities, start, path, best_cost, best_path, max_distance);
    auto end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end_time - start_time;
    cout << "\tDFS: \n\tMax distance = " << max_distance << "\n";
    cout << "\tPath: ";
    for (int city : best_path) {
        cout << city << " ";
    }
    cout << "\n";
    cout << "\tDFS took " << duration.count() << " seconds\n";
    cout << "-------------------------------------------------------\n";


    // Least-Cost Search
    start_time = chrono::high_resolution_clock::now();
    LeastCostSearch(cities, start);
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    cout << "\tLeast-Cost Search took " << duration.count() << " seconds\n";
    cout << "-------------------------------------------------------\n";

    // A* Search
    start_time = chrono::high_resolution_clock::now();
    AStarSearch(cities, start);
    end_time = chrono::high_resolution_clock::now();
    duration = end_time - start_time;
    cout << "\tA* Search took " << duration.count() << " seconds\n";
    cout << "-------------------------------------------------------\n";

    return 0;
}