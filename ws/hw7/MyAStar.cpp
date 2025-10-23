/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- ChatGPT misc debugging
- ChatGPT for psuedocode and starter code
*/

#include "MyAStar.h"
#include <queue>
#include <unordered_map>
#include <limits>
#include <functional>

struct NodeCost {
    amp::Node node;
    double f_cost;  // f = g + h
    double g_cost;  // cost from start to this node

    bool operator>(const NodeCost& other) const {
        return f_cost > other.f_cost; // for min-heap
    }
};


// Implement the search method for the A* algorithm
MyAStarAlgo::GraphSearchResult MyAStarAlgo::search(const amp::ShortestPathProblem& problem, const amp::SearchHeuristic& heuristic) {
    std::cout << "Starting A* Graph Search: Init --> goal | " << problem.init_node << " --> " << problem.goal_node << std::endl;
    GraphSearchResult result = {false, {}, 0.0}; // initialize the results object

    // Min-heap priority queue
    std::priority_queue<NodeCost, std::vector<NodeCost>, std::greater<NodeCost>> open_set;

    // Maps nodes to their cheapest g_cost so far
    std::unordered_map<amp::Node, double> g_cost_map;
    g_cost_map[problem.init_node] = 0.0;

    // Map to reconstruct path
    std::unordered_map<amp::Node, amp::Node> came_from;

    // Start node
    open_set.push({problem.init_node, heuristic(problem.init_node), 0.0});

    // 
    int iter = 1;
    while (!open_set.empty()) {
        //std::cout << iter << std::endl;
        NodeCost current = open_set.top();
        open_set.pop();

        // Goal check
        if (current.node == problem.goal_node) {
            result.success = true;
            result.path_cost = current.g_cost;

            // Reconstruct path
            amp::Node n = current.node;
            while (came_from.find(n) != came_from.end()) {
                result.node_path.insert(result.node_path.begin(), n);
                n = came_from[n];
            }
            result.node_path.insert(result.node_path.begin(), problem.init_node);

            result.print();
            return result;
        }

        // Explore neighbors
        const auto& neighbors = problem.graph->children(current.node);
        const auto& edges = problem.graph->outgoingEdges(current.node);

        for (size_t i = 0; i < neighbors.size(); ++i) {
            amp::Node neighbor = neighbors[i];
            double edge_cost = edges[i]; // assumes outgoingEdges() returns weights as doubles

            double tentative_g = current.g_cost + edge_cost;
            if (g_cost_map.find(neighbor) == g_cost_map.end() || tentative_g < g_cost_map[neighbor]) {
                came_from[neighbor] = current.node;
                g_cost_map[neighbor] = tentative_g;

                double f = tentative_g + heuristic(neighbor);
                open_set.push({neighbor, f, tentative_g});
            }
        }

        iter = iter + 1;
    }

    result.print();
    return result;
}
