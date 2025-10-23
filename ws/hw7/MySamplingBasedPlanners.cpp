/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- ChatGPT misc debugging
- Lecture slides
- ChatGPT for random sampling and nearest neighbors
- ChatGPT for general assistance with C++
*/


# include "MySamplingBasedPlanners.h"

#include <iostream>
#include <memory>
#include <vector>
#include <map>
#include <tuple>
#include <random>
#include <algorithm>
#include <Eigen/Dense>
#include "MyAStar.h"



/*
A point in Polygon algorithm that returns true if the point is in the polygon
*/
bool pip(std::vector<Eigen::Vector2d> vertices, Eigen::Vector2d x) {
    // Initialize a count variable
    int count = 0;

    // Happy direction vector!
    Eigen::Vector2d d(1.0, 0.0);

    // Loop through the verticies
    for (size_t j = 0; j < vertices.size(); ++j) {
        Eigen::Vector2d p1 = vertices[j];
        Eigen::Vector2d p2 = vertices[(j + 1) % vertices.size()];
        Eigen::Vector2d edge = p2 - p1;

        Eigen::Matrix2d A;
        A << d, -edge;

        double det = A.determinant();
        if (std::abs(det) < 1e-12) continue; // Parallel ray/edge

        Eigen::Vector2d sol = A.inverse() * (p1 - x);
        double t = sol[0], u = sol[1];

        if (t >= 0 && u >= 0 && u <= 1) count++;
    }

    return (count % 2 == 1);
}


bool segmentIntersect(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                      const Eigen::Vector2d& q1, const Eigen::Vector2d& q2)
{
    Eigen::Vector2d r = p2 - p1;
    Eigen::Vector2d s = q2 - q1;

    double rxs = r.x() * s.y() - r.y() * s.x();
    double qpxr = (q1 - p1).x() * r.y() - (q1 - p1).y() * r.x();

    if (rxs == 0 && qpxr == 0) return true; // collinear
    if (rxs == 0 && qpxr != 0) return false; // parallel

    double t = ((q1 - p1).x() * s.y() - (q1 - p1).y() * s.x()) / rxs;
    double u = ((q1 - p1).x() * r.y() - (q1 - p1).y() * r.x()) / rxs;

    return (t >= 0 && t <= 1 && u >= 0 && u <= 1);
}


bool lineIntersectsPolygon(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2, const std::vector<Eigen::Vector2d>& poly)
{
    int N = poly.size();
    for (int i = 0; i < N; ++i)
    {
        Eigen::Vector2d q1 = poly[i];
        Eigen::Vector2d q2 = poly[(i+1) % N];  // wrap around
        if (segmentIntersect(p1, p2, q1, q2))
            return true;
    }
    return false;
}

struct EuclideanHeuristic2D : public amp::SearchHeuristic {
    // Reference to node positions
    const std::map<amp::Node, Eigen::Vector2d>& nodes;

    // Goal node coordinates
    Eigen::Vector2d goal;

    // Constructor
    EuclideanHeuristic2D(const std::map<amp::Node, Eigen::Vector2d>& nodes_, const Eigen::Vector2d& goal_)
        : nodes(nodes_), goal(goal_) {}

    // Override the virtual operator()
    virtual double operator()(amp::Node node) const override {
        auto it = nodes.find(node);
        if (it == nodes.end()) {
            throw std::runtime_error("Node not found in heuristic!");
        }
        const Eigen::Vector2d& pos = it->second;

        // Euclidean distance to goal
        return (pos - goal).norm();
    }
};

// Implement your PRM algorithm here
amp::Path2D MyPRM::plan(const amp::Problem2D& problem) {
    // Random number generation 
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<double> distX(problem.x_min, problem.x_max); 
    std::uniform_real_distribution<double> distY(problem.y_min, problem.y_max); 
    // std::uniform_real_distribution<double> distX(-1.0, 11.0); 
    // std::uniform_real_distribution<double> distY(-3.0, 3.0); 

    // Sample points 
    Eigen::Vector2d q(0.0,0.0); 
    for (amp::Node i = 0; i < N; ++i) { 
        // Get a sample point 
        q = Eigen::Vector2d(distX(gen), distY(gen)); 

        // Determine if valid 
        bool collision = false; 
        for (const auto& obs : problem.obstacles) { 
            if (pip(obs.verticesCCW(), q)) { 
                collision = true; break; 
            } 
        } 

        // If not in collision, add to nodes 
        if (!collision) {
             nodes[nodes.size()+1] = q; 
        } 
    } 

    // Add starting and ending nodes 
    amp::Node start_idx = nodes.size();
    nodes[start_idx] = problem.q_init;

    amp::Node goal_idx = nodes.size();
    ++goal_idx; 
    nodes[goal_idx] = problem.q_goal;
        
    // Connect nodes f
    for (amp::Node i = 0; i < nodes.size(); ++i) {
        for (amp::Node j = i + 1; j < nodes.size(); ++j) { // i+1 avoids duplicate edges 
            double d = (nodes[i] - nodes[j]).norm(); if (d <= r) { 

            // Determine if valid 
            bool collision = false; 
            for (const auto& obs : problem.obstacles) {
                if (lineIntersectsPolygon(nodes[i], nodes[j], obs.verticesCCW()))
                {
                    collision = true;
                    break;
                }
            }
             // If not in collision, add connection 
            if (!collision) { 
                graphPtr->connect(i, j, d); 
                graphPtr->connect(j, i, d);
            } 
            
            } 
        }
    } 


    // Output the path
    amp::Path2D path;

    auto graph_nodes = graphPtr->nodes();
    bool start_exists = std::find(graph_nodes.begin(), graph_nodes.end(), start_idx) != graph_nodes.end();
    bool goal_exists = std::find(graph_nodes.begin(), graph_nodes.end(), goal_idx) != graph_nodes.end();

    if (!start_exists || !goal_exists) {
        amp::Path2D path;
        path.waypoints.push_back(problem.q_init);
        path.waypoints.push_back(problem.q_goal);
        return path;
    }
                
    // Use A* to find the shortest 
    amp::ShortestPathProblem sp = amp::ShortestPathProblem(graphPtr, start_idx, goal_idx); 
    EuclideanHeuristic2D heuristic(nodes, nodes[goal_idx]); 
    MyAStarAlgo algo; 
    MyAStarAlgo::GraphSearchResult result = algo.search(sp, heuristic);


    if (!result.success)
    {
        path.waypoints.push_back(problem.q_init);
        path.waypoints.push_back(problem.q_goal);
        return path;
    } else {
        plan_success = true;
    }
    
    if (!smooth)
    {
        for (const amp::Node& node : result.node_path) {
            path.waypoints.push_back(nodes[node]);
        }
    } else { 
        std::random_device rd;
        std::mt19937 gen(rd());
        const size_t N = result.node_path.size();

        std::vector<amp::Node> path_vec(result.node_path.begin(), result.node_path.end());

        for (size_t k = 0; k < 10; ++k) {
            if (path_vec.size() < 2) break;

            std::uniform_int_distribution<size_t> dist(0, path_vec.size() - 1);
            size_t i = dist(gen);
            size_t j = dist(gen);
            if (i == j) continue;
            if (i > j) std::swap(i, j);

            // collision check
            bool collision = false;
            for (const auto& obs : problem.obstacles) {
                if (lineIntersectsPolygon(nodes[path_vec[i]], nodes[path_vec[j]], obs.verticesCCW())) {
                    collision = true;
                    break;
                }
            }

            if (!collision) {
                std::vector<amp::Node> new_path;
                new_path.reserve(path_vec.size() - (j - i - 1)); // reserve estimated size

                // copy up to and including i
                new_path.insert(new_path.end(), path_vec.begin(), path_vec.begin() + (i + 1));

                // add j
                new_path.push_back(path_vec[j]);

                // copy remaining nodes after j (if any)
                if (j + 1 < path_vec.size()) {
                    new_path.insert(new_path.end(), path_vec.begin() + (j + 1), path_vec.end());
                }

                // swap into path_vec
                path_vec.swap(new_path);

                // restart attempts so we re-check shortcuts on the shortened path
                k = 0;

            }
        }

        // Rebuild list
        result.node_path.assign(path_vec.begin(), path_vec.end());
        
    }

    // Push waypoints
    for (const amp::Node& node : result.node_path) {
        path.waypoints.push_back(nodes[node]);
    }

    if (path.waypoints.size() == 0)
    {
        path.waypoints.push_back(problem.q_init);
        path.waypoints.push_back(problem.q_goal);
    }
    

    
    return path;
}














// Implement your RRT algorithm here
amp::Path2D MyRRT::plan(const amp::Problem2D& problem) {
    // Random number generation 
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<double> distX(problem.x_min, problem.x_max); 
    std::uniform_real_distribution<double> distY(problem.y_min, problem.y_max);
    std::uniform_real_distribution<double> dis(0.0, 1.0);


    // Create the base node 
    int node_idx = 0; int goal_idx = 0;
    nodes[node_idx] = problem.q_init;

    // Loop until soultion found 
    int safety = 0; Eigen::Vector2d q(0.0,0.0);
    while (true)
    {

        // Get a random sample with goal biasing
        if (dis(gen) < p_goal)
        {
            q = problem.q_goal;
        } else {
            q = Eigen::Vector2d(distX(gen), distY(gen)); 
        }

        // Determine closest node in the tree
        int i_close = 0; float dist_close = 1000.0;
        float dist;
        for (size_t i = 0; i < nodes.size(); i++)
        {
            dist = (nodes[i] - q).norm(); 

            if (dist < dist_close)
            {
                dist_close = dist; 
                i_close = i;
            }
        }

        // Take a step towars that point
        Eigen::Vector2d dir = q - nodes[i_close];
        if (dir.norm() > r) dir.normalize();
        Eigen::Vector2d q_new = nodes[i_close] + dir * std::min(static_cast<double>(r), dir.norm());
        
        // Determine if collison free
        bool collision = false; 
        for (const auto& obs : problem.obstacles) {
            if (lineIntersectsPolygon(nodes[i_close], q_new, obs.verticesCCW()))
            {
                collision = true;
                break;
            }
        }

        // Add to path if valid 
        if (!collision)
        {

            // Create the new node and add it to the graph
            node_idx = node_idx + 1;
            nodes[node_idx] = q_new;
            graphPtr->connect(i_close, node_idx, r); 
            graphPtr->connect(node_idx, i_close, r); 

            // See if a solution is found 
            if ((q_new - problem.q_goal).norm() <= epsilon)
            {
                // Add goal node
                node_idx = node_idx + 1;
                node_idx = node_idx;
                nodes[node_idx] = problem.q_goal;
                graphPtr->connect(node_idx-1, node_idx, r); 
                graphPtr->connect(node_idx, node_idx-1, r); 
                goal_idx = node_idx;

                break;
            } 
        }
    
        // Safety mechanism
        safety = safety + 1;
        if (safety >= 10000) //10000
        {
            std::cout << "Safety" << std::endl;
            break;
        } 
    }

    // Use A* to find the shortest 
    amp::ShortestPathProblem sp = amp::ShortestPathProblem(graphPtr, 0, goal_idx); 
    EuclideanHeuristic2D heuristic(nodes, nodes[goal_idx]); 
    MyAStarAlgo algo; 
    MyAStarAlgo::GraphSearchResult result = algo.search(sp, heuristic);

    // Output the path
    amp::Path2D path;
    if (!result.success)
    {
        path.waypoints.push_back(problem.q_init);
        path.waypoints.push_back(problem.q_goal);
        return path;
    } else {
        plan_success = true;
    }


    for (const amp::Node& node : result.node_path) {
        path.waypoints.push_back(nodes[node]);
    }

    return path;
}