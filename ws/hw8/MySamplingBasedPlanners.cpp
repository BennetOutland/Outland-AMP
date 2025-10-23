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



// Compute squared distance from point p to segment [a,b]
double pointToSegmentDist2(const Eigen::Vector2d& p,
                           const Eigen::Vector2d& a,
                           const Eigen::Vector2d& b)
{
    Eigen::Vector2d ab = b - a;
    double t = (p - a).dot(ab) / ab.squaredNorm();
    t = std::clamp(t, 0.0, 1.0);
    Eigen::Vector2d proj = a + t * ab;
    return (p - proj).squaredNorm();
}

// Compute minimum distance squared between two segments [p1,p2] and [q1,q2]
double segmentToSegmentDist2(const Eigen::Vector2d& p1, const Eigen::Vector2d& p2,
                             const Eigen::Vector2d& q1, const Eigen::Vector2d& q2)
{
    // Brute-force check both endpoints projected onto the other segment
    double d2 = std::min({
        pointToSegmentDist2(p1, q1, q2),
        pointToSegmentDist2(p2, q1, q2),
        pointToSegmentDist2(q1, p1, p2),
        pointToSegmentDist2(q2, p1, p2)
    });
    return d2;
}


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

// Check if a disk robot moving along [p1,p2] collides with polygon verticesCCW
bool diskPathIntersectsPolygon(const Eigen::Vector2d& p1,
                               const Eigen::Vector2d& p2,
                               const std::vector<Eigen::Vector2d>& vertices,
                               double robot_radius)
{
    int n = vertices.size();
    double r2 = robot_radius * robot_radius;

    for (int i = 0; i < n; i++) {
        const Eigen::Vector2d& q1 = vertices[i];
        const Eigen::Vector2d& q2 = vertices[(i+1) % n];

        if (segmentToSegmentDist2(p1, p2, q1, q2) <= r2)
            return true;
    }

    // Also check if start or end point is *inside* the polygon (robot starts within obstacle)
    if (pip(vertices, p1) || pip(vertices, p2))
        return true;

    return false;
}



bool segmentAgentCollision(const Eigen::Vector2d& q1,
                           const Eigen::Vector2d& q2,
                           double r1,
                           const std::vector<amp::Path2D>& agent_paths,
                           double step_size)
{
    int samples = std::max(2, (int)std::ceil((q2 - q1).norm() / step_size));
    for (int s = 0; s <= samples; ++s) {
        double tau = double(s) / samples;
        Eigen::Vector2d p = q1 + tau * (q2 - q1);
        double cost_here = tau * (q2 - q1).norm();  // local cost along segment

        // Check all agents’ interpolated positions
        for (size_t id = 0; id < agent_paths.size(); ++id) {
            const auto& wps = agent_paths[id].waypoints;
            if (wps.size() < 2) continue;

            double t_cont = cost_here / step_size;  // same "time" scale
            int idx = std::clamp<int>(std::floor(t_cont),
                                      0, (int)wps.size() - 2);
            double frac = t_cont - idx;
            Eigen::Vector2d p_other =
                (1 - frac) * wps[idx] + frac * wps[idx + 1];

            if ((p - p_other).norm() <= (r1 + r1))
                return true; // collision
        }
    }
    return false;
}

bool agentsCollide(const Eigen::Vector2d& p1, double r1,
                   const Eigen::Vector2d& p2, double r2)
{
    double dist = (p1 - p2).norm();
    return dist <= (r1 + r2);
}






// Implement your RRT algorithm here
// amp::Path2D MyRRT::plan(const amp::Problem2D problem) 
amp::Path2D MyRRT::plan(const amp::Problem2D& problem) {
    // Random number generation 
    std::random_device rd; 
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<double> distX(problem.x_min, problem.x_max); 
    std::uniform_real_distribution<double> distY(problem.y_min, problem.y_max);
    std::uniform_real_distribution<double> dis(0.0, 1.0);


    // Create the base node 
    std::map<int, double> cost_to_come;
    int node_idx = 0; int goal_idx = 0;
    cost_to_come[0] = 0.0; // start node
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

        // Take a step towards that point
        Eigen::Vector2d dir = q - nodes[i_close];
        if (dir.norm() > r) dir.normalize();
        Eigen::Vector2d q_new = nodes[i_close] + dir * std::min(static_cast<double>(r), dir.norm());
        

        // Check for obstacle collisions 
        bool collision = false;
        for (const auto& obs : problem.obstacles) {
            if (diskPathIntersectsPolygon(nodes[i_close], q_new, obs.verticesCCW(), robot_radius)) {
                collision = true;
                break;
            }
        }


        // Determine the shortest path back to the start node to determine the time 
        if (!collision && (nodes.size() > 1) && (agent_paths.size() >= 1))
        {
            // Estimate time step from cumulative path length
            double edge_cost = (q_new - nodes[i_close]).norm();
            double new_cost = cost_to_come[i_close] + edge_cost;
            int t_idx = static_cast<int>(std::round(new_cost / r));
            std::vector<int> stencil = {-1, 0, 1};

            // Check against all previously planned agents
            for (size_t id = 0; id < agent_paths.size(); id++)
            {
                const auto& wp = agent_paths[id].waypoints;
                if (wp.empty()) continue;

                // Clamp index to valid range
                for (size_t i = 0; i < stencil.size(); i++)
                {
                    int idx = std::min(t_idx, static_cast<int>(wp.size()) - 1) + stencil[i];
                    if (agentsCollide(q_new, robot_radius+0.25, wp[idx], robot_radius))
                    {
                        collision = true;
                        break;
                    }
                }

                // If we prev broke, do so again
                if (collision)
                {
                    break;
                }
            }
        }

        // Add to path if valid 
        if (!collision)
        {

            // Create the new node and add it to the graph
            node_idx = node_idx + 1;
            nodes[node_idx] = q_new;
            double edge_cost = (q_new - nodes[i_close]).norm();
            cost_to_come[node_idx] = cost_to_come[i_close] + edge_cost;
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
        if (safety >= 15000) //10000
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








struct EuclideanHeuristicND : public amp::SearchHeuristic {
    // Reference to node positions in N-D
    const std::map<amp::Node, Eigen::VectorXd>& nodes;

    // Goal node coordinates in N-D
    Eigen::VectorXd goal;

    // Constructor
    EuclideanHeuristicND(const std::map<amp::Node, Eigen::VectorXd>& nodes_,
                         const Eigen::VectorXd& goal_)
        : nodes(nodes_), goal(goal_) {}

    // Override the virtual operator()
    virtual double operator()(amp::Node node) const override {
        auto it = nodes.find(node);
        if (it == nodes.end()) {
            throw std::runtime_error("Node not found in heuristic!");
        }
        const Eigen::VectorXd& pos = it->second;

        // Euclidean distance in N-D to goal
        return (pos - goal).norm();
    }
};



bool agentsCollideAtConfiguration(const Eigen::VectorXd& q, int N, double agent_radius) {
    for (int i = 0; i < N; ++i) {
        Eigen::Vector2d pos_i = q.segment<2>(i*2);
        for (int j = i+1; j < N; ++j) {
            Eigen::Vector2d pos_j = q.segment<2>(j*2);
            if ((pos_i - pos_j).norm() <= 2*agent_radius) {
                return true; // collision between agent i and j
            }
        }
    }
    return false;
}

bool agentsCollideAlongSegment(const Eigen::VectorXd& q_old,
                               const Eigen::VectorXd& q_new,
                               int N,
                               double agent_radius,
                               int n_steps = 4) // 4
{
    for (int step = 0; step <= n_steps; ++step) {
        double alpha = static_cast<double>(step) / n_steps;
        Eigen::VectorXd q_interp = q_old + alpha * (q_new - q_old);

        if (agentsCollideAtConfiguration(q_interp, N, agent_radius))
            return true;
    }
    return false;
}



amp::MultiAgentPath2D MyCentralRRT::plan(const amp::MultiAgentProblem2D& problem) {
    // Get number of agents
    int N = problem.numAgents();

    // Random number generation
    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<std::uniform_real_distribution<double>> dists;
    for (size_t n = 0; n < N; ++n) {
        dists.emplace_back(problem.x_min, problem.x_max);
        dists.emplace_back(problem.y_min, problem.y_max);
    }
    std::uniform_real_distribution<double> dis(0.0, 1.0);

    // Get init and goal 
    Eigen::VectorXd joint_q_init(2 * N);
    Eigen::VectorXd joint_q_goal(2 * N);

    for (size_t i = 0; i < N; ++i) {
        // Copy each agent's 2D initial position into the N-D vector
        joint_q_init.segment<2>(i * 2) = problem.agent_properties[i].q_init;
        joint_q_goal.segment<2>(i * 2) = problem.agent_properties[i].q_goal;
    }

    // Initialize
    nodes.clear();
    graphPtr->clear();
    int node_idx = 0, goal_idx = 0;
    nodes[node_idx] = joint_q_init;

    int safety = 0;
    Eigen::VectorXd q(2*N);

    while (true) {
        // Sample with goal bias
        if (dis(gen) < p_goal) {
            q = joint_q_goal;
        } else {
            for (size_t d = 0; d < 2*N; ++d)
                q[d] = dists[d](gen);
        }

        // Find closest node
        int i_close = 0;
        double dist_close = std::numeric_limits<double>::max();
        for (const auto& [idx, pos] : nodes) {
            double d = (pos - q).norm();
            if (d < dist_close) {
                dist_close = d;
                i_close = idx;
            }
        }

        // Step towards sample
        Eigen::VectorXd dir = q - nodes[i_close];
        if (dir.norm() > r) dir.normalize();
        Eigen::VectorXd q_new = nodes[i_close] + dir * std::min(static_cast<double>(r), dir.norm());

        // Obstacle collision (generalized for N-D)
        bool collision = false;
        for (const auto& obs : problem.obstacles) {
            for (size_t i = 0; i < N; i++)
            {

                Eigen::Vector2d q_old_agent = nodes[i_close].segment<2>(i * 2);
                Eigen::Vector2d q_new_agent = q_new.segment<2>(i * 2);


                if (diskPathIntersectsPolygon(q_old_agent, q_new_agent, obs.verticesCW(), robot_radius+0.25)) {
                    collision = true;
                    break;
                }
            }

            if (collision)
            {
                break;
            } 
        }

        // Check agent collisions (same index = same time step)
        if (!collision) {
            for (int i = 0; i < N; ++i) {
                Eigen::Vector2d pi = q_new.segment<2>(i * 2); // agent i position

                for (int j = i + 1; j < N; ++j) {
                    Eigen::Vector2d pj = q_new.segment<2>(j * 2); // agent j position
                    if (!collision && agentsCollideAlongSegment(nodes[i_close], q_new, N, robot_radius + 0.25, 10)) {
                        collision = true;
                    }
                }
            }
        }

   

        // Add node if valid
        if (!collision) {
            node_idx++;
            nodes[node_idx] = q_new;
            graphPtr->connect(i_close, node_idx, r);
            graphPtr->connect(node_idx, i_close, r);

            // Check if goal reached
            if ((q_new - joint_q_goal).norm() <= epsilon * std::sqrt(2 * N)) {
                node_idx++;
                nodes[node_idx] = joint_q_goal;
                graphPtr->connect(node_idx - 1, node_idx, r);
                graphPtr->connect(node_idx, node_idx - 1, r);
                goal_idx = node_idx;
                break;
            }
        }

        // Safety break
        if (++safety >= 50000) {
            std::cout << "Safety" << std::endl;
            break;
        }
    }

    // Extract shortest path
    amp::ShortestPathProblem sp(graphPtr, 0, goal_idx);
    EuclideanHeuristicND heuristic(nodes, nodes[goal_idx]);
    MyAStarAlgo algo;
    auto result = algo.search(sp, heuristic);
    //result.print();


    amp::MultiAgentPath2D path;
    if (!result.success)
    {
        for (const amp::CircularAgentProperties& agent : problem.agent_properties) {
            amp::Path2D agent_path;
            agent_path.waypoints = {agent.q_init, agent.q_goal};
            path.agent_paths.push_back(agent_path);
        }

        return path;
    } 


    for (size_t i = 0; i < N; ++i) {
        amp::Path2D agent_path;

        // Extract agent i's 2D slice at each joint configuration
        for (const amp::Node& node : result.node_path) {
            const Eigen::VectorXd& q_joint = nodes.at(node);
            Eigen::Vector2d q_agent = q_joint.segment<2>(i * 2);
            agent_path.waypoints.push_back(q_agent);
        }

        path.agent_paths.push_back(agent_path);
    }


    return path;
}