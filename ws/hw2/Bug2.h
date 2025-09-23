#pragma once

#include "AMPCore.h"
#include "hw/HW2.h"

struct obsMemb {
    bool valid;
    int obs;
};


struct m2g {
    amp::Path2D path;
    Eigen::Vector2d dir;
    int idx;
};

/// @brief Declare your bug algorithm class here. Note this class derives the bug algorithm class declared in HW2.h
class Bug2 : public amp::BugAlgorithm {
    public:
        // Override and implement the bug algorithm in the plan method. The methods are declared here in the `.h` file
        virtual amp::Path2D plan(const amp::Problem2D& problem) override;

        // Class Methods 
        obsMemb in_obstacle(amp::Problem2D problem, Eigen::Vector2d x);
        m2g move_to_goal(amp::Problem2D problem, amp::Path2D path, Eigen::Vector2d q_start, float dx, float r_g);
        Eigen::Vector2d get_to_boundary(amp::Problem2D problem, int idx, Eigen::Vector2d wpt, Eigen::Vector2d dir, float dx, float e_o);
        std::vector<std::vector<int>> connected_obstacles(amp::Problem2D problem, float r_obs);
        amp::Path2D move_along_boundary(amp::Path2D path, amp::Problem2D problem, Eigen::Vector2d q, Eigen::Vector2d dir, int idx, std::vector<std::vector<int>> adj, float dx, float e_o);
        bool pip(std::vector<Eigen::Vector2d> vertices, Eigen::Vector2d x);
        Eigen::Vector2d find_a_boundary(amp::Problem2D problem, int idx, Eigen::Vector2d q, Eigen::Vector2d dir, float dx, float e_o);
    
    private:
        // Add any member variables here...
};