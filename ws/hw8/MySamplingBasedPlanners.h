#pragma once

// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework headers
#include "hw/HW7.h"


class MyRRT : public amp::GoalBiasRRT2D {
    public:
        // Constructor 
        MyRRT(double p_goal_in, float r_in, float epsilon_in, double robot_radius_in, std::vector<amp::Path2D> agent_paths_in) : p_goal(p_goal_in), r(r_in), epsilon(epsilon_in), robot_radius(robot_radius_in), agent_paths(agent_paths_in) {}

        // Parameters 
        double robot_radius;
        double p_goal;
        std::vector<amp::Path2D> agent_paths;
        float r;
        float epsilon;
        bool plan_success = false;

        // Graph 
        std::shared_ptr<amp::Graph<double>> graphPtr = std::make_shared<amp::Graph<double>>();
        std::map<amp::Node, Eigen::Vector2d> nodes;

        // Planning function
        virtual amp::Path2D plan(const amp::Problem2D& problem) override; 
};



class MyCentralRRT {
    public:
        // Constructor 
        MyCentralRRT(double p_goal_in, float r_in, float epsilon_in, double robot_radius_in) : p_goal(p_goal_in), r(r_in), epsilon(epsilon_in), robot_radius(robot_radius_in) {}

        // Parameters 
        double robot_radius;
        double p_goal;
        float r;
        float epsilon;
        bool plan_success = false;

        // Graph 
        std::shared_ptr<amp::Graph<double>> graphPtr = std::make_shared<amp::Graph<double>>();
        std::map<amp::Node, Eigen::VectorXd> nodes;

        // Planning function
        amp::MultiAgentPath2D plan(const amp::MultiAgentProblem2D& problem); 
};