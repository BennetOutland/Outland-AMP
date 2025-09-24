/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- Lecture slides
- Roboscience book for the IK
- ChatGPT for the random library usage
*/


#include "ManipulatorSkeleton.h"
#include <random>

MyManipulator2D::MyManipulator2D(std::vector<double> links)
    : LinkManipulator2D(links) // Default to a 2-link with all links of 1.0 length
{}

// Override this method for implementing forward kinematics
Eigen::Vector2d MyManipulator2D::getJointLocation(const amp::ManipulatorState& state, uint32_t joint_index) const {
    // Implement forward kinematics to calculate the joint position given the manipulator state (angles)

    // Get link lengths
    std::vector<double> a = LinkManipulator2D::getLinkLengths();

    // Calculate the joint locations
    if (nLinks() == 2) {
        Eigen::Vector2d p0 = Eigen::Vector2d(0, 
                                            0);
        Eigen::Vector2d p1 = Eigen::Vector2d(a[0] * cos(state[0]), 
                                            a[0] * sin(state[0]));
        Eigen::Vector2d p2 = Eigen::Vector2d(a[0] * cos(state[0]) + a[1] * cos(state[0] + state[1]), 
                                            a[0] * sin(state[0]) + a[1] * sin(state[0] + state[1]));

        // Package up
        std::vector<Eigen::Vector2d> joint_positions = {p0, p1, p2};

        // Ship out
        return joint_positions[joint_index];
    } else {
        Eigen::Vector2d p0 = Eigen::Vector2d(0, 
                                            0);
        Eigen::Vector2d p1 = Eigen::Vector2d(a[0] * cos(state[0]), 
                                            a[0] * sin(state[0]));
        Eigen::Vector2d p2 = Eigen::Vector2d(a[0] * cos(state[0]) + a[1] * cos(state[0] + state[1]), 
                                            a[0] * sin(state[0]) + a[1] * sin(state[0] + state[1]));
        Eigen::Vector2d p3 = Eigen::Vector2d(a[0] * cos(state[0]) + a[1] * cos(state[0] + state[1]) + a[2] * cos(state[0] + state[1] + state[2]), 
                                            a[0] * sin(state[0]) + a[1] * sin(state[0] + state[1]) + a[2] * sin(state[0] + state[1] + state[2]));

        // Package up
        std::vector<Eigen::Vector2d> joint_positions = {p0, p1, p2, p3};

        // Ship out
        return joint_positions[joint_index];
    }
}

// Override this method for implementing inverse kinematics
amp::ManipulatorState MyManipulator2D::getConfigurationFromIK(const Eigen::Vector2d& end_effector_location) const {
    // Implement inverse kinematics here

    // Get link lengths
    std::vector<double> a = LinkManipulator2D::getLinkLengths();

    // For better readability
    double x = end_effector_location[0];
    double y = end_effector_location[1];
    
    // If you have different implementations for 2/3/n link manipulators, you can separate them here
    if (nLinks() == 2) {
        // Get the trig funcs of the angles
        double cos_theta_2 = ((1.0 / (2*a[0]*a[1])) * (pow(x, 2) + pow(y, 2) - pow(a[0], 2) - pow(a[1], 2)));
        double sin_theta_2 = sqrt(1 - pow(cos_theta_2, 2)); // choosing positive bc i happy :(

        // D constant 
        double D = (pow(x, 2) + pow(y, 2) - pow(a[0], 2) - pow(a[1], 2)) / (2.0 * a[0] * a[1]);

        // Get specific angles 
        double theta_1 = atan2(y, x) - atan2(a[1]*sin_theta_2, a[0] + a[1]*cos_theta_2);
        double theta_2 = atan2(sqrt(1 - pow(D, 2)), D);

        // Apply
        Eigen::Vector2d joint_angles = Eigen::Vector2d(theta_1, theta_2);

        return joint_angles;
    } else if (nLinks() == 3) {

        std::random_device rd;
        std::mt19937 gen(rd());

        for (size_t i = 0; i < 100; i++)
        {
            // Pick a theta_1 in the direction of the point
            double pi = 3.14159;
            std::uniform_real_distribution<float> dist(-pi/2, pi/2);
            double theta_1 = atan2(y, x) + dist(gen);

            // offset the end effector location 
            double x_o = x * cos(theta_1) + y * sin(theta_1) - a[0];
            double y_o = -1*x * sin(theta_1) + y * cos(theta_1);

            // Determine theta_3
            double cos_theta_3 = (pow(x_o, 2) + pow(y_o, 2) - pow(a[1], 2) - pow(a[2], 2)) / (2.0 * a[1] * a[2]);
            double theta_3 = acos(cos_theta_3);

            // Determine theta_2
            double c1 = a[1] + a[2] * cos(theta_3);
            double c2 = a[2] * sin(theta_3);
            double theta_2 = atan2(y_o, x_o) - atan2(c2, c1);

            // Collect 
            Eigen::Vector3d joint_angles = Eigen::Vector3d(theta_1, theta_2, theta_3);

            // Check if the solution is good
            Eigen::Vector2d pos = MyManipulator2D::getJointLocation(joint_angles, 3);
            if ((end_effector_location - pos).norm() < 1e-6){
                return joint_angles;
            }

        }
        
        
        
    } 

    return  Eigen::Vector3d(0, 0, 0);
}