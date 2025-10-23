/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None
License: NA

Resources:
- ChatGPT for debigging
- ChatGPT for polygon centroid
- Lecture slides
*/


#include "MyGDAlgorithm.h"
#include <Eigen/Dense>


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

// Get the centroid of the polygon
Eigen::Vector2d polygonCentroid(const std::vector<Eigen::Vector2d>& vertices) {

    double signedArea = 0.0;
    double Cx = 0.0;
    double Cy = 0.0;

    for (size_t i = 0; i < vertices.size(); ++i) {
        const Eigen::Vector2d& p0 = vertices[i];
        const Eigen::Vector2d& p1 = vertices[(i + 1) % vertices.size()]; // wrap-around
        double a = p0[0] * p1[1] - p1[0] * p0[1];
        signedArea += a;
        Cx += (p0[0] + p1[0]) * a;
        Cy += (p0[1] + p1[1]) * a;
    }

    signedArea *= 0.5;

    Cx /= (6.0 * signedArea);
    Cy /= (6.0 * signedArea);

    return Eigen::Vector2d(Cx, Cy);
}


//=============================================================================================//
//                                    POTENTIAL FUNCTION
//=============================================================================================//

// Class constructor
MyPotentialFunction::MyPotentialFunction(double d_star, double zeta, double Q_star, double eta) : d_star(d_star),
    zeta(zeta),
    Q_star(Q_star),
    eta(eta)  {
}

// Add the boundary conditions 
void MyPotentialFunction::add_boundary_values(amp::Problem2D problem) {
    MyPotentialFunction::q_init = problem.q_init;
    MyPotentialFunction::q_goal = problem.q_goal;
    MyPotentialFunction::problem = problem;
}

// Returns the potential function value (height) for a given 2D point. 
double MyPotentialFunction::operator()(const Eigen::Vector2d& q) const {
    // Initialize the potential
    double potential;

    // Attractive
    double dist = (q - q_goal).norm();
    if (dist <= d_star)
    {
        potential = 0.5 * zeta * dist * dist;
    } else {
        potential = d_star * zeta * dist - (0.5 * zeta * d_star * d_star);
    }

    // Repulsive
    double divisions = 1000.0;
    for (size_t i = 0; i < problem.obstacles.size(); i++)
    {
        // Obstacle parameters
        auto obs = problem.obstacles[i];
        Eigen::Vector2d q_obs = polygonCentroid(obs.verticesCCW());

        // std::cout << q_obs << "\n";

        // Get the straightline distance to the obstacle 
        double d_i = 0.0;
        for (int l = 0; l <= divisions; ++l) {
            double t = static_cast<double>(l) / divisions;
            Eigen::Vector2d test_pt = (1 - t) * q + t * q_obs; 
            if (pip(obs.verticesCCW(), test_pt)) {
                t = static_cast<double>(l-1) / divisions;
                test_pt = (1 - t) * q + t * q_obs; 
                d_i = (q - test_pt).norm();
                break;
            }
        }


        // Add potential if too close
        if (d_i <= Q_star)
        {
            potential = potential + 0.5 * eta * pow(((1.0 / d_i) - (1.0 / Q_star)), 2.0);
        }
    }


    return potential;
}

Eigen::Vector2d MyPotentialFunction::getGradient(const Eigen::Vector2d& q) const {
    // Initialize the potential
    Eigen::Vector2d grad;

    // Attractive
    double dist = (q - q_goal).norm();
    if (dist <= d_star)
    {
        grad = 0.5 * zeta * (q - q_goal);
    } else {
        grad = (d_star * zeta * (q - q_goal)) / dist;
    }

    // Repulsive
    double divisions = 1000.0;
    for (size_t i = 0; i < problem.obstacles.size(); i++)
    {
        // Obstacle parameters
        auto obs = problem.obstacles[i];
        Eigen::Vector2d q_obs = polygonCentroid(obs.verticesCCW());

        // Get the straightline distance to the obstacle 
        double d_i = 0.0;
        Eigen::Vector2d c;
        for (int l = 0; l <= divisions; ++l) {
            double t = static_cast<double>(l) / divisions;
            Eigen::Vector2d test_pt = (1 - t) * q + t * q_obs; 
            if (pip(obs.verticesCCW(), test_pt)) {
                d_i = (q - test_pt).norm();
                c = test_pt;
                break;
            }
        }

        // Add potential if too close
        if (d_i <= Q_star)
        {
            // Distance gradient approximation
            Eigen::Vector2d dd_i = (q - c) / d_i;

            // Update to the potential gradient
            grad = grad + eta * ((1.0 / Q_star) - (1.0 / d_i)) * (dd_i / (d_i * d_i));
        }
    }

    return grad;
}


//=============================================================================================//
//                                    GRADIENT DESCENT
//=============================================================================================//

// Assign params nd create the potential function
MyGDAlgorithm::MyGDAlgorithm(double d_star, double zeta, double Q_star, double eta) :
    d_star(d_star),
    zeta(zeta),
    Q_star(Q_star),
    eta(eta),
    pot_func(d_star, zeta, Q_star, eta) {
}


// Implement your plan method here, similar to HW2:
amp::Path2D MyGDAlgorithm::plan(const amp::Problem2D& problem) {
    // Define the boundary values
    pot_func.add_boundary_values(problem);

    // Define and start the path
    amp::Path2D path;
    path.waypoints.push_back(problem.q_init);

    // Follow the gradient until close enough to the goal
    for (size_t k = 0; k < 5000; k++)
    {
        path.waypoints.push_back(path.waypoints.back() - 0.01*pot_func.getGradient(path.waypoints.back()));

        if ((path.waypoints.back() - problem.q_goal).norm() <= 0.25) {
            path.waypoints.push_back(problem.q_goal);
            break;
        }

    }
    
    return path;
}
