/*
Author: Bennet Outland
Affiliation: CU Boulder
Control: None, University Product

Resources:
- https://www.w3schools.com/cpp
- https://libeigen.gitlab.io/eigen/docs-nightly/group__QuickRefPage.html
- https://en.wikipedia.org/wiki/Line_(geometry)#Linear_equation
- ChatGPT for the point in polygon implementation
- Brin "Tea Time Numerical Analysis" 
*/

// Includes
#include "Bug1.h"
#include <Eigen/Dense>
#include <Eigen/Core>


// Raycasting algorithm to determine if in the obstacle
obsMemb Bug1::in_obstacle(amp::Problem2D problem, Eigen::Vector2d x) {
    // Determine if there is an intersection
    //obsMemb obstacleMemebership = obsMemb(true, -1);
    bool valid = true;
    int idx = -1;

    // Happy direction vector!
    Eigen::Vector2d d(1.0, 0.0);

    // Iterate through each obstacle
    for (size_t i = 0; i < problem.obstacles.size(); i++)
    {
        // Initialize a count variable
        int count = 0;

        // Define a verticies variable for conviencece 
        std::vector<Eigen::Vector2d> vertices = problem.obstacles[i].verticesCCW();

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

        if (count % 2 == 1)
        {
            valid = false;
            idx = i;
            break;
        }
    }

    // Determine if an intersection was found
    return obsMemb(valid, idx);
    
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


/* 
Determines which obstacles are connected in an adjacency list
*/
std::vector<std::vector<int>> Bug1::connected_obstacles(amp::Problem2D problem, float r_obs) {
    // Adjacency list
    std::vector<std::vector<int>> adj;
    adj.resize(problem.obstacles.size());

    // Loop through all of the obstacles 
    for (size_t i = 0; i < problem.obstacles.size(); i++) {
        for (size_t j = 0; j < problem.obstacles.size(); j++) {

            // Exclude self assesment
            if (i != j)
            {
                // Intersection
                bool intersect = false;

                // Determine if the vertices or centroid of obstacle j is in obstacle i
                std::vector<Eigen::Vector2d> v_i = problem.obstacles[i].verticesCCW();
                std::vector<Eigen::Vector2d> v_j = problem.obstacles[j].verticesCCW();

                // Test verticies and add up for centroid
                Eigen::Vector2d c_j(0.0, 0.0);
                bool flag = false;

                for (size_t k = 0; k < v_j.size(); k++)
                {
                    if (pip(v_i, v_j[k]))
                    {
                        // Add to polygon to adjacency matrix
                        adj[i].push_back(j);
                        flag = true;
                    }

                    // Add up for centroid 
                    c_j = c_j + v_j[k];
                }

                // Test the centroid if not already added
                if (!flag && pip(v_i, c_j / v_j.size()))
                {
                    adj[i].push_back(j);
                }
                
            }  
        }
    }

    // Reduce down the adjacency list
    return adj;
}




/* 
Gets within a desired distance of the boundary. Make sure offset distance sign is correct good
*/ 
Eigen::Vector2d Bug1::get_to_boundary(amp::Problem2D problem, int idx, Eigen::Vector2d wpt, Eigen::Vector2d dir, float dx, float e_o) {
    // Variables
    Eigen::Vector2d x = wpt;

    // Determine if the point is inside the polygon or not
    bool inside = pip(problem.obstacles[idx].verticesCCW(), x);   
    //std::cout << inside; 

    if (inside)
    {

        Eigen::Vector2d t_dir = dir; // dir


        // Domain where we can find the boundary
        std::vector<float> dom = {-2*dx, 0}; 

        // Determine the boundary location by bisection 
        int N = ceil((log(dom[1] - dom[0]) - log(0.005)) / log(2)) + 1;
        float m = (dom[1] + dom[0])/2;
        for (size_t i = 0; i < N; i++)
        {
            // Set the bisection point
            m = (dom[1] + dom[0])/2;
            bool M = pip(problem.obstacles[idx].verticesCCW(), x + m*t_dir);

            if (M)
            {
                dom = {dom[0], m};
            } else {
                dom = {m, dom[1]};
            }
            
        }
        m = (dom[1] + dom[0])/2;

        // Return the optimized distance; 
        return x + m*t_dir - e_o*dir; 

    } 
    
}


Eigen::Vector2d find_a_boundary(amp::Problem2D problem, int idx, Eigen::Vector2d q, Eigen::Vector2d dir, float dx, float e_o) {

    // Defining opt variables 
    float pi = 3.14;
    float theta = pi/2; 
    float theta_prev = 0.0;
    float dtheta = -0.01;
    int N = ceil(abs((3*pi/2) / dtheta));
    Eigen::Matrix<double, 2, 2> R, R_safe; R << cos(theta), -sin(theta), sin(theta), cos(theta);
    Eigen::Vector2d t_dir, t_dir_safe, q_b0, q_b1, n_obs;

    // Rotate until within e_o of boundary 
    for (size_t i = 0; i < N; i++)
    {
        // Get the vector
        R << cos(theta), -sin(theta), sin(theta), cos(theta);
        t_dir = R * dir;

        // See if it falls into the polygon 
        if (pip(problem.obstacles[idx].verticesCCW(), q + dx*t_dir + e_o*t_dir)) {
            // Take the previous safe angle
            R_safe << cos(theta_prev), -sin(theta_prev), sin(theta_prev), cos(theta_prev);
            t_dir_safe = R_safe * dir;

            // Get the direction normal to the obstacle 
            n_obs = {-t_dir_safe[1], t_dir_safe[0]};

            return q + dx*t_dir_safe + n_obs * e_o;
        } else {
            // Save as a safe theta value
            theta_prev = theta;
            theta = theta + dtheta;
        }
    }

    return q;
    
}


/*
Used to determine the path for moving around the boundary of an object for a left turning robot. 
*/
amp::Path2D Bug1::move_along_boundary(amp::Path2D path, amp::Problem2D problem, Eigen::Vector2d q, Eigen::Vector2d dir, int idx, std::vector<std::vector<int>> adj, float dx, float e_o) {
    // Register the hit point and prev point
    Eigen::Vector2d q_H = q;
    Eigen::Vector2d q_p = q;
    Eigen::Vector2d tv;

    // First turn is left
    Eigen::Vector2d dv(dir[1], -dir[0]); // get direction vector
    q = q + dx*dv; // get a canidate waypoint
    q = get_to_boundary(problem, idx, q, dv, dx, e_o); // do one correction 

    // Update the path 
    path.waypoints.push_back(q);

    // Save closest to goal and indicies
    int q_H_idx = path.waypoints.size() - 1; // - 1 off by one?
    float q_L_dist = (problem.q_goal - q_H).norm(); // Initialize as hit point
    int q_L_idx = path.waypoints.size() - 1; // - 1 off by one?

    // std::cout << q_L_idx << " \n";

    float dist;
    std::vector<Eigen::Vector2d> wpts;

    // Go along the boundary 
    for (size_t i = 1; i < 2000; i++)
    {
        // Get previous direction vector 
        Eigen::Vector2d tv = (q - q_p) / (q - q_p).norm();
        q_p = q;

        // Find the next point on the boundary
        q = find_a_boundary(problem, idx, q, tv, dx, e_o);

        // Get distance to goal and index and see if closer 
        dist = (problem.q_goal - q).norm();
        if (dist < q_L_dist)
        {
            q_L_dist = dist; 
            q_L_idx = path.waypoints.size() + 1;
        }
        

        // Determine if it is close-ish to the hit point 
        if ((q - q_H).norm() < 2*dx)
        {
            // Return to hit point
            path.waypoints.push_back(q_H);

            // Refollow the path to q_L in reverse order 
            for (size_t j = path.waypoints.size()-1; j >= q_L_idx; j--)
            {
                path.waypoints.push_back(path.waypoints[j]);
            }

            break;
        } else {
            path.waypoints.push_back(q);
        }
        
    }


    return path;
}




// Move to the goal 
m2g Bug1::move_to_goal(amp::Problem2D problem, amp::Path2D path, Eigen::Vector2d q_start, float dx, float r_g) {
    // Det start and goal values 
    Eigen::Vector2d q_init = q_start;
    Eigen::Vector2d q_goal = problem.q_goal;

    // Include the initial point
    path.waypoints.push_back(q_init);

    // Determine the line parameters 
    Eigen::Vector2d dir = (q_goal - q_init) / (q_goal - q_init).norm();

    // Current location and waypoints
    Eigen::Vector2d wpt = q_init;

    // Safety counter 
    int safety = 0;

    // Obstacle index 
    int obs = -1;

    while ((path.waypoints.back() - q_goal).norm() > r_g)
    {
        // Create the next waypoint 
        wpt = wpt + dx*dir;

        // Check to see if a waypoint is in an obstacale 
        obsMemb wpt_info = in_obstacle(problem, wpt);
        bool valid_wpt = wpt_info.valid;
        obs = wpt_info.obs;

        if (valid_wpt)
        {
            // Push the waypoint
            path.waypoints.push_back(wpt);
        } else {
            // Get close to the boundary without entering
            Eigen::Vector2d b_wpt = get_to_boundary(problem, wpt_info.obs, wpt, dir, dx, 0.01);
            path.waypoints.push_back(b_wpt);
            break;
        }

        // Loop end condition if something bad happens
        safety = safety + 1;
        if (safety >= 10000) {
            std::cout << "safety!";
            break;
        }

    }

    // If really close to the goal, go to the goal
    if ((path.waypoints.back() - q_goal).norm() <= r_g) {
        path.waypoints.push_back(q_goal);
    }
    
    // Build the object 
    m2g m2g_obj(path, dir, obs);

    return m2g_obj;
}




// Implement your methods in the `.cpp` file, for example:
amp::Path2D Bug1::plan(const amp::Problem2D& problem) {

    // Your algorithm solves the problem and generates a path.
    amp::Path2D path;

    // Parameters 
    float dx = 0.05;
    float dx_b = 0.005; // 0.005
    float r_g = 0.2;
    float e_o = 0.01;

    // Leave point 
    Eigen::Vector2d q_L = problem.q_init;
    Eigen::Vector2d q_H;


    // Determine connected obstacles 
    std::vector<std::vector<int>> adj_list = connected_obstacles(problem, 0.01);

    int safety = 0;
    while (safety < 10) {
    
        // Move towards the goal until stopped by an obstacle
        m2g m2g_obj = move_to_goal(problem, path, q_L, dx, r_g);
        path = m2g_obj.path;

        // set a hit point 
        if (path.waypoints.back() == problem.q_goal)
        {
            std::cout << "Done!";
            //break;
        } else {
            // set a hit point 
            Eigen::Vector2d q_H = path.waypoints.back();

            path = move_along_boundary(path, problem, q_H, m2g_obj.dir, m2g_obj.idx, adj_list, dx_b, e_o);
        }

        q_L = path.waypoints.back();


        safety = safety + 1;
    }

    return path;
}
