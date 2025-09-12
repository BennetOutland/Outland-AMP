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
Eigen::Vector2d get_to_boundary_first(amp::Problem2D problem, int idx, Eigen::Vector2d wpt, Eigen::Vector2d dir, float dx, float e_o) {
    Eigen::Vector2d x = wpt;

    // Domain where we can find the boundary
    std::vector<float> dom = {-dx, 0}; 

    // Determine the boundary location by bisection 
    int N = ceil((log(dom[1] - dom[0]) - log(0.005)) / log(2)) + 1;
    float m = (dom[1] + dom[0])/2;
    for (size_t i = 0; i < N; i++)
    {
        // Set the bisection point
        m = (dom[1] + dom[0])/2;
        bool M = pip(problem.obstacles[idx].verticesCCW(), x + m*dir);

        if (M)
        {
            dom = {dom[0], m};
        } else {
            dom = {m, dom[1]};
        }
        
    }
    m = (dom[1] + dom[0])/2;

    Eigen::Vector2d t_dir = {-dir[1], dir[0]};

    // Return the optimized distance; 
    return x + m*dir; //+ t_dir*0.001; 
}



/* 
Gets within a desired distance of the boundary. Make sure offset distance sign is correct good
*/ 
Eigen::Vector2d Bug1::get_to_boundary(amp::Problem2D problem, int idx, Eigen::Vector2d wpt, Eigen::Vector2d dir, float dx, float e_o) {
    // Variables
    Eigen::Vector2d x = wpt;

    // Determine if the point is inside the polygon or not
    bool inside = pip(problem.obstacles[idx].verticesCCW(), x);   
    std::cout << inside; 

    if (inside)
    {

        // Rotate direction vector to the left of local frame 
        // double theta = 1.5707;
        // Eigen::Matrix<double, 2, 2> R; R << cos(theta), -sin(theta), sin(theta), cos(theta);
        Eigen::Vector2d t_dir = dir; //R * dir;


        // Domain where we can find the boundary
        std::vector<float> dom = {-dx, 0}; 

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

        // Eigen::Vector2d t_dir = {-dir[1], dir[0]};

        // Return the optimized distance; 
        return x + m*t_dir; //+ t_dir*0.001; 

    } else {
        
        // Rotate direction vector to the right of local frame 
        double theta = -1.5707;
        Eigen::Matrix<double, 2, 2> R; R << cos(theta), -sin(theta), sin(theta), cos(theta);
        Eigen::Vector2d t_dir = R * dir;

        // Domain where we can find the boundary
        std::vector<float> dom = {0, dx}; 

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



        // Return the optimized distance
        return x + m*t_dir;// - e_o*t_dir;
    }

    
}



/*
Used to determine the path for moving around the boundary of an object for a left turning robot. 
*/
amp::Path2D Bug1::move_along_boundary(amp::Path2D path, amp::Problem2D problem, Eigen::Vector2d q, Eigen::Vector2d dir, int idx, std::vector<std::vector<int>> adj, float dx, float e_o) {
    // Register the hit point and prev point
    Eigen::Vector2d q_H = q;
    Eigen::Vector2d q_p = q;

    // First turn is left
    Eigen::Vector2d dv(-dir[1], dir[0]); // get direction vector
    q = q + dx*dv; // get a canidate waypoint
    q = get_to_boundary_first(problem, idx, q, dv, dx, e_o); // do one correction 

    // Update the path 
    path.waypoints.push_back(q);

    // Go along the boundary 
    Eigen::Vector2d q_c = path.waypoints.back();
    Eigen::Vector2d t_obs = (q_c - q_p) / (q_c - q_p).norm();
    bool inside_prev = true;

    for (size_t i = 1; i < 3; i++)
    {
        // Get the direction tangent to the surface 
        q_c = path.waypoints.back();
        t_obs = (q_c - q_p) / (q_c - q_p).norm();
        q_p = q_c;

        // Get a canidate 
        q = q + dx*t_obs; 

        // Get it to a boundary 
        q = get_to_boundary(problem, idx, q, t_obs, dx, e_o);



        // // Get the normal vector 
        // double theta = 1.5707;
        // Eigen::Matrix<double, 2, 2> R; R << cos(theta), -sin(theta), sin(theta), cos(theta);
        // Eigen::Vector2d n_obs = R * t_obs;


        // // Adding some stability
        // if (inside_prev)
        // {
        //     q = q + 0.05*dx*n_obs;
        //     inside_prev = false;
        // }

        


        // Update the path 
        path.waypoints.push_back(q);
    }
    

    


    return path;
}




// Move to the goal 
m2g Bug1::move_to_goal(amp::Problem2D problem, float dx, float r_g) {
    // Det start and goal values 
    Eigen::Vector2d q_init = problem.q_init;
    Eigen::Vector2d q_goal = problem.q_goal;

    // Define an initial path
    amp::Path2D path;

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


    // Determine connected obstacles 
    std::vector<std::vector<int>> adj_list = connected_obstacles(problem, 0.01);

    int safety = 0;
    // while (safety < 10000) {
    
    // Move towards the goal until stopped by an obstacle
    m2g m2g_obj = move_to_goal(problem, dx, r_g);
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





    //     safety = safety + 1;
    // }
        

    // go around the boundary and calculate distances etc.

    // determine best point and travel to it 

    // continue moving to goal

    return path;
}
