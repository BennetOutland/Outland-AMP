/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- ChatGPT misc debugging
- Lecture slides
- ChatGPT for psuedocode
*/

#include <Eigen/Dense>
#include <algorithm>
#include <Eigen/Core>
#include <queue>
#include <vector>


#include "MyCSConstructors.h"

////////////////////// THIS IS FROM HW4 //////////////////////

/* You can just move these classes to shared folder and include them instead of copying them to hw6 project*/


std::pair<std::size_t, std::size_t> MyGridCSpace2D::getCellFromPoint(double x0, double x1) const {

    double pi = 3.14159;
    double x_min = x0_min;
    double x_max = x0_max;
    double y_min = x1_min;
    double y_max = x1_max;

    float m_cells_per_dim = float(N);  // TODO ===================================================================== Magic number

    std::size_t cell_x = static_cast<std::size_t>(
        (x0 - x_min) / ((x_max - x_min) / m_cells_per_dim)
    );
    std::size_t cell_y = static_cast<std::size_t>(
        (x1 - y_min) / ((y_max - y_min) / m_cells_per_dim)
    );

    // Clamp to valid range
    if (cell_x >= m_cells_per_dim) cell_x = m_cells_per_dim - 1;
    if (cell_y >= m_cells_per_dim) cell_y = m_cells_per_dim - 1;

    return {cell_x, cell_y};
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



std::unique_ptr<amp::GridCSpace2D> MyManipulatorCSConstructor::construct(const amp::LinkManipulator2D& manipulator, const amp::Environment2D& env) {
    // Create an object of my custom cspace type (e.g. MyGridCSpace2D) and store it in a unique pointer. 
    // Pass the constructor parameters to std::make_unique()

    // Values
    double pi = 3.14159;
    double x_min = 0;
    double x_max = 2*pi;
    double y_min = 0;
    double y_max = 2*pi;


    std::unique_ptr<MyGridCSpace2D> cspace_ptr = std::make_unique<MyGridCSpace2D>(m_cells_per_dim, m_cells_per_dim, x_min, x_max, y_min, y_max);
    // In order to use the pointer as a regular GridCSpace2D object, we can just create a reference
    MyGridCSpace2D& cspace = *cspace_ptr;

    // Determine if each cell is in collision or not, and store the values the cspace. This `()` operator comes from DenseArray base class
    double dx = (x_max - x_min) / m_cells_per_dim;
    double dy = (y_max - y_min) / m_cells_per_dim;

    // Line increment parameter 
    int divisions = 250; // 50

    for (int kx = 0; kx < m_cells_per_dim; ++kx) {
        double x_i = x_min + (kx + 0.5) * dx;
        for (int ky = 0; ky < m_cells_per_dim; ++ky) {
            double y_i = y_min + (ky + 0.5) * dy;

            Eigen::Vector2d pt(x_i, y_i);
            Eigen::Vector2d j_1 = manipulator.getJointLocation(pt, 1);
            Eigen::Vector2d j_2 = manipulator.getJointLocation(pt, 2);

            //std::cout << j_1;

            bool collision = false;

            for (const auto& obs : env.obstacles) {
                // Sample first link
                for (int l = 0; l <= divisions; ++l) {
                    double t = static_cast<double>(l) / divisions;
                    Eigen::Vector2d test_pt = (1 - t) * Eigen::Vector2d(0,0) + t * j_1;
                    if (pip(obs.verticesCCW(), test_pt)) {
                        collision = true;
                        break;
                    }
                }
                if (collision) break;

                // Sample second link
                for (int l = 0; l <= divisions; ++l) {
                    double t = static_cast<double>(l) / divisions;
                    Eigen::Vector2d test_pt = (1 - t) * j_1 + t * j_2;
                    if (pip(obs.verticesCCW(), test_pt)) {
                        collision = true;
                        break;
                    }
                }
                if (collision) break;
            }

            cspace(kx, ky) = collision;
            //cspace(ky, kx) = collision;
        }
    }

    return cspace_ptr;

}

//////////////////////////////////////////////////////////////

// Override this method for computing all of the boolean collision values for each cell in the cspace
std::unique_ptr<amp::GridCSpace2D> MyPointAgentCSConstructor::construct(const amp::Environment2D& env) {
    // Create an object of my custom cspace type (e.g. MyGridCSpace2D) and store it in a unique pointer. 
    // Pass the constructor parameters to std::make_unique()
    std::unique_ptr<MyGridCSpace2D> cspace_ptr = std::make_unique<MyGridCSpace2D>(m_cells_per_dim, m_cells_per_dim, env.x_min, env.x_max, env.y_min, env.y_max);
    // In order to use the pointer as a regular GridCSpace2D object, we can just create a reference
    MyGridCSpace2D& cspace = *cspace_ptr;

    // Values 
    double x_min = env.x_min;
    double x_max = env.x_max;
    double y_min = env.y_min;;
    double y_max = env.y_max;;
  
    // Determine if each cell is in collision or not, and store the values the cspace. This `()` operator comes from DenseArray base class
    double dx = (x_max - x_min) / m_cells_per_dim;
    double dy = (y_max - y_min) / m_cells_per_dim;

    // Line increment parameter 
    int divisions = 250; // 50

    for (int kx = 0; kx < m_cells_per_dim; ++kx) {
        double x_i = x_min + (kx + 0.5) * dx;
        for (int ky = 0; ky < m_cells_per_dim; ++ky) {
            double y_i = y_min + (ky + 0.5) * dy;

            Eigen::Vector2d pt(x_i, y_i);
    

            bool collision = false;

            for (const auto& obs : env.obstacles) {

                if (pip(obs.verticesCCW(), pt)) {
                    collision = true;
                    break;
                }

            }

            cspace(kx, ky) = collision;
        }
    }

    return cspace_ptr;

}



Eigen::Vector2d getPointFromCell(int x, int y, const amp::GridCSpace2D& grid_cspace) {
    // Hard code x_min.x_max etc. to be in the c space not the workspace

    double x_min = grid_cspace.x0Bounds().first;
    double x_max = grid_cspace.x0Bounds().second;
    double y_min = grid_cspace.x1Bounds().first;
    double y_max = grid_cspace.x1Bounds().second;

    float m_cells_per_dim = float(grid_cspace.size().first); 

    //double dx = (x_max - x_min) / m_cells_per_dim;
    //double dy = (y_max - y_min) / m_cells_per_dim;
    double dx = (x_max - x_min) / grid_cspace.size().first;
    double dy = (y_max - y_min) / grid_cspace.size().second;

    return Eigen::Vector2d(x_min + double(x) * dx + dx/2, y_min + double(y)*dy + dy/2);
}



amp::Path2D MyWaveFrontAlgorithm::planInCSpace(const Eigen::Vector2d& q_init, const Eigen::Vector2d& q_goal, const amp::GridCSpace2D& grid_cspace, bool isManipulator) {
    // Implement your WaveFront algorithm here
    amp::Path2D path;

    // Get init and goal cells
    std::pair<std::size_t, std::size_t> q_init_cell_pair = grid_cspace.getCellFromPoint(q_init[0], q_init[1]);
    std::pair<std::size_t, std::size_t> q_goal_cell_pair = grid_cspace.getCellFromPoint(q_goal[0], q_goal[1]);
    Eigen::Vector2i q_init_cell(q_init_cell_pair.first, q_init_cell_pair.second);
    Eigen::Vector2i q_goal_cell(q_goal_cell_pair.first, q_goal_cell_pair.second);

    // Make the grid
    int H = grid_cspace.size().second;
    int W = grid_cspace.size().first;
    Eigen::MatrixXi grid = Eigen::MatrixXi::Zero(H, W);

    // Initialize grid values: obstacle: 1, goal: 2, others: 0
    for (int i = 0; i < W; i++)
    {
        for (int j = 0; j < H; j++)
        {
            if (grid_cspace(i, j))
            {
                grid(j, i) = 1;
            } else {
                grid(j, i) = 0;
            }            
        }   
    }
    grid(q_goal_cell_pair.second, q_goal_cell_pair.first) = 2;

    // BFS queue
    std::queue<Eigen::Vector2i> q;
    q.push(q_goal_cell);
    
    // Connected Neighbors 
    const int dx[4] = { 1,  0, -1,  0};
    const int dy[4] = { 0,  1,  0, -1};

    // Flood-fill 
    while (!q.empty()) {
        // Get the current cell
        Eigen::Vector2i c = q.front(); q.pop();
        int cval = grid(c.y(), c.x());

        // Search over the neighbors
        for (int k = 0; k < 4; ++k) {
            int nx = c.x() + dx[k];
            int ny = c.y() + dy[k];
            
            // Wrapping if needed
            if (isManipulator)
            {
                nx = (nx + W) % W;
                ny = (ny + H) % H;
            } else {
                if (ny < 0 || ny >= H || nx < 0 || nx >= W)
                {
                    continue;
                }  
            }

            // See if this cell has been visited yet
            int& cellVal = grid(ny, nx);
            if (cellVal == 0) { 
                cellVal = cval + 1;
                q.push(Eigen::Vector2i(nx, ny));
            }
        }
    }

    //std::cout << grid;

    // There is no path
    if (grid(q_init_cell_pair.second, q_init_cell_pair.first) == 0) {
        std::cout << "No path exists";
        return path;
    }


    // Create the c-space path
    std::vector<Eigen::Vector2i> cs_path;
    cs_path.push_back(q_init_cell);

    // Search for the path
    int safety = 0;
    while (cs_path.back() != q_goal_cell)
    {
        // Get current node
        Eigen::Vector2i c = cs_path.back();
        int& currVal = grid(c.y(), c.x());
        
        // Search neighbors
        for (int k = 0; k < 4; ++k) {
            int nx = c.x() + dx[k];
            int ny = c.y() + dy[k];
            
            // Wrapping if needed
            if (isManipulator)
            {
                nx = (nx + W) % W;
                ny = (ny + H) % H;
            } else {
                if (ny < 0 || ny >= H || nx < 0 || nx >= W)
                {
                    continue;
                }  
            }

            // See if neighbor is closer
            int& cellVal = grid(ny, nx);
            if (cellVal != 1 && cellVal < currVal) { 
                cs_path.push_back(Eigen::Vector2i(nx, ny));
            }
        }

        // Prevent infintie loop
        safety = safety + 1;
        if (safety == 50000)
        {
            std::cout << "Failure";
            break;
        }
        
    }

    // Get configurations for each cell and push to path
    path.waypoints.push_back(q_init);
    for (size_t i = 0; i < cs_path.size(); i++)
    {
        //path.waypoints.push_back(getPointFromCell(cs_path[i].x(), cs_path[i].y(), grid_cspace.x1Bounds().second));
        path.waypoints.push_back(getPointFromCell(cs_path[i].x(), cs_path[i].y(), grid_cspace));
    }
    path.waypoints.push_back(q_goal);
    


    if (isManipulator) {
        Eigen::Vector2d bounds0 = Eigen::Vector2d(0.0, 0.0);
        Eigen::Vector2d bounds1 = Eigen::Vector2d(2*M_PI, 2*M_PI);
        amp::unwrapWaypoints(path.waypoints, bounds0, bounds1);
    }

    return path;
}
