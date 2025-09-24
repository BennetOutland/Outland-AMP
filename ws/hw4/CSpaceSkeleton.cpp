/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- ChatGPT misc debugging
*/

#include <Eigen/Dense>
#include <Eigen/Core>
#include "CSpaceSkeleton.h"

// Override this method for returning whether or not a point is in collision

std::pair<std::size_t, std::size_t> MyGridCSpace2D::getCellFromPoint(double x0, double x1) const {
    // Implment your discretization procedure here, such that the point (x0, x1) lies within the returned cell
    std::size_t cell_x = 0; // x index of cell
    std::size_t cell_y = 0; // y index of cell
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

// Override this method for computing all of the boolean collision values for each cell in the cspace
std::unique_ptr<amp::GridCSpace2D> MyManipulatorCSConstructor::construct(const amp::LinkManipulator2D& manipulator, const amp::Environment2D& env) {
    // Create an object of my custom cspace type (e.g. MyGridCSpace2D) and store it in a unique pointer. 
    // Pass the constructor parameters to std::make_unique()
    std::unique_ptr<MyGridCSpace2D> cspace_ptr = std::make_unique<MyGridCSpace2D>(m_cells_per_dim, m_cells_per_dim, env.x_min, env.x_max, env.y_min, env.y_max);
    // In order to use the pointer as a regular GridCSpace2D object, we can just create a reference
    MyGridCSpace2D& cspace = *cspace_ptr;

    // Determine if each cell is in collision or not, and store the values the cspace. This `()` operator comes from DenseArray base class

    // Get a point from the cell
    float dx = (env.x_max - env.x_min) / (m_cells_per_dim);
    float dy = (env.y_max - env.y_min) / (m_cells_per_dim);

    // Check if collision
    float x_i = env.x_min + (dx / 2.0);
    float y_i = env.y_min + (dy / 2.0);

    // Line increment parameter 
    int divisions = 50;

    // Loop through each grid point
    int kx = 0;
    int ky = 0;
    while (x_i < env.x_max)
    {
        // Loop trhough
        while (y_i < env.y_max)
        {
            // Collect the cell point
            Eigen::Vector2d pt = Eigen::Vector2d(x_i, y_i);

            // Get the coordinates of the arm points
            Eigen::Vector2d j_1 = manipulator.getJointLocation(pt, 1);
            Eigen::Vector2d j_2 = manipulator.getJointLocation(pt, 2);

            // Determine if they are in any obstacles 
            bool collision = false;
            for (size_t i = 0; i < env.obstacles.size(); i++)
            {
                // Create first joint line 
                for (size_t l = 0; l <= divisions; l++)
                {
                    // Test points along the joint
                    double t = static_cast<double>(l) / divisions;
                    Eigen::Vector2d test_pt = (t) * Eigen::Vector2d(0, 0) + (1 - t) * j_1;
                    if (pip(env.obstacles[i].verticesCCW(), test_pt)) {
                        collision = true;
                        break;
                    }
                }
                if (collision)
                {
                    break;
                }
                
                
                // Create second joint line 
                            // Create first joint line 
                for (size_t l = 0; l <= divisions; l++)
                {
                    // Test points along the joint
                    double t = static_cast<double>(l) / divisions;
                    Eigen::Vector2d test_pt = (t) * j_1 + (1 - t) * j_2;
                    if (pip(env.obstacles[i].verticesCCW(), test_pt)) {
                        collision = true;
                        break;
                    }
                }
                if (collision)
                {
                    break;
                }
             
            }
            
            // Update the c-space
            cspace(kx, ky) = collision;

            // Update y info 
            y_i = y_i + dy;
            ky += 1;
        }

        // Reset y_i 
        y_i = env.y_min + (dy / 2.0);
        ky = 0;


        // Update x info 
        x_i = x_i + dx;
        kx += 1;
    }
    

    // Returning the object of type std::unique_ptr<MyGridCSpace2D> can automatically cast it to a polymorphic base-class pointer of type std::unique_ptr<amp::GridCSpace2D>.
    // The reason why this works is not super important for our purposes, but if you are curious, look up polymorphism!
    return cspace_ptr;
}
