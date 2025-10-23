/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- ChatGPT misc debugging
*/

#include <Eigen/Dense>
#include <algorithm>
#include <Eigen/Core>
#include "CSpaceSkeleton.h"

// Override this method for returning whether or not a point is in collision

std::pair<std::size_t, std::size_t> MyGridCSpace2D::getCellFromPoint(double x0, double x1) const {
    // Hard code x_min.x_max etc. to be in the c space not the workspace
    double pi = 3.14159;
    double x_min = 0;
    double x_max = 2*pi;
    double y_min = 0;
    double y_max = 2*pi;
    float m_cells_per_dim = 250.0;

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

    // TODO: hard code x_min.x_max etc. to be in the c space not the workspace
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