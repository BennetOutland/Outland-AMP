/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- ChatGPT for randomized testing
*/



// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW4.h"

// Include the headers for HW4 code
#include "CSpaceSkeleton.h"
#include "ManipulatorSkeleton.h"

// Include the header of the shared class
#include "HelpfulClass.h"
#include <random>

using namespace amp;



Eigen::Vector2d randomVector2dInDisk(double max_radius = 3.5) {
    static std::random_device rd;
    static std::mt19937 gen(rd());
    static std::uniform_real_distribution<> dist_angle(0.0, 2.0 * M_PI);
    static std::uniform_real_distribution<> dist_radius(0.0, 1.0);

    double theta = dist_angle(gen);
    double r = max_radius * std::sqrt(dist_radius(gen)); // sqrt ensures uniform distribution in area

    return Eigen::Vector2d(r * std::cos(theta), r * std::sin(theta));
}

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // Part 2
    //MyManipulator2D manipulator(std::vector<double> {0.5, 1.0, 0.5}); // 2a
    // MyManipulator2D manipulator(std::vector<double> {1.5, 1.0, 1.0}); // 2b

    // // You can visualize your manipulator given an angle state like so:
    // amp::ManipulatorState test_state;
    // double pi = 3.14195;
    // //test_state = Eigen::Vector3d(pi/6, pi/3, 7*pi/4); // 2a
    // test_state = manipulator.getConfigurationFromIK(Eigen::Vector2d(2.0, 0.0)); // 2b

    // // The visualizer uses your implementation of forward kinematics to show the joint positions so you can use that to test your FK algorithm
    // Visualizer::makeFigure(manipulator, test_state); 
    // Visualizer::saveFigures();

    // Testing suite
    // int passes = 0;
    // int iters = 1000;
    // for (size_t i = 0; i < iters; i++)
    // {
    //     MyManipulator2D manipulator(std::vector<double> {1.0, 1.0});
    //     Eigen::Vector2d pt = randomVector2dInDisk(2.0);
    //     amp::ManipulatorState test_state = manipulator.getConfigurationFromIK(pt);

    //     if ((pt - manipulator.getJointLocation(test_state, 2)).norm() <= 0.0001)
    //     {
    //         passes += 1;
    //     }
        
    // }

    // std::cout << "Results \n";
    // std::cout << passes;

    // 2 link FK test 
    // MyManipulator2D manipulator(std::vector<double> {1.0, 1.0});
    // amp::ManipulatorState test_state = Eigen::Vector2d(pi/2, 0.0);
    // Visualizer::makeFigure(manipulator, test_state); 
    // Visualizer::saveFigures();
    //amp::ManipulatorState test_state = manipulator.getConfigurationFromIK();
    



    // Part 3
    // MyManipulator2D manipulator(std::vector<double> {1.0, 1.0});

    // // Create the collision space constructor
    // std::size_t n_cells = 250;
    // MyManipulatorCSConstructor cspace_constructor(n_cells);

    // // Create the collision space using a given manipulator and environment
    // std::unique_ptr<amp::GridCSpace2D> cspace = cspace_constructor.construct(manipulator, HW4::getEx3Workspace2());

    // // You can visualize your cspace 
    // //Visualizer::makeFigure(*cspace);

    // // Visulaize the workspace 
    // Visualizer::makeFigure(HW4::getEx3Workspace3());

    // Visualizer::saveFigures();

    // Grade method
    MyManipulator2D manipulator();
    std::size_t n_cells = 250; // 250
    MyManipulatorCSConstructor cspace_constructor(n_cells);

    amp::HW4::grade<MyManipulator2D>(cspace_constructor, "bennet.outland@colorado.edu", argc, argv);
    return 0;
}