// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW4.h"

// Include the headers for HW4 code
#include "CSpaceSkeleton.h"
#include "ManipulatorSkeleton.h"

// Include the header of the shared class
#include "HelpfulClass.h"

using namespace amp;

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // Part 2
    //MyManipulator2D manipulator(std::vector<double> {0.5, 1.0, 0.5}); // 2a
    //MyManipulator2D manipulator(std::vector<double> {1.5, 1.0, 1.0}); // 2b

    // You can visualize your manipulator given an angle state like so:
    amp::ManipulatorState test_state;
    double pi = 3.14195;
    //test_state = Eigen::Vector3d(pi/6, pi/3, 7*pi/4); // 2a
    //test_state = manipulator.getConfigurationFromIK(Eigen::Vector2d(2, 0)); // 2b

    // The visualizer uses your implementation of forward kinematics to show the joint positions so you can use that to test your FK algorithm
    //Visualizer::makeFigure(manipulator, test_state); 



    // Part 3
    MyManipulator2D manipulator(std::vector<double> {1.0, 1.0});

    // Create the collision space constructor
    std::size_t n_cells = 250;
    MyManipulatorCSConstructor cspace_constructor(n_cells);

    // Create the collision space using a given manipulator and environment
    std::unique_ptr<amp::GridCSpace2D> cspace = cspace_constructor.construct(manipulator, HW4::getEx3Workspace1());

    // You can visualize your cspace 
    Visualizer::makeFigure(*cspace);

    Visualizer::saveFigures();

    // Grade method
    // amp::HW4::grade<MyManipulator2D>(cspace_constructor, "nonhuman.biologic@myspace.edu", argc, argv);
    return 0;
}