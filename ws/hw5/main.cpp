/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None
License: NA

Resources:
- 
*/


// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"
#include "hw/HW2.h"

// Include any custom headers you created in your workspace
#include "MyGDAlgorithm.h"

using namespace amp;

int main(int argc, char** argv) {
    /* Include this line to have different randomized environments every time you run your code (NOTE: this has no affect on grade()) */
    amp::RNG::seed(amp::RNG::randiUnbounded());

    // Part a
    // MyGDAlgorithm algo(1.0, 1.0, 0.5, 2.0);
    // Problem2D problem = HW5::getWorkspace1();
    // Path2D path = algo.plan(problem);
    // LOG("path length: " << path.length());
    // Visualizer::makeFigure(problem, path);
    // Visualizer::makeFigure(algo.pot_func, problem, 60);
    // Visualizer::saveFigures();

    // Part b
    // MyGDAlgorithm algo(1.0, 0.25, 1.0, 1.0);
    // Problem2D problem = HW2::getWorkspace1();
    // Path2D path = algo.plan(problem);
    // LOG("path length: " << path.length());
    // Visualizer::makeFigure(problem, path);
    // Visualizer::makeFigure(algo.pot_func, problem, 60);
    // Visualizer::saveFigures();

 
    // Arguments following argv correspond to the constructor arguments of MyGDAlgorithm:
    // MyGDAlgorithm algo(1.0, 1.0, 0.5, 2.0);
    // Path2D path;
    // Problem2D prob;
    // bool success = HW5::generateAndCheck(algo, path, prob);
    // Visualizer::makeFigure(prob, path);

    HW5::grade<MyGDAlgorithm>("bennet.outland@colorado.edu", argc, argv, 1.0, 1.0, 0.25, 0.5);
    return 0;
}