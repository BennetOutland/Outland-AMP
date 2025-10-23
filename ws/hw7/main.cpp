/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None

Resources:
- Used ChatGPT to help set up the trials
- Collaborated with Emily Maxwell Outland on the benchmarking code
*/



// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"
#include "hw/HW2.h"
#include "hw/HW5.h"
#include "MySamplingBasedPlanners.h"

using namespace amp;

#include <string>
#include <chrono>


int main(int argc, char** argv) {
    //HW7::hint(); // Consider implementing an N-dimensional planner 

    // Example of creating a graph and adding nodes for visualization
    std::shared_ptr<amp::Graph<double>> graphPtr = std::make_shared<amp::Graph<double>>();
    std::map<amp::Node, Eigen::Vector2d> nodes;
    

    // Test PRM on Workspace1 of HW2
    // Problem2D problem = HW2::getWorkspace1();
    // MyPRM prm(1000, 1.0, true);
    // Path2D path = prm.plan(problem);
    // std::cout << prm.plan_success;
    // LOG("path length: " << path.length());
    // Visualizer::makeFigure(problem, path, *prm.graphPtr, prm.nodes);
    // Visualizer::saveFigures();



    // // Define the (n, r) test pairs
    // std::vector<std::pair<int, double>> param_pairs = {
    //     {200, 0.5}, {200, 1.0}, {200, 1.5}, {200, 2.0},
    //     {500, 0.5}, {500, 1.0}, {500, 1.5}, {500, 2.0}
    // };
    std::vector<std::pair<int, double>> param_pairs = {
        {200, 1.0}, {200, 2.0}, {500, 1.0}, {500, 2.0},
        {1000, 1.0}, {1000, 2.0}
    };
    std::list<std::vector<double>> time_all_runtimes;
    std::list<std::vector<double>> valid_all_runtimes;
    std::list<std::vector<double>> path_all_runtimes;
    std::vector<std::string> labels;


    // Loop over all (n, r) parameter combinations
    for (const auto& [n, r] : param_pairs) {
        std::vector<double> time_runtimes;
        time_runtimes.reserve(100);

        std::vector<double> valid_runtimes;
        valid_runtimes.reserve(100);

        std::vector<double> path_runtimes;
        path_runtimes.reserve(100);

        for (size_t i = 0; i < 100; ++i) {
            Problem2D problem = HW2::getWorkspace2();
            MyPRM prm(n, r, true);

            // Timing
            auto start = std::chrono::high_resolution_clock::now();
            Path2D path = prm.plan(problem);
            auto end = std::chrono::high_resolution_clock::now();
            double elapsed_ms = duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
            time_runtimes.push_back(elapsed_ms);

            // Valid
            try
            {
                valid_runtimes.push_back(static_cast<double>(prm.plan_success));
            }
            catch(...)
            {
                valid_runtimes.push_back(0.0);
            }
            

            // Path length 
            try
            {
                path_runtimes.push_back(path.length());
            }
            catch(const std::exception& e)
            {
                path_runtimes.push_back(0.0);
            }
            

        }

        // Make the label
        std::ostringstream label_stream;
        label_stream << "(" << n << "," << std::fixed << std::setprecision(1) << r << ")";
        std::string label = label_stream.str();
        labels.push_back(label);

        // Push the runtimes
        time_all_runtimes.push_back(time_runtimes);
        valid_all_runtimes.push_back(valid_runtimes);
        path_all_runtimes.push_back(path_runtimes);


    }


    // Plot the box plots
    amp::Visualizer::makeBoxPlot(
        time_all_runtimes, 
        labels, 
        "PRM Benchmark Runtimes", 
        "(n, r) Parameters", 
        "Computation Time [ms]"
    );

    amp::Visualizer::makeBoxPlot(
        path_all_runtimes, 
        labels, 
        "PRM Benchmark Runtimes", 
        "(n, r) Parameters", 
        "Path Length"
    );

    amp::Visualizer::makeBoxPlot(
        valid_all_runtimes, 
        labels, 
        "PRM Benchmark Runtimes", 
        "(n, r) Parameters", 
        "Valid Solutions"
    );



    //Generate a random problem and test RRT
    // Problem2D problem = HW2::getWorkspace2();
    // MyRRT rrt(0.05, 0.5, 0.25);
    // Path2D path = rrt.plan(problem);
    // LOG("path length: " << path.length());
    // Visualizer::makeFigure(problem, path, *rrt.graphPtr, rrt.nodes);



    // // Problem vector 
    // std::vector<Problem2D> problem_vec = {HW5::getWorkspace1(), HW2::getWorkspace1(), HW2::getWorkspace2()};
    // std::list<std::vector<double>> time_all_runtimes;
    // std::list<std::vector<double>> valid_all_runtimes;
    // std::list<std::vector<double>> path_all_runtimes;
    // std::vector<std::string> labels = {"HW2W1", "HW2W2", "HW5W1"};

    // // Loop over all (n, r) parameter combinations
    // for (const auto& prob: problem_vec) {
    //     std::vector<double> time_runtimes;
    //     time_runtimes.reserve(100);

    //     std::vector<double> valid_runtimes;
    //     valid_runtimes.reserve(100);

    //     std::vector<double> path_runtimes;
    //     path_runtimes.reserve(100);

    //     for (size_t i = 0; i < 100; ++i) {
    //         MyRRT rrt(0.05, 0.5, 0.25);

    //         // Timing
    //         auto start = std::chrono::high_resolution_clock::now();
    //         Path2D path = rrt.plan(prob);
    //         auto end = std::chrono::high_resolution_clock::now();
    //         double elapsed_ms = duration_cast<std::chrono::microseconds>(end - start).count() / 1000.0;
    //         time_runtimes.push_back(elapsed_ms);

    //         // Valid
    //         try
    //         {
    //             valid_runtimes.push_back(static_cast<double>(rrt.plan_success));
    //         }
    //         catch(...)
    //         {
    //             valid_runtimes.push_back(0.0);
    //         }
            

    //         // Path length 
    //         try
    //         {
    //             path_runtimes.push_back(path.length());
    //         }
    //         catch(const std::exception& e)
    //         {
    //             path_runtimes.push_back(0.0);
    //         }
            

    //     }


    //     // Push the runtimes
    //     time_all_runtimes.push_back(time_runtimes);
    //     valid_all_runtimes.push_back(valid_runtimes);
    //     path_all_runtimes.push_back(path_runtimes);

    // }


    // // Plot the box plots
    // amp::Visualizer::makeBoxPlot(
    //     time_all_runtimes, 
    //     labels, 
    //     "RRT Benchmark Runtimes", 
    //     "Problems", 
    //     "Computation Time [ms]"
    // );

    // amp::Visualizer::makeBoxPlot(
    //     path_all_runtimes, 
    //     labels, 
    //     "RRT Benchmark Runtimes", 
    //     "Problems", 
    //     "Path Length"
    // );

    // amp::Visualizer::makeBoxPlot(
    //     valid_all_runtimes, 
    //     labels, 
    //     "RRT Benchmark Runtimes", 
    //     "Problems", 
    //     "Valid Solutions"
    // );




    Visualizer::saveFigures();


    //HW7::generateAndCheck(rrt, path, problem);

    // Grade method
    //HW7::grade<MyPRM, MyRRT>("bennet.outland@colorado.edu", argc, argv, std::make_tuple(1000, 2.0, true), std::make_tuple(0.05, 2.0, 0.35));
    return 0;
}