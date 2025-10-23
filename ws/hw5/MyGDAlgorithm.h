/*
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None
License: NA

Resources:
- Lecture Slides
*/

#pragma once

// This includes all of the necessary header files in the toolbox
#include "AMPCore.h"

// Include the correct homework header
#include "hw/HW5.h"



class MyPotentialFunction : public amp::PotentialFunction2D {
    public:
		// Initialize the field with the passed in params
		MyPotentialFunction(double d_star, double zeta, double Q_star, double eta);

		// Adds the boundary conditions for the trajectory
		void add_boundary_values(amp::Problem2D problem);

		// Returns the potential function value (height) for a given 2D point. 
		double operator()(const Eigen::Vector2d& q) const;

		// Returns the gradient at a point
		Eigen::Vector2d getGradient(const Eigen::Vector2d& q) const;

		double d_star, zeta, Q_star, eta;
		Eigen::Vector2d q_init, q_goal;
		amp::Problem2D problem;
};





class MyGDAlgorithm : public amp::GDAlgorithm {
	public:
		// Constructor
		MyGDAlgorithm(double d_star, double zeta, double Q_star, double eta);

		// Override this method to solve a given problem.
		virtual amp::Path2D plan(const amp::Problem2D& problem) override;
		// private:
		double d_star, zeta, Q_star, eta, alpha;
		MyPotentialFunction pot_func;
		// Add additional member variables here...
};

