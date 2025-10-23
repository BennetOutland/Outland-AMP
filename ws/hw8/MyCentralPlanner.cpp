#include "MyMultiAgentPlanners.h"
#include "MySamplingBasedPlanners.h"










// ProblemND makeProblemNDFromMultiAgent2D(const amp::MultiAgentProblem2D& multi_problem) {
//     size_t n_agents = multi_problem.numAgents();
//     ProblemND probND(n_agents);

//     for (size_t i = 0; i < n_agents; ++i) {
//         // Flatten q_init and q_goal
//         probND.q_init.segment<2>(i*2) = multi_problem.agent_properties[i].q_init;
//         probND.q_goal.segment<2>(i*2) = multi_problem.agent_properties[i].q_goal;

//         // Repeat workspace bounds for each agent
//         probND.bounds_min.segment<2>(i*2) << multi_problem.x_min, multi_problem.y_min;
//         probND.bounds_max.segment<2>(i*2) << multi_problem.x_max, multi_problem.y_max;
//     }

//     // Copy the 2D obstacles (broadcasted conceptually to each agent)
//     probND.obstacles = multi_problem.obstacles;

//     return probND;
// }





amp::MultiAgentPath2D MyCentralPlanner::plan(const amp::MultiAgentProblem2D& problem) {
    amp::MultiAgentPath2D path;

    // Determine the max radius 
    double max_r = 0;
    for (size_t i = 0; i < problem.numAgents(); i++)
    {
        if (problem.agent_properties[0].radius > max_r)
        {
            max_r = problem.agent_properties[0].radius;
        }
    }

    MyCentralRRT rrt(0.01, 1.5, 1.0, max_r); // 0.25, 2.0, 2.0 | 0.05, 1.0, 1.0,| 0.025, 1.0, 1.0
    path = rrt.plan(problem);

    return path;
}






amp::Problem2D makeSingleAgentProblem(
    const amp::MultiAgentProblem2D& multi_problem,
    const amp::CircularAgentProperties& agent)
{
    amp::Problem2D p;

    // Copy environment bounds and obstacles
    p.x_min = multi_problem.x_min;
    p.x_max = multi_problem.x_max;
    p.y_min = multi_problem.y_min;
    p.y_max = multi_problem.y_max;
    p.obstacles = multi_problem.obstacles;

    // Assign agent-specific start and goal
    p.q_init = agent.q_init;
    p.q_goal = agent.q_goal;

    return p;
}





amp::MultiAgentPath2D MyDecentralPlanner::plan(const amp::MultiAgentProblem2D& problem) {
    amp::MultiAgentPath2D path;

    // Break into sub-problems and solve
    for (const amp::CircularAgentProperties& agent : problem.agent_properties) {

        amp::Problem2D single_problem = makeSingleAgentProblem(problem, agent);

        
        MyRRT rrt(0.05, 0.5, 0.25, agent.radius, path.agent_paths);
        amp::Path2D agent_path = rrt.plan(single_problem);

        path.agent_paths.push_back(agent_path);
    }


    return path;
}