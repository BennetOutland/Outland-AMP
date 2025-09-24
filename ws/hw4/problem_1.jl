"""
Author: Bennet Outland
Affiliation: University of Colorado Boulder
Control: None
License: NA

Resouces:
- LazySets Docs: https://juliareach.github.io/LazySets.jl/dev/
- ChatGPT for need to overapproximate
- ChatGPT for type debugging assistance
- ChatGPT for plotting debugging
"""

# Usings
using LazySets
using Plots
using LinearAlgebra
using Polyhedra
using GLMakie

# Define the obstacle and robot sets 
O = VPolygon([[0.0, 0.0], [1.0, 2.0], [0.0, 2.0]])
A = VPolygon([[0.0, 0.0], [1.0, 2.0], [0.0, 2.0]])
A_reflected = LinearMap(-I, A)


# Part A 
C_obs_lazy = MinkowskiSum(O, A_reflected)
C_obs = overapproximate(C_obs_lazy, VPolygon)
vertices(C_obs)
# plot(C_obs)


# Part B 

"""
2D Rotation matrix given an angle
"""
function Rot(θ)
    return [cos(θ) -sin(θ); sin(θ) cos(θ)]
end

# Loo through and collect obstacles
C_obs_vec = []
dθ = 2π/11
θ_vec = 0.0:dθ:2π

for θ ∈ θ_vec
    local A_rot = LinearMap(Rot(θ), A_reflected)
    local C_obs_lazy = MinkowskiSum(O, A_rot)
    local C_obs = overapproximate(C_obs_lazy, VPolygon)

    push!(C_obs_vec, C_obs)
end


# - Make a new polygon with the angle as third coordinate
# - CH two polygons
# - union the hulls together
obs_3d = [] 

for i ∈ 1:length(C_obs_vec)-1
    # Bottom plane
    vs_1 = vertices(C_obs_vec[i])
    for (j, v1) ∈ enumerate(vs_1)
        vs_1[j] = [v1[1], v1[2], θ_vec[i]]
    end

    # Top plane 
    vs_2 = vertices(C_obs_vec[i+1])
    for (j, v2) ∈ enumerate(vs_2)
        vs_2[j] = [v2[1], v2[2], θ_vec[i+1]]
    end

    # Collect into a singular object
    vs_collected = [vs_1..., vs_2...]

    # Make a convex hull polygon
    hull = convex_hull(vs_collected)
    poly_hull = VPolytope(hull)

    push!(obs_3d, poly_hull)
end



# Create the figure
fig = Figure()
ax = Axis3(fig[1, 1], xlabel="x", ylabel="y", zlabel="θ")
for i in 1:length(obs_3d)
    LazySets.plot3d!(obs_3d[i])  
end
fig
