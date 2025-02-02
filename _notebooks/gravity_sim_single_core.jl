using Pkg
Pkg.activate(".")

using LinearAlgebra
using Plots


# Structure to represent a point mass
mutable struct Body
    position::Vector{Float64}  # [x, y]
    velocity::Vector{Float64}  # [vx, vy]
    mass::Float64
    force::Vector{Float64}     # [fx, fy]
    radius::Float64           # For collision detection
    active::Bool             # Track if body is still active (not merged)
end

# Quadtree node structure
mutable struct QuadTreeNode
    center::Vector{Float64}    # Center of the quadrant
    size::Float64             # Side length of the quadrant

    mass::Float64             # Total mass in this node
    com::Vector{Float64}      # Center of mass
    body::Union{Body, Nothing} # Single body (for leaf nodes)

    children::Vector{Union{QuadTreeNode, Nothing}}  # Four children nodes
    isLeaf::Bool              # Is this a leaf node?

    # Constructor for empty node
    function QuadTreeNode(center::Vector{Float64}, size::Float64)
        new(center, size, 0.0, [0.0, 0.0], nothing, 
            [nothing, nothing, nothing, nothing], true)
    end
end

# Constants
const G = 1.0  # Gravitational constant

const θ = 0.5  # Barnes-Hut threshold (smaller = more accurate)
const dt = 0.1 # Time step

# Insert a body into the quadtree
function insert!(node::QuadTreeNode, body::Body)
    if !body.active
        return  # Don't insert inactive bodies
    end
    
    if node.isLeaf && node.body === nothing
        # Empty leaf node - add the body here
        node.body = body
        node.mass = body.mass
        node.com = body.position
        return
    end
    
    if node.isLeaf && node.body !== nothing

        # Need to split this node
        old_body = node.body
        node.body = nothing
        node.isLeaf = false
        
        # Create four children
        size = node.size / 2
        for i in 1:4
            dx = (i in [2, 4]) ? size/2 : -size/2
            dy = (i in [3, 4]) ? size/2 : -size/2
            center = node.center + [dx, dy]
            node.children[i] = QuadTreeNode(center, size)
        end
        
        # Insert the old body into appropriate child
        insert!(node, old_body)
    end
    
    # Update mass and center of mass
    node.mass += body.mass
    node.com = (node.com * node.mass + body.position * body.mass) / (node.mass + body.mass)
    
    # Find appropriate quadrant and insert
    if !node.isLeaf
        quad = get_quadrant(node, body.position)
        insert!(node.children[quad], body)
    end
end


# Helper function to determine which quadrant a position belongs to
function get_quadrant(node::QuadTreeNode, pos::Vector{Float64})
    if pos[1] >= node.center[1]
        return pos[2] >= node.center[2] ? 4 : 2
    else
        return pos[2] >= node.center[2] ? 3 : 1
    end
end

# Calculate forces on a body using Barnes-Hut approximation
function calculate_force!(node::QuadTreeNode, body::Body)
    if !body.active || node.body === body || node.mass == 0
        return
    end
    

    d = norm(node.com - body.position)
    if d == 0
        return
    end
    
    if node.isLeaf || (node.size / d < θ)
        # Either a leaf node or far enough to approximate
        direction = normalize(node.com - body.position)

        magnitude = G * body.mass * node.mass / (d * d)

        body.force += direction * magnitude

    else
        # Node is too close, need to check children
        for child in node.children
            if child !== nothing
                calculate_force!(child, body)
            end
        end
    end
end

# Handle collision between two bodies (perfectly inelastic collision)
function merge_bodies!(body1::Body, body2::Body)
    # Calculate new mass
    new_mass = body1.mass + body2.mass
    

    # Calculate new position (center of mass)
    new_position = (body1.position * body1.mass + body2.position * body2.mass) / new_mass
    
    # Calculate new velocity (conservation of momentum)
    new_velocity = (body1.velocity * body1.mass + body2.velocity * body2.mass) / new_mass
    
    # Calculate new radius (assuming density remains constant)

    new_radius = (body1.radius^3 + body2.radius^3)^(1/3)
    
    # Update primary body

    body1.mass = new_mass
    body1.position = new_position
    body1.velocity = new_velocity
    body1.radius = new_radius

    
    # Deactivate secondary body
    body2.active = false
end

# Check and handle collisions between all active bodies
function handle_collisions!(bodies::Vector{Body})
    n = length(bodies)
    for i in 1:n
        if !bodies[i].active
            continue
        end
        for j in (i+1):n
            if !bodies[j].active
                continue
            end
            
            # Calculate distance between bodies
            distance = norm(bodies[i].position - bodies[j].position)
            
            # Check for collision
            if distance < (bodies[i].radius + bodies[j].radius)

                merge_bodies!(bodies[i], bodies[j])
            end
        end
    end
end

# Update positions and velocities using basic Euler integration
function update!(bodies::Vector{Body}, dt::Float64)
    for body in bodies
        if !body.active
            continue
        end
        
        # Update velocity and position using current force
        body.velocity += body.force / body.mass * dt
        body.position += body.velocity * dt
        
        # Reset force for next iteration
        body.force = [0.0, 0.0]
    end
end

# Modified simulation function with visualization
function simulate_with_vis(n_bodies::Int, n_steps::Int)
    # Initialize random bodies
    bodies = [Body(
        rand(2) * 100,           # Random position in [0, 100]
        (rand(2) .- 0.5) * 10,   # Random velocity in [-5, 5]
        rand() * 10 + 1,         # Random mass in [1, 11]
        [0.0, 0.0],             # Initial force
        (rand() * 2 + 1),       # Random radius between 1 and 3
        true                    # Initially active
    ) for _ in 1:n_bodies]
    
    # Setup plot
    gr()  # Use GR backend for better performance
    
    # Create color gradient for masses
    color_range = range(0, 1, length=100)
    
    # Animation loop
    @animate for step in 1:n_steps
        # Create root node of quadtree
        root = QuadTreeNode([50.0, 50.0], 100.0)
        
        # Build tree (only with active bodies)
        for body in bodies
            insert!(root, body)
        end
        
        # Calculate forces
        for body in bodies

            calculate_force!(root, body)
        end
        

        # Handle collisions
        handle_collisions!(bodies)

        
        # Update positions and velocities

        update!(bodies, dt)
        
        # Extract data for plotting
        active_bodies = filter(b -> b.active, bodies)
        positions = getindex.(getfield.(active_bodies, :position), [1, 2])
        x = getindex.(positions, 1)
        y = getindex.(positions, 2)
        
        # Get masses and radii for active bodies
        masses = getfield.(active_bodies, :mass)

        radii = getfield.(active_bodies, :radius)

        
        # Normalize masses for color mapping
        max_mass = maximum(masses)
        normalized_masses = masses ./ max_mass
        
        # Create plot
        scatter(x, y,
            marker_z=masses,  # Color by mass
            color=:thermal,
            markersize=radii .* 2,  # Scale points by radius
            legend=false,
            xlim=(0, 100),

            ylim=(0, 100),
            title="N-body Simulation\nStep: $step",
            xlabel="X",
            ylabel="Y",
            aspect_ratio=:equal,
            background_color=:black,
            grid=false,
            framestyle=:box
        )
    end every 5  # Only create animation frame every 5 steps for better performance
end

# Function to run simulation and save animation
function run_simulation(n_bodies::Int, n_steps::Int, filename::String="nbody_sim.gif")
    anim = simulate_with_vis(n_bodies, n_steps)
    gif(anim, filename, fps=30)
end

# Example usage:
n_bodies = 100
n_steps = 500
run_simulation(n_bodies, n_steps)
