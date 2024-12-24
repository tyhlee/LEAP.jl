# Import required packages
using Asthma_Julia
using JLD, JLD2, CSV
using Setfield, DataFrames
using Statistics, Distributions, StatsBase, SpecialFunctions
using Distributed

# Set up the simulation :
# max_age, province, starting_year, time_horizon, n (proportion, integer, or "full"), population_growth_type
simulation = Asthma_Julia.set_up(111,"CA",2001,30,"full","M3");
# Run the simulation: 
# Takes about ~5 hours on 
# Platform: aarch64-apple-darwin20
# Running under: macOS 15.2
simulation_output = Asthma_Julia.process(simulation,1,false,true);
# Save the output
JLD2.save("internal-validation/LEAP_output.jld","output",simulation_output)
