using Revise
using Test
using LEAP
using JLD2, CSV
using Setfield, DataFrames
using Statistics, Distributions, StatsBase, SpecialFunctions
using Distributed

simulation = LEAP.set_up(111,"BC",2001,18,50,"M3");
run_test= LEAP.process(simulation,1,false,true);

