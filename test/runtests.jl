using Revise
using Test
using JLD, JLD2, CSV
using Setfield, DataFrames
using Statistics, Distributions, StatsBase, SpecialFunctions
using Distributed

simulation = Asthma_Julia.set_up(111,"CA",2001,40,100,"M3");
run_test= Asthma_Julia.process(simulation,1,false,true);