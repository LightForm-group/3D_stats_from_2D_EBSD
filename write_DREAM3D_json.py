import os
import json
import numpy as np
import sys

# Load in variables from text files
num_bins                = np.loadtxt('temp/num_bins.txt',delimiter=',')
bin_number              = np.loadtxt('temp/bin_number.txt',delimiter=',')
Feature_Diameter_Info   = np.loadtxt('temp/Feature_Diameter_Info.txt',delimiter=',')
mu_grain_size           = np.loadtxt('temp/mu_grain_size.txt',delimiter=',')
sigma_grain_size        = np.loadtxt('temp/sigma_grain_size.txt',delimiter=',')
alpha_b_over_a          = np.loadtxt('temp/alpha_b_over_a.txt',delimiter=',')
beta_b_over_a           = np.loadtxt('temp/beta_b_over_a.txt',delimiter=',')
alpha_c_over_a          = np.loadtxt('temp/alpha_c_over_a.txt',delimiter=',')
beta_c_over_a           = np.loadtxt('temp/beta_c_over_a.txt',delimiter=',')
mu_neighbors            = np.loadtxt('temp/mu_neighbors.txt',delimiter=',')
sigma_neighbors         = np.loadtxt('temp/sigma_neighbors.txt',delimiter=',')
Axis_ODF                = np.loadtxt('temp/Axis_ODF.txt',delimiter=',')
sampled_odf             = np.loadtxt('temp/sampled_odf.txt',delimiter=',')

# Set current working directory
cwd = os.getcwd() 

# Edit the json file containing creation of the microstructure
# Opening JSON file
f = open('EBSD2RVE.json',)
 
# returns JSON object as a dictionary
data = json.load(f)

# Close the file
f.close()

# Set grain size, it is stored in a log format in the json file
data["00"]["StatsDataArray"]["1"]["FeatureSize Distribution"]["mu"] = float(mu_grain_size)
data["00"]["StatsDataArray"]["1"]["FeatureSize Distribution"]["Standard Deviation"] = float(sigma_grain_size)

# Set B/A & C/A distributions
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs B Over A Distributions"]["Alpha"] = list(alpha_b_over_a)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs B Over A Distributions"]["Beta"] = list(beta_b_over_a)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs C Over A Distributions"]["Alpha"] = list(alpha_c_over_a)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs C Over A Distributions"]["Beta"] = list(beta_c_over_a)

# Neighbor distributions
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs Neighbors Distributions"]["Average"] = list(mu_neighbors)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs Neighbors Distributions"]["Standard Deviation"] = list(sigma_neighbors)

# Axis ODF
data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Euler 1"] = list(Axis_ODF[:,0])
data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Euler 2"] = list(Axis_ODF[:,1])
data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Euler 3"] = list(Axis_ODF[:,2])
data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Sigma"  ] = list(Axis_ODF[:,3])
data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Weight" ] = list(Axis_ODF[:,4])

# ODF
data["00"]["StatsDataArray"]["1"]["ODF-Weights"]["Euler 1"] = list(sampled_odf[:,0])
data["00"]["StatsDataArray"]["1"]["ODF-Weights"]["Euler 2"] = list(sampled_odf[:,1])
data["00"]["StatsDataArray"]["1"]["ODF-Weights"]["Euler 3"] = list(sampled_odf[:,2])
data["00"]["StatsDataArray"]["1"]["ODF-Weights"]["Sigma"  ] = list(sampled_odf[:,3])
data["00"]["StatsDataArray"]["1"]["ODF-Weights"]["Weight" ] = list(sampled_odf[:,4])

# Binning
data["00"]["StatsDataArray"]["1"]["Feature_Diameter_Info"] = list(Feature_Diameter_Info)
data["00"]["StatsDataArray"]["1"]["BinNumber"] = list(bin_number)
data["00"]["StatsDataArray"]["1"]["Bin Count"] = float(num_bins)

# Output the files to the current working directory
ouput_path = cwd
data["10"]["OutputPath"] = ouput_path
data["10"]["OutputFilePath"] = ouput_path

# write updates to file and close it
with open('EBSD2RVE.json', 'w') as outfile:
    json.dump(data, outfile)