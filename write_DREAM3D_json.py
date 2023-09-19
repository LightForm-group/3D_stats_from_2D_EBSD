import os
import json
import numpy as np
import sys

# Read the input JSON file
with open('ebsd_parameters.json', 'r') as f:
    data = json.load(f)

# Access the variables
num_bins = np.array(data['num_bins'])
bin_number = np.array(data['bin_number'])
feature_diameter_info = np.array(data['feature_diameter_info'])

mu_grain_size = np.array(data['mu_grain_size'])
sigma_grain_size = np.array(data['sigma_grain_size'])

alpha_b_over_a = np.array(data['alpha_b_over_a'])
beta_b_over_a = np.array(data['beta_b_over_a'])
alpha_c_over_a = np.array(data['alpha_c_over_a'])
beta_c_over_a = np.array(data['beta_c_over_a'])

mu_neighbors = np.array(data['mu_neighbors'])
sigma_neighbors = np.array(data['sigma_neighbors'])

Axis_ODF = np.array(data['Axis_ODF'])
sampled_odf = np.array(data['sampled_odf'])

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
data["00"]["StatsDataArray"]["1"]["FeatureSize Distribution"]["Average"] = float(mu_grain_size)
data["00"]["StatsDataArray"]["1"]["FeatureSize Distribution"]["Standard Deviation"] = float(sigma_grain_size)

# Set B/A & C/A distributions
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs B Over A Distributions"]["Alpha"] = list(alpha_b_over_a)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs B Over A Distributions"]["Beta"] = list(beta_b_over_a)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs C Over A Distributions"]["Alpha"] = list(alpha_c_over_a)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs C Over A Distributions"]["Beta"] = list(beta_c_over_a)

# Neighbor distributions
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs Neighbors Distributions"]["Average"] = list(mu_neighbors)
data["00"]["StatsDataArray"]["1"]["FeatureSize Vs Neighbors Distributions"]["Standard Deviation"] = list(sigma_neighbors)

# Axis ODF - this seems to have issues and so 
if Axis_ODF.ndim == 1:
    Axis_ODF = Axis_ODF.reshape(5, 1)
    data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Euler 1"] = list(Axis_ODF[0])
    data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Euler 2"] = list(Axis_ODF[1])
    data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Euler 3"] = list(Axis_ODF[2])
    data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Sigma"  ] = list(Axis_ODF[3])
    data["00"]["StatsDataArray"]["1"]["AxisODF-Weights"]["Weight" ] = list(Axis_ODF[4])
else:
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
data["00"]["StatsDataArray"]["1"]["Feature_Diameter_Info"] = list(feature_diameter_info)
data["00"]["StatsDataArray"]["1"]["BinNumber"] = list(bin_number)
data["00"]["StatsDataArray"]["1"]["Bin Count"] = int(num_bins)

# Output the files to the current working directory
ouput_path = cwd
data["10"]["OutputPath"] = ouput_path
data["10"]["OutputFilePath"] = ouput_path

def numpy_to_native(d):
    if isinstance(d, dict):
        return {key: numpy_to_native(value) for key, value in d.items()}
    elif isinstance(d, (list, tuple)):
        return [numpy_to_native(item) for item in d]
    elif isinstance(d, np.ndarray):
        return d.tolist()
    elif isinstance(d, (np.int32, np.int64, np.float32, np.float64)):
        return d.item()
    else:
        return d

# Convert NumPy types to native types
data_native = numpy_to_native(data)

# write updates to file and close it
with open('EBSD2RVE_pancake.json', 'w') as outfile:
    json.dump(data_native, outfile)