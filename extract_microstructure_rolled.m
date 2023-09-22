%% Requirments
% The Statistics and Machine Learning Toolbox
% Automatic Saltykov - https://uk.mathworks.com/matlabcentral/fileexchange/73726-automatic-saltykov
close all
clear 
clc

%% Specify Dream3D statistics generation inputs
% For the size distribution - used for binning as described here
% (http://127.0.0.1:32456/Filters/SyntheticBuildingFilters/StatsGeneratorFilter/)
bin_size = 25;
min_sigma_cutoff = 3;
max_sigma_cutoff = 2;

%% Specify Crystal, Grain Structure and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98]),...
  'notIndexed'};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

%% Load in Data 
% NDRD - fibrous structure
fname       = 'test_data\Fibrous model\ND\ND.ctf'; % Path to the NDRD EBSD map
grains_ND  = load_ebsd_data(fname,CS);

% RDTD - equiaxed structure
fname       = 'test_data\Fibrous model\TD\TD.ctf'; % Path to the RDTD EBSD map
grains_TD = load_ebsd_data(fname,CS);

%% Process 
% For rolled grains we can assume that they are roughly pancake shaped
% cylinders with two EBSD maps taken, one orthogonal to the short axis and 
% one with the short axis in the plane

% Find the axis and centroids of the ND structure - this one will look
% more equiaxed
[~,a,b] = grains_ND.fitEllipse; % coincides with the actual grain area
[x,y]           = centroid(grains_ND); % find the grain centres
C_full_ND      = cat(2,x,y); % Combine

% Find the axis and centroids of the TD structure - this one will look
% more fibrous
[~,~,c] = grains_TD.fitEllipse; % coincides with the actual grain area
[x,y]             = centroid(grains_TD); % find the grain centres
C_full_TD       = cat(2,x,y); % Combine

% Find the distribution of the equivelant radius
% Here the volume is found from taking the short axis of the RDTD and randomly 
% assigning it to the long axis from the NDRD distribution while keeping the 
% distribution consistent. this is subsequently used for further 
% calculations.
sampled_c = exp(normrnd(mean(log(c)), std(log(c)), [1 size(b,1)]))';

% These values then need to be reordered and potentially changed to ensure
% the condition of a>b>c is maintained

% Sort b and sampled_c while keeping track of the original indices
[sorted_b, original_indices] = sort(b);
sorted_sampled_c = sort(sampled_c);

% Find the elements in sorted_sampled_c that are greater than corresponding elements in sorted_b
invalid_indices = find(sorted_sampled_c > sorted_b, 1);

% If any such elements are found, replace them and all subsequent elements with corresponding sorted_b values
if ~isempty(invalid_indices)
    sorted_sampled_c(invalid_indices:end) = sorted_b(invalid_indices:end);
end

% Reorder sorted_sampled_c back to match the original sequence of b
sampled_c(original_indices) = sorted_sampled_c;

% We will also need to create a list of sampled_a and sampled_b
sampled_a = exp(normrnd(mean(log(a)), std(log(a)), [1 size(c,1)]))';
sampled_b = exp(normrnd(mean(log(b)), std(log(b)), [1 size(c,1)]))';

% Sort sampled_b and c while keeping track of the original indices
[sorted_c, original_indices] = sort(c);
sorted_sampled_b = sort(sampled_b);

% Find the elements in sorted_sampled_c that are greater than corresponding elements in sorted_b
invalid_indices = find(sorted_sampled_b > sorted_c, 1);

% If any such elements are found, replace them and all subsequent elements with corresponding sorted_b values
if ~isempty(invalid_indices)
    sorted_c(invalid_indices:end) = sorted_sampled_b(invalid_indices:end);
end

% Reorder sorted_sampled_c back to match the original sequence of c
sorted_c(original_indices) = sorted_c;

%% Find Size Distribution Statistics 
% Find an estimation for the equivelant radius distribution
equi_D = 2*nthroot(0.75*a.*b.*sampled_c,3); % for the ND grains
equi_D_TD = 2*nthroot(0.75*sampled_a.*sampled_b.*sorted_c,3); % for the TD grains

% Fit to a lognormal distribution
pd_grains        = fitdist(log(equi_D),'Normal');
mu_grain_size    = pd_grains.mu;
sigma_grain_size = pd_grains.sigma;

%% Binning

% find the bin edges
% Define cutoff_min and cutoff_max
cutoff_min = exp(mu_grain_size_Saltykov - min_sigma_cutoff * sigma_grain_size_Saltykov);
cutoff_max = exp(mu_grain_size_Saltykov + max_sigma_cutoff * sigma_grain_size_Saltykov);

% Define the number of bins
num_bins = round((cutoff_max - cutoff_min) / bin_size + 1);

% Calculate the bin size based on the new number of bins and the range
new_bin_size = (cutoff_max - cutoff_min) / (num_bins - 1);

% Generate the bin edges
bin_edges = linspace(cutoff_min, cutoff_max, num_bins + 1);

% Create a list of what bin the elements belong to using the ND grains
diameters_binned = discretize(equi_D,bin_edges);

% we will also need one for the TD grains for further calculations
diameters_binned_TD = discretize(equi_D_TD,bin_edges);

%% Find B/A & C/A Ratios and Neighbor Distribution

% Initalise
alpha_b_over_a       = zeros(1,num_bins);
beta_b_over_a        = zeros(1,num_bins);
alpha_c_over_a       = zeros(1,num_bins);
beta_c_over_a        = zeros(1,num_bins);
mu_neighbors         = zeros(1,num_bins);
sigma_neighbors      = zeros(1,num_bins);

% Loop over each bin to find the statistics
for bin = 1:num_bins
    % B/A & C/A Ratios
    b_over_a = b(diameters_binned==bin)./a(diameters_binned==bin);
    c_over_a = sampled_c(diameters_binned==bin)./a(diameters_binned==bin);
    
    pd_b_over_a         = fitdist(b_over_a,'Beta');
    alpha_b_over_a(bin) = pd_b_over_a.a;
    beta_b_over_a(bin)  = pd_b_over_a.b;
    pd_c_over_a         = fitdist(c_over_a,'Beta');
    alpha_c_over_a(bin) = pd_c_over_a.a;
    beta_c_over_a(bin)  = pd_c_over_a.b;

    % Neighbor distributions
    % Here as a=b>>c we can assume that in a sphere of one diameter the
    % only neighbours that need to be considered are the ones that would be
    % within one radius of the ND map in 2D.
    C = C_full_TD(diameters_binned_TD==bin,:); 
    D_binned  = equi_D_TD(diameters_binned_TD==bin);
    
    distances = pdist2(C_full_TD(:,1:2), C(:,1:2)); % calculate distances between all pairs of points
    neighbors = sum(distances < D_binned', 1)'; % count neighbors for each point
    
    % Fit to a lognormal distribution
    pd_neighbors = fitdist(neighbors,'Lognormal');
    mu_neighbors(bin)    = pd_neighbors.mu;
    sigma_neighbors(bin) = pd_neighbors.sigma;

end

%% Axis ODF
% this can be assumed from the short axis being orientated along the RD
% direction
Axis_ODF = [pi/2,0,0,1,50000];

%% ODF & MDF
% Find the ODF
odf = calcDensity(grains_ND.meanOrientation);

% sample 5000 points from this
sampled_odf = odf.discreteSample(5000);
sampled_odf = [sampled_odf.phi1,sampled_odf.Phi,sampled_odf.phi2,...
    ones(5000,1),ones(5000,1)];

% CalcOrientations - better as input to DREAM3D

% Find the MDF - can't currently be read into dream3d
%mdf = calcMDF(odf);

%% Write to JSON file to be read by python 

% Create a structure to hold the variables
data = struct();

% binning
data.num_bins = num_bins;
data.bin_number = bin_edges(1:end-1);
data.feature_diameter_info = [bin_size, cutoff_max, cutoff_min];

% grain sizes
data.mu_grain_size = mu_grain_size;
data.sigma_grain_size = sigma_grain_size;

% Axis ratios
data.alpha_b_over_a = alpha_b_over_a;
data.beta_b_over_a = beta_b_over_a;
data.alpha_c_over_a = alpha_c_over_a;
data.beta_c_over_a = beta_c_over_a;

% Neighbor distributions
data.mu_neighbors = mu_neighbors;
data.sigma_neighbors = sigma_neighbors;

% Axis ODFs
data.Axis_ODF = Axis_ODF;

% Sampled ODFs
data.sampled_odf = sampled_odf;

% Save the structure to a JSON file
json_data = jsonencode(data);
fid = fopen('ebsd_parameters.json', 'w');
fprintf(fid, json_data);
fclose(fid);

%% Functions

function grains = load_ebsd_data(fname,CS)

% Import and create an EBSD variable containing the data
ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertSpatial2EulerReferenceFrame');

% Al phase selection
ebsd_Al = ebsd('Aluminium');

% Al grains
[grains,ebsd_Al.grainId,ebsd_Al.mis2mean] = calcGrains(ebsd_Al,'angle',10*degree);
ebsd_Al(grains(grains.grainSize<10))      = []; % clean up the grains smaller than 10 pixels
[grains,ebsd_Al.grainId,ebsd_Al.mis2mean] = calcGrains(ebsd_Al,'angle',10*degree);

F = meanFilter; % filling grains
F.numNeighbours = 4;
grains = smooth(grains);
end