%% Requirments
% The Statistics and Machine Learning Toolbox
% Automatic Saltykov - https://uk.mathworks.com/matlabcentral/fileexchange/73726-automatic-saltykov
% TODO: Add as seperate folder
%% Data Paths
close all
clear 
clc

% Add the path to the data
fname = 'test_data\test_data.ctf';
output_name = 'test_structure';

%% Specify Dream3D statistics generation inputs
% For the size distribution - used for binning as described here
% (http://127.0.0.1:32456/Filters/SyntheticBuildingFilters/StatsGeneratorFilter/)
num_bins = 6; 
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

%% Find Size Distribution Statistics 
% For equiaxed grains we must account for taking 2D slices of a 3D
% structure which roughly represents the corpuscule problem. The Saltykov
% algorithm can be used to account for this and adjust the size distribution. 

% Fit to a normal distribution
% Calculated to compare to the adjusted distribution
pd_grains        = fitdist(log(2*grains.equivalentRadius),'Normal'); % fit to a normal distribution
mu_grain_size    = pd_grains.mu; % Mean log grain size
sigma_grain_size = pd_grains.sigma; % Standard deviation in log grain size

% Find the Saltykov distribution
[freq, centers]           = autoSaltykov(log(2*grains.equivalentRadius)); % Run Saltykov algorithm
eq_di_Satlykov            = repelem(exp(centers),round(freq/min(freq(freq>0)))); % Extract diameters
pd_grains_Saltykov        = fitdist(log(eq_di_Satlykov'),'Normal'); % fit to a normal distribution
mu_grain_size_Saltykov    = pd_grains_Saltykov.mu; 
sigma_grain_size_Saltykov = pd_grains_Saltykov.sigma;

%% Binning

% Define cutoff_min and cutoff_max
cutoff_min = exp(mu_grain_size_Saltykov - min_sigma_cutoff * sigma_grain_size_Saltykov);
cutoff_max = exp(mu_grain_size_Saltykov + max_sigma_cutoff * sigma_grain_size_Saltykov);

% Calculate the new bin size based on the predefined number of bins and the range
bin_size = (cutoff_max - cutoff_min) / (num_bins - 1);

% Generate the bin edges
bin_edges = linspace(cutoff_min, cutoff_max, num_bins + 1);


% Create a list of what bin the elements belong to, this requires adjusting
% the origional diameters to fit the distribution calculated using the
% Saltykov algorithm above
D = 2*grains.equivalentRadius;

% Compute z-scores of each value in log(D)
z_scores_2D = (log(D) - mu_grain_size) / sigma_grain_size;

% Transform z-scores to those of the Saltykov distribution
log_D_saltykov = z_scores_2D * sigma_grain_size_Saltykov + mu_grain_size_Saltykov;

% Convert back to non log values
D_saltykov = exp(log_D_saltykov);

% Find the bin labels for all the grains
diameters_binned = discretize(D_saltykov,bin_edges);

%% Find B/A & C/A Ratios and Neighbor Distribution
% Fit an ellipse to the grain area and find the ellipsis axis lengths,
% angle to horizontal and centre coordinates. 
[omega,a,b] = grains.fitEllipse; % coincides with the actual grain area
[x,y]  = centroid(grains); % find the grain centres
C_full = cat(2,x,y); % Combine

% Initalise vectors
alpha_b_over_a  = zeros(1,num_bins);
beta_b_over_a   = zeros(1,num_bins);
alpha_c_over_a  = zeros(1,num_bins);
beta_c_over_a   = zeros(1,num_bins);
mu_neighbors    = zeros(1,num_bins);
sigma_neighbors = zeros(1,num_bins);

% Loop over each bin to find the statistics
for bin = 1:num_bins
    % Here assume that B=C in the case of only one slice of data
    b_over_a = b(diameters_binned==bin)./a(diameters_binned==bin);
    
    % Fit to a beta distribution and extract alpha and beta values
    pd_b_over_a         = fitdist(b_over_a,'Beta');
    alpha_b_over_a(bin) = pd_b_over_a.a;
    beta_b_over_a(bin)  = pd_b_over_a.b;
    alpha_c_over_a(bin) = pd_b_over_a.a;
    beta_c_over_a(bin)  = pd_b_over_a.b;

    % Find the diameters and coordinates of the grains in the bin
    C = C_full(diameters_binned==bin,:); 
    D_binned  = D(diameters_binned==bin);
    
    % Compute distances between all pairs of points
    distances = sqrt((C_full(:,1)-C(:,1)').^2 + (C_full(:,2)-C(:,2)').^2);
    
    % Compute the number of neighbors in 2D and then assume they are
    % perfectly packed as circles of radius r in a larger circle of radius
    % D_binned.
    r = sqrt((0.785*(D_binned/2).^2./sum(distances'<D_binned,2))); 

    % Compute estimate of neighbors in 3D assuming a perfect packing of 
    % spheres with radius r in a sphere of diameter D_binned
    neighbors_3D = 0.74*(D_binned./(2*r)).^3; 
    
    % Fit to a lognormal distribution
    pd_neighbors = fitdist(neighbors_3D,'Lognormal');
    mu_neighbors(bin)    = pd_neighbors.mu;
    sigma_neighbors(bin) = pd_neighbors.sigma;
end

%% Axis ODF
% I am currently unable to figure out how to do this. This will have to
% come from an assumption that the 2D and 3D distributions have the
% same form. This is based on a distribution with n boxes

% The number of segments to split each 180 degree section into
n = 12; 

% Find the relative probabilities of th angles
hist_outcome     = histogram(omega*(180/pi),n).Values;
probabilities    = hist_outcome/sum(hist_outcome);
probabilities    = repmat(probabilities,1,2); % expand to 360 degrees

% Build this into 3D with the weighting used in DREAM3D given by the
% relative probability
probabilities_3D = repelem(probabilities,(2*n)^2).*repmat(repelem(probabilities,2*n),1,2*n).*repmat(probabilities,1,(2*n)^2);

% Build a 5x(2n)^3 dataset with the euler angles, relative weight and a
% sigma of 1
eulers           = 0:2*pi/(2*n-1):2*pi;
Axis_ODF         = [repelem(eulers,(2*n)^2)',repmat(repelem(eulers,(2*n)),1,(2*n))',repmat(eulers,1,(2*n)^2)',(probabilities_3D/min(probabilities_3D))',(ones(1,(2*n)^3))'];

%% ODF & MDF
% Find the ODF
odf = calcDensity(grains.meanOrientation);

% sample 5000 points from this
sampled_odf = odf.discreteSample(5000);
sampled_odf = [sampled_odf.phi1,sampled_odf.Phi,sampled_odf.phi2,...
    ones(5000,1),ones(5000,1)];

% Find the MDF - can't currently be read into dream3d
%mdf = calcMDF(odf);

%% Write to JSON file to be read by python 

% Create a structure to hold the variables
data = struct();

% binning
data.num_bins = num_bins;
data.bin_number = bin_edges(1:end-1);
data.feature_diameter_info = [bin_size, max_sigma_cutoff, min_sigma_cutoff];

% grain sizes
data.mu_grain_size = mu_grain_size_Saltykov;
data.sigma_grain_size = sigma_grain_size_Saltykov;

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



















