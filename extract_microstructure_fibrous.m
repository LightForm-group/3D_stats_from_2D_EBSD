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
fname       = 'test_data\NDRD.ctf'; % Path to the NDRD EBSD map
grains_fib  = load_ebsd_data(fname,CS);

% RDTD - equiaxed structure
fname       = 'test_data\test_data.ctf'; % Path to the RDTD EBSD map
grains_equi = load_ebsd_data(fname,CS);

%% Process 
% For elongated grains we can assume that they are roughly elongated boxes
% with two EBSD maps taken, one orthogonal to the long axis and one with
% the long axis in the plane

% Find the axis and centroids of the NDRD structure
[~,a,~] = grains_fib.fitEllipse; % coincides with the actual grain area
[x,y]           = centroid(grains_fib); % find the grain centres
C_full_fib      = cat(2,x,y); % Combine

% Find the axis and centroids of the RDTD structure
[~,b,c] = grains_equi.fitEllipse; % coincides with the actual grain area
[x,y]             = centroid(grains_equi); % find the grain centres
C_full_equi       = cat(2,x,y); % Combine

% Find the distribution of the equivelant radius
% Here the volume is found from taking the axis of the RDTD and randomly 
% selecting a long axis from the NDRD distribution while keeping the 
% distribution consistent. this is subsequently used for further 
% calculations.
sampled_a = exp(normrnd(mean(log(a)), std(log(a)), [1 size(b,1)]))';

%% Find Size Distribution Statistics 
% Find an estimation for the equivelant radius distribution
equi_D = 2*nthroot(sampled_a.*b.*c,3);

% Fit to a lognormal distribution
pd_grains        = fitdist(log(equi_D),'Normal');
mu_grain_size    = pd_grains.mu;
sigma_grain_size = pd_grains.sigma;

%% Binning

% find the bin edges
cutoff_min = exp(mu_grain_size-min_sigma_cutoff*sigma_grain_size);
cutoff_max = exp(mu_grain_size+max_sigma_cutoff*sigma_grain_size);
num_bins   = round((cutoff_max-cutoff_min)/bin_size + 1);
bin_edges  = linspace(cutoff_min,cutoff_max,num_bins+1);

% Create a list of what bin the elements belong to
diameters_binned = discretize(equi_D,bin_edges);

%% Find B/A & C/A Ratios and Neighbor Distribution

% Initalise
alpha_b_over_a       = zeros(1,num_bins);
beta_b_over_a        = zeros(1,num_bins);
alpha_c_over_a       = zeros(1,num_bins);
beta_c_over_a        = zeros(1,num_bins);
mu_neighbors(bin)    = zeros(1,num_bins);
sigma_neighbors(bin) = zeros(1,num_bins);

% Loop over each bin to find the statistics
for bin = 1:num_bins
    % B/A & C/A Ratios
    b_over_a = b(diameters_binned==bin)./sampled_a(diameters_binned==bin);
    c_over_a = c(diameters_binned==bin)./sampled_a(diameters_binned==bin);
    
    pd_b_over_a         = fitdist(b_over_a,'Beta');
    alpha_b_over_a(bin) = pd_b_over_a.a;
    beta_b_over_a(bin)  = pd_b_over_a.b;
    pd_c_over_a         = fitdist(c_over_a,'Beta');
    alpha_c_over_a(bin) = pd_c_over_a.a;
    beta_c_over_a(bin)  = pd_c_over_a.b;

    % Neighbor distributions
    % Here as a>>>b>c we can assume that in a sphere of one diameter the 
    % structure is approximatly an extruded 2D structure and we can use 
    % the 2D neighbor distribution from the RDTD map
    C = C_full_equi(diameters_binned==bin,:); 
    D_binned  = equi_D(diameters_binned==bin);
    
    distances = pdist2(C_full(:,1:2), C(:,1:2)); % calculate distances between all pairs of points
    neighbors = sum(distances < D_binned, 1)'; % count neighbors for each point
    
    % Fit to a lognormal distribution
    pd_neighbors = fitdist(neighbors,'Lognormal');
    mu_neighbors(bin)    = pd_neighbors.mu;
    sigma_neighbors(bin) = pd_neighbors.sigma;

end

%% Axis ODF
% this can be assumed from the long axis being orientated along the RD
% direction
Axis_ODF = [90,0,0,50000,1];

%% ODF & MDF
% Find the ODF
odf = calcDensity(grains.meanOrientation);

% sample 5000 points from this
sampled_odf = odf.discreteSample(5000);
sampled_odf = [sampled_odf.phi1,sampled_odf.Phi,sampled_odf.phi2,...
    ones(5000,1),ones(5000,1)];

% CalcOrientations - better as input to DREAM3D

% Find the MDF - can't currently be read into dream3d
%mdf = calcMDF(odf);

%% Write to txt files to be read by python 

% binning
writematrix(num_bins, 'temp/num_bins.txt')
writematrix(bin_edges(1:end-1), 'temp/bin_number.txt')
writematrix([bin_size,cutoff_max,cutoff_min], 'temp/Feature_Diameter_Info.txt')

% grain sizes
writematrix(mu_grain_size, 'temp/mu_grain_size.txt')
writematrix(sigma_grain_size, 'temp/sigma_grain_size.txt')

% Axis ratios
writematrix(alpha_b_over_a, 'temp/alpha_b_over_a.txt')
writematrix(beta_b_over_a, 'temp/beta_b_over_a.txt')
writematrix(alpha_c_over_a, 'temp/alpha_c_over_a.txt')
writematrix(beta_c_over_a, 'temp/beta_c_over_a.txt')

% Neighbor distributions
writematrix(mu_neighbors, 'temp/mu_neighbors.txt')
writematrix(sigma_neighbors, 'temp/sigma_neighbors.txt')

% Axis ODFs
writematrix(Axis_ODF, 'temp/Axis_ODF.txt')

% Sampled ODFs
writematrix(sampled_odf, 'temp/sampled_odf.txt')

%% Run python script
pyrunfile('write_DREAM3D_json.py')

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