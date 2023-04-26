%% Requirments
% The Statistics and Machine Learning Toolbox
% Automatic Saltykov - https://uk.mathworks.com/matlabcentral/fileexchange/73726-automatic-saltykov

%% Data Paths
close all
clear 
clc

% Add the path to the data
fname = 'test_data\test_data.ctf';

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
pd_grains        = fitdist(log(2*grains.equivalentRadius),'Normal');
mu_grain_size    = pd_grains.mu;
sigma_grain_size = pd_grains.sigma;

% Find the Saltykov distribution
[freq, centers]           = autoSaltykov(log(2*grains.equivalentRadius));
eq_di_Satlykov            = repelem(exp(centers),round(freq/min(freq(freq>0))));
pd_grains_Saltykov        = fitdist(log(eq_di_Satlykov'),'Normal');
mu_grain_size_Saltykov    = pd_grains_Saltykov.mu;
sigma_grain_size_Saltykov = pd_grains_Saltykov.sigma;

%% Binning
%%% Find the parameters of the bins that will be used in DREAM3D
cutoff_min = exp(mu_grain_size_Saltykov-min_sigma_cutoff*sigma_grain_size_Saltykov);
cutoff_max = exp(mu_grain_size_Saltykov+max_sigma_cutoff*sigma_grain_size_Saltykov);
num_bins   = round((cutoff_max-cutoff_min)/bin_size + 1);
bin_edges  = linspace(cutoff_min,cutoff_max,num_bins+1);

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
    
    % Fit to a beta distribution and extract stats
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
    
    % Compute radii and number of neighbors in 3D
    r = sqrt((0.785*(D_binned/2).^2./sum(distances'<D_binned,2))); 
    neighbors_3D = 0.74*(D_binned./(2*r)).^3; 

    pd_neighbors = fitdist(neighbors_3D,'Lognormal');
    mu_neighbors(bin)    = pd_neighbors.mu;
    sigma_neighbors(bin) = pd_neighbors.sigma;
end

%% Axis ODF
% I am currently unable to figure out how to do this. This will have to
% come from an assumption that the 2D and 3D distributions have the
% same form. This is based on a distribution with n boxes
n = 12;
hist_outcome     = histogram(omega*(180/pi),n).Values;
probabilities    = hist_outcome/sum(hist_outcome);
probabilities    = repmat(probabilities,1,2);
probabilities_3D = repelem(probabilities,(2*n)^2).*repmat(repelem(probabilities,2*n),1,2*n).*repmat(probabilities,1,(2*n)^2);
eulers           = 0:360/(2*n-1):360;
Axis_ODF         = [repelem(eulers,(2*n)^2)',repmat(repelem(eulers,(2*n)),1,(2*n))',repmat(eulers,1,(2*n)^2)',(probabilities_3D/min(probabilities_3D))',(ones(1,(2*n)^3))'];

%% ODF & MDF
% Find the ODF
odf = calcDensity(grains.meanOrientation);

% Find the MDF
mdf = calcMDF(odf);
