%% Requirments
close all
clear all
clc

%% Load in Euler angles text file
Eulers = 180/pi*importdata('EulerAngles_rolled.txt');
matrix_of_matrices = num2cell(Eulers, 2);

%% Reshape into 3D structure for sectioning
% Matrix dimensions
rve_dims = [256,256,256];

% Reshape to 3D structure
reshaped_3D_cell = reshape(matrix_of_matrices, [rve_dims(1), rve_dims(2), rve_dims(3)]);

%% Take slices of the planes under consideration
% Take a numberof slices from the structure and build synthetic EBSD data
% from these
for ind = linspace(16, 256, 16)
extract_plane(reshaped_3D_cell,'xy',ind) % Equiaxed
extract_plane(reshaped_3D_cell,'zy',ind) % Elongated
end

% Free up memory
clear

%% Load synthetic data into Mtex and visualise

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98]),...
  'notIndexed'};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

% Extract sybnthetic ebsd data from text files and load them into MTex
grains_synthetic_combined_xy = combine_ebsd_slices(CS,'xy');
grains_synthetic_combined_zy = combine_ebsd_slices(CS,'zy');

%% Load experimental data into Mtex

%%% ND Data
fname_ND = 'test_data\Fibrous model\ND\ND.ctf';

ebsd_exp_ND = EBSD.load(fname_ND,CS,'interface','ctf',...
  'convertSpatial2EulerReferenceFrame');

% Al phase selection
ebsd_Al_exp_ND = ebsd_exp_ND('Aluminium');

% Al grains
[grains_ND,ebsd_Al_exp_ND.grainId,ebsd_Al_exp_ND.mis2mean] = calcGrains(ebsd_Al_exp_ND,'angle',10*degree);
ebsd_Al_exp_ND(grains_ND(grains_ND.grainSize<10))      = []; % clean up the grains smaller than 10 pixels
[grains_ND,ebsd_Al_exp_ND.grainId,ebsd_Al_exp_ND.mis2mean] = calcGrains(ebsd_Al_exp_ND,'angle',10*degree);

%%% TD Data
fname_TD = 'test_data\Fibrous model\TD\TD.ctf';

ebsd_exp_TD = EBSD.load(fname_TD,CS,'interface','ctf',...
  'convertSpatial2EulerReferenceFrame');

% Al phase selection
ebsd_Al_exp_TD = ebsd_exp_TD('Aluminium');

% Al grains
[grains_TD,ebsd_Al_exp_TD.grainId,ebsd_Al_exp_TD.mis2mean] = calcGrains(ebsd_Al_exp_TD,'angle',10*degree);
ebsd_Al_exp_TD(grains_TD(grains_TD.grainSize<10))      = []; % clean up the grains smaller than 10 pixels
[grains_TD,ebsd_Al_exp_TD.grainId,ebsd_Al_exp_TD.mis2mean] = calcGrains(ebsd_Al_exp_TD,'angle',10*degree);

%% Extract synthetic data
%%% ND
% Radii
radii_synthetic_ND = cellfun(@(x) x.equivalentRadius, grains_synthetic_combined_xy(1:end), 'UniformOutput', false);
radii_synthetic_ND = cat(1, radii_synthetic_ND{:});

% Neighbor number
neighbors_synthetic_ND = cellfun(@(x) numNeighbors(x), grains_synthetic_combined_xy(1:end), 'UniformOutput', false);
neighbors_synthetic_ND = cat(1, neighbors_synthetic_ND{:});

% Ellipse parameters
[omega_synth_ND, a_synth_ND, b_synth_ND] = cellfun(@(x) x.fitEllipse, grains_synthetic_combined_xy(1:end), 'UniformOutput', false);
omega_synth_ND = cat(1, omega_synth_ND{:});
a_synth_ND = cat(1, a_synth_ND{:});
b_synth_ND = cat(1, b_synth_ND{:});

%%% TD
% Radii
radii_synthetic_TD = cellfun(@(x) x.equivalentRadius, grains_synthetic_combined_zy(1:end), 'UniformOutput', false);
radii_synthetic_TD = cat(1, radii_synthetic_TD{:});

% Neighbor number
neighbors_synthetic_TD = cellfun(@(x) numNeighbors(x), grains_synthetic_combined_zy(1:end), 'UniformOutput', false);
neighbors_synthetic_TD = cat(1, neighbors_synthetic_TD{:});

% Ellipse parameters
[omega_synth_TD, a_synth_TD, b_synth_TD] = cellfun(@(x) x.fitEllipse, grains_synthetic_combined_zy(1:end), 'UniformOutput', false);
omega_synth_TD = cat(1, omega_synth_TD{:});
a_synth_TD = cat(1, a_synth_TD{:});
b_synth_TD = cat(1, b_synth_TD{:});

%% Compare grain shape parameters

%%% ND
% Sizes
fig1 = figure;
hold on
hist_exp   = histogram(grains_ND.equivalentRadius,Normalization="probability");
hist_synth = histogram(radii_synthetic_ND, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('Radius / \mum');
ylim([0,0.25])
xlim([0,80])
hold off
saveas(fig1,'comparisons/radii_rolled_ND.png')

% experimental ellipse parameters
[omega_exp_ND,a_exp_ND,b_exp_ND] = grains_ND.fitEllipse; % coincides with the actual grain area

% B/A ratio
fig2 = figure;
hold on
hist_exp   = histogram(b_exp_ND./a_exp_ND,Normalization="probability");
hist_synth = histogram(b_synth_ND./a_synth_ND, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('B/A Axis Ratios');
ylim([0,0.25])
xlim([0,1])
hold off
saveas(fig2,'comparisons/b_over_a_rolled_ND.png')

% Neighbors
fig4 = figure;
hold on
hist_exp   = histogram(numNeighbors(grains_ND),Normalization="probability");
hist_synth = histogram(neighbors_synthetic_ND, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('Number of Neighbours');
ylim([0,0.25])
xlim([0,25])
hold off
saveas(fig4,'comparisons/neighbors_rolled_ND.png')

%%% TD
% Sizes
fig1 = figure;
hold on
hist_exp   = histogram(grains_TD.equivalentRadius,Normalization="probability");
hist_synth = histogram(radii_synthetic_TD, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('Radius / \mum');
ylim([0,0.25])
xlim([0,80])
hold off
saveas(fig1,'comparisons/radii_rolled_TD.png')

% experimental ellipse parameters
[omega_exp_TD,a_exp_TD,b_exp_TD] = grains_TD.fitEllipse; % coincides with the actual grain area

% B/A ratio
fig2 = figure;
hold on
hist_exp   = histogram(b_exp_TD./a_exp_TD,Normalization="probability");
hist_synth = histogram(b_synth_TD./a_synth_TD, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('C/A Axis Ratios');
ylim([0,0.25])
xlim([0,1])
hold off
saveas(fig2,'comparisons/b_over_a_rolled_TD.png')

% Neighbors
fig4 = figure;
hold on
hist_exp   = histogram(numNeighbors(grains_TD),Normalization="probability");
hist_synth = histogram(neighbors_synthetic_TD, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('Number of Neighbours');
ylim([0,0.25])
xlim([0,25])
hold off
saveas(fig4,'comparisons/neighbors_rolled_TD.png')

%% functions

function extract_plane(reshaped_3D_cell,slice_plane,slice_ind)

% Selct a plane
if strcmp( slice_plane,'xy' )
    slice = reshaped_3D_cell(:, :, slice_ind);  
elseif strcmp( slice_plane,'xz' )
    slice = reshaped_3D_cell(:, slice_ind, :);  
else
    slice = reshaped_3D_cell(slice_ind, :, :);  
end

% Reshape into NxN matrix
slice = squeeze(slice);

% Add coordinates and phase label to each row
phase = 1;
[rows, cols] = size(slice);
plane_with_extra = cell(rows, cols);
for i = 1:rows
    for j = 1:cols
        plane_with_extra{i, j} = [slice{i, j}, phase, j, i];
    end
end

% Write to a text file
filename = sprintf('temp/synthetic_ebsd_plane=%s_slice_ind=%g.txt', slice_plane, slice_ind);
fid = fopen(filename, 'w');
for i = 1:rows
    for j = 1:cols
        fprintf(fid, '%f %f %f %d %d %d\n', plane_with_extra{i, j});
    end
end
fclose(fid);
end


function grains_synthetic_combined = combine_ebsd_slices(CS,slice_plane)

grains_synthetic_combined = cell(1, size(linspace(16, 256, 16),2));
i = 1;
for ind = linspace(16, 256, 16)
% Load into Mtex
filename = sprintf('temp/synthetic_ebsd_plane=%s_slice_ind=%g.txt', slice_plane, ind);
ebsd_synthetic   = loadEBSD_generic(filename,'ColumnNames',{'Euler1','Euler2','Euler3','phase','x','y'},'cs',CS);

% Al phase selection
ebsd_Al_synthetic = ebsd_synthetic('Aluminium');

% Al grains
[grains_synthetic,ebsd_Al_synthetic.grainId,ebsd_Al_synthetic.mis2mean] = calcGrains(ebsd_Al_synthetic,'angle',10*degree);

grains_synthetic_combined{i} = grains_synthetic;
i = i+1;
end

end