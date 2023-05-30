%% Requirments
close all
clear all
clc

%% Load in Euler angles text file
Eulers = 180/pi*importdata('EulerAngles.txt');
matrix_of_matrices = num2cell(Eulers, 2);

%% Reshape into 3D structure for sectioning
% Matrix dimensions
rve_dims = [256,256,256];

% Reshape to 3D structure
reshaped_3D_cell = reshape(matrix_of_matrices, [rve_dims(1), rve_dims(2), rve_dims(3)]);

%% Take slices in all 3 planes 
% Take a numberof slices from the strctur and build synthetic EBSD data
% from these
for ind = linspace(16, 256, 16)
extract_plane(reshaped_3D_cell,'xy',ind)
extract_plane(reshaped_3D_cell,'xz',ind)
extract_plane(reshaped_3D_cell,'zy',ind)
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
grains_synthetic_combined_xz = combine_ebsd_slices(CS,'xz');
grains_synthetic_combined_zy = combine_ebsd_slices(CS,'zy');

%% Load experimental data into Mtez

fname = 'test_data\test_data.ctf';

ebsd_exp = EBSD.load(fname,CS,'interface','ctf',...
  'convertSpatial2EulerReferenceFrame');

% Al phase selection
ebsd_Al_exp = ebsd_exp('Aluminium');

% Al grains
[grains,ebsd_Al_exp.grainId,ebsd_Al_exp.mis2mean] = calcGrains(ebsd_Al_exp,'angle',10*degree);
ebsd_Al_exp(grains(grains.grainSize<10))      = []; % clean up the grains smaller than 10 pixels
[grains,ebsd_Al_exp.grainId,ebsd_Al_exp.mis2mean] = calcGrains(ebsd_Al_exp,'angle',10*degree);


%% Extract synthetic data
% Radii
radii_synthetic = [cellfun(@(x) x.equivalentRadius, grains_synthetic_combined_xy(1:end), 'UniformOutput', false);
         cellfun(@(x) x.equivalentRadius, grains_synthetic_combined_xz(1:end), 'UniformOutput', false);
         cellfun(@(x) x.equivalentRadius, grains_synthetic_combined_zy(1:end), 'UniformOutput', false)];
radii_synthetic = cat(1, radii_synthetic{:});

% Neighbor number
neighbors_synthetic = [cellfun(@(x) numNeighbors(x), grains_synthetic_combined_xy(1:end), 'UniformOutput', false);
         cellfun(@(x) numNeighbors(x), grains_synthetic_combined_xz(1:end), 'UniformOutput', false);
         cellfun(@(x) numNeighbors(x), grains_synthetic_combined_zy(1:end), 'UniformOutput', false)];
neighbors_synthetic = cat(1, neighbors_synthetic{:});

% Ellipse parameters
[omega_synth_xy, a_synth_xy, b_synth_xy] = cellfun(@(x) x.fitEllipse, grains_synthetic_combined_xy(1:end), 'UniformOutput', false);
[omega_synth_xz, a_synth_xz, b_synth_xz] = cellfun(@(x) x.fitEllipse, grains_synthetic_combined_xz(1:end), 'UniformOutput', false);
[omega_synth_zy, a_synth_zy, b_synth_zy] = cellfun(@(x) x.fitEllipse, grains_synthetic_combined_zy(1:end), 'UniformOutput', false);

omega_synth = [cat(1, omega_synth_xy{:});cat(1, omega_synth_xz{:});cat(1, omega_synth_zy{:})];
a_synth     = [cat(1, a_synth_xy{:});cat(1, a_synth_xz{:});cat(1, a_synth_zy{:})];
b_synth     = [cat(1, b_synth_xy{:});cat(1, b_synth_xz{:});cat(1, b_synth_zy{:})];

%% Compare grain shape parameters

% Sizes
fig1 = figure;
hold on
hist_exp   = histogram(grains.equivalentRadius,Normalization="probability");
hist_synth = histogram(radii_synthetic, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('Radius / \mum');
ylim([0,0.25])
xlim([0,80])
hold off
saveas(fig1,'comparisons/radii.png')

% experimental ellipse parameters
[omega_exp,a_exp,b_exp] = grains.fitEllipse; % coincides with the actual grain area

% B/A ratio
fig2 = figure;
hold on
hist_exp   = histogram(b_exp./a_exp,Normalization="probability");
hist_synth = histogram(b_synth./a_synth, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('B/A Axis Ratios');
ylim([0,0.25])
%xlim([0,100])
hold off
saveas(fig2,'comparisons/b_over_a.png')

% Omega
fig3 = figure;
hold on
hist_exp   = histogram(180/pi*omega_exp,12,Normalization="probability");
hist_synth = histogram(180/pi*omega_synth, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('\omega / degrees');
ylim([0,0.25])
xlim([0,180])
hold off
saveas(fig3,'comparisons/omega.png')

% Neighbors
fig4 = figure;
hold on
hist_exp   = histogram(numNeighbors(grains),Normalization="probability");
hist_synth = histogram(neighbors_synthetic, hist_exp.BinEdges,Normalization="probability");
legend({'Experimental','Synthetic'},'Location','northeast')
legend boxoff 
ylabel('Probability');
xlabel('Number of Neighbours');
ylim([0,0.25])
xlim([0,25])
hold off
saveas(fig4,'comparisons/neighbors.png')

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


