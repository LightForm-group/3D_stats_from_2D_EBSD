%% Requirments
close all
clear all
clc

%% Specify Crystal, Grain Structure and Specimen Symmetries

% crystal symmetry
CS = {... 
  'notIndexed',...
  crystalSymmetry('m-3m', [4 4 4], 'mineral', 'Aluminium', 'color', [0.53 0.81 0.98]),...
  'notIndexed'};

% plotting convention
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','outOfPlane');

fname = 'test_data\NDRD.ctf';

ebsd = EBSD.load(fname,CS,'interface','ctf',...
  'convertSpatial2EulerReferenceFrame');
%ebsd = loadEBSD_crc(fname);

% Al phase selection
ebsd_Al = ebsd('Aluminium');

% Al grains
[grains,ebsd_Al.grainId,ebsd_Al.mis2mean] = calcGrains(ebsd_Al,'angle',10*degree);
ebsd_Al(grains(grains.grainSize<10))      = []; % clean up the grains smaller than 10 pixels
[grains,ebsd_Al.grainId,ebsd_Al.mis2mean] = calcGrains(ebsd_Al,'angle',10*degree);

% Plot
ipfKey = ipfColorKey(ebsd);
colors = ipfKey.orientation2color(ebsd.orientations);
plot(ebsd,colors)


