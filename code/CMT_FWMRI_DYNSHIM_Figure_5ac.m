
%%% CODE FOR RECREATING CMT PAPER FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%%% This m file contains the script to recreate Figure 5 in the paper.
%%% The script reads reconstructed data from the ../data_output/ folder and
%%% puts the tif image in the  ../Figures folder.
%%% FIGURE 5 :  Dynamic Shim coronal FSF and R2Star Maps for subject 1

% clean slate
clear all; close all; clc;

code_path = fileparts(mfilename('fullpath'));
data_path = sprintf('%s/../data_output', code_path);

n = zeros(4);


%%%%% Cropping Range %%%%%%%

P = 320;
crop_range_LR = P/2-(1.5*P/4):P/2+(1.5*P/4);
crop_range_AP = P/2-(P/4):P/2+(P/4);


%%%% 256 Recons : Dynamically Shimmed Python Recons 4 echoes

mat_file{1} = 'Sub1_DS_outParamsQPBO_PYTHON';
load( sprintf('%s/%s.mat', data_path, mat_file{1}) );

outParams.species(1).amps = outParams.species(1).amps(crop_range_LR,crop_range_AP,:);
outParams.species(2).amps = outParams.species(2).amps(crop_range_LR,crop_range_AP,:);
outParams.r2starmap = outParams.r2starmap(crop_range_LR,crop_range_AP,:);


%%%%%%% Coronal
n =  78;
x = 30:220;
y = 1:720;

im = zeros(size(x,2),size(y,2),2);
mask = zeros(size(x,2),size(y,2),2);
    
R2s = squeeze(outParams.r2starmap(x,n,y));
water = squeeze(abs(outParams.species(1).amps(x,n,y)));
fat = squeeze(abs(outParams.species(2).amps(x,n,y)));

R2s_Final = flip(permute(R2s,[2 1 3]),2);
FSF = fat./(water+fat);
FSF_Final = flip(permute(FSF,[2 1 3]),2);

figure; imagesc(FSF_Final,[ 0 1]); axis image off; colormap(jet); colorbar;
F = getframe;
outfile = sprintf('%s/../figures/Figure_5a.tif', code_path);
imwrite(F.cdata,outfile,'tif');

figure; imagesc(R2s_Final,[ 0 120]); axis image off; colormap(jet); colorbar;
F = getframe;
outfile = sprintf('%s/../figures/Figure_5b.tif', code_path);
imwrite(F.cdata,outfile,'tif');
