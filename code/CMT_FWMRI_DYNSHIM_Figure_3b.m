
%%% CODE FOR RECREATING CMT PAPER FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%%% This m file contains the script to recreate Figure 3b in the paper.
%%% The script reads reconstructed data from the ../data_output/ folder and
%%% puts the tif image in the  ../Figures folder.
%%% FIGURE 3b :  Subject 1 Axial Fat/Water Image

% clean slate
clear all; close all; clc;

code_path = fileparts(mfilename('fullpath'));
data_path = sprintf('%s/../data_output', code_path);

mat_file{1} = 'Sub1_DS_outParamsQPBO_PYTHON';

%%%%% Cropping Range %%%%%%%

P = 320;
crop_range_LR = P/2-(1.5*P/4):P/2+(1.5*P/4);
crop_range_AP = P/2-(P/4):P/2+(P/4);

%%% Axial Slices %%%%%%%%%%%%%%%%%
n = [ 50, 130 , 200 , 290, 450,500]  ;
x = 45:200;
y = 20:160;
im = zeros(size(x,2),size(y,2),6);


load( sprintf('%s/%s.mat', data_path, mat_file{1}) );
outParams.species(1).amps = outParams.species(1).amps(crop_range_LR,crop_range_AP,:);
outParams.species(2).amps = outParams.species(2).amps(crop_range_LR,crop_range_AP,:);

    
water = squeeze(abs(outParams.species(1).amps(x,y,n)));
fat = squeeze(abs(outParams.species(2).amps(x,y,n)));

water = flip(flip(permute(water,[2 1 3]),1),2);
fat = flip(flip(permute(fat,[2 1 3]),1),2);

figure; imagesc([ water(:,:,1) fat(:,:,1); water(:,:,2) fat(:,:,2); water(:,:,3) fat(:,:,3); ...
    water(:,:,4) fat(:,:,4);water(:,:,5) fat(:,:,5); water(:,:,6) fat(:,:,6)]); axis image off; colormap(gray);

F = getframe;
outfile = sprintf('%s/../figures/Figure_3b.tif', code_path);
imwrite(F.cdata,outfile,'tif');


