
%%% CODE FOR RECREATING CMT PAPER FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%%% This m file contains the script to recreate Figure 4b in the paper.
%%% The script reads reconstructed data from the ../data_output/ folder and
%%% puts the tif image in the  ../Figures folder.
%%% FIGURE 4b:  Subject 2 Sagittal Fat/Water Image

% clean slate
clear all; close all; clc;

code_path = fileparts(mfilename('fullpath'));
data_path = sprintf('%s/../data_output', code_path);

mat_file{1} = 'Sub2_DS_outParamsQPBO_PYTHON';

%%%%% Cropping Range %%%%%%%

P = 320;
crop_range_LR = P/2-(1.5*P/4):P/2+(1.5*P/4);
crop_range_AP = P/2-(P/4):P/2+(P/4);

%%%%%%% Sagittal Slice %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n =  110;
x = 1:161;
y = 1:720;

im = zeros(size(x,2),size(y,2),2);

load( sprintf('%s/%s.mat', data_path, mat_file{1}) );
outParams.species(1).amps = outParams.species(1).amps(crop_range_LR,crop_range_AP,:);
outParams.species(2).amps = outParams.species(2).amps(crop_range_LR,crop_range_AP,:);
    
im(:,:,1) = squeeze(abs(outParams.species(1).amps(n,x,y)));
im(:,:,2) = squeeze(abs(outParams.species(2).amps(n,x,y)));

A = flip(permute(im,[2 1 3]),2);

A_Final = [ padarray(A(:,:,1),[5, 8])   padarray(A(:,:,2),[5, 8])  ];  
figure; imagesc(A_Final); axis image off; colormap(gray);

F = getframe;
outfile = sprintf('%s/../figures/Figure_4b.tif', code_path);
imwrite(F.cdata,outfile,'tif');


