
%%% CODE FOR RECREATING CMT PAPER FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%%% This m file contains the script to recreate Figure 2f in the paper.
%%% The script reads reconstructed data from the ../data_output/ folder and
%%% puts the tif image in the  ../Figures folder.
%%% FIGURE 2f :  No Shim and Dynamic Shim sagittal fieldmaps for Subject 2

% clean slate
clear all; close all; clc;

code_path = fileparts(mfilename('fullpath'));
data_path = sprintf('%s/../data_output', code_path);

n = zeros(4);

%%%% 256 Recons : Dynamically Shimmed Python Recons 4 echoes

mat_file{1} = 'Sub2_NS_outParamsQPBO_PYTHON';
mat_file{2} = 'Sub2_DS_outParamsQPBO_PYTHON';


%%%%% Cropping Range %%%%%%%

P = 320;
crop_range_LR = P/2-(1.5*P/4):P/2+(1.5*P/4);
crop_range_AP = P/2-(P/4):P/2+(P/4);


%%%% Sagittal
n =  40:140;
x = 107;
y = 1:720;
im = zeros(size(n,2),size(y,2),2);
mask = zeros(size(n,2),size(y,2),2);

for i = 1 : 2
    i
    load( sprintf('%s/%s.mat', data_path, mat_file{i}) );
    outParams.fieldmap = outParams.fieldmap(crop_range_LR,crop_range_AP,:);
    clearvars -except outParams i n data_path mat_file im x y mask code_path crop_range_LR crop_range_AP
    
    im(:,:,i) = squeeze(-outParams.fieldmap(x,n,y));
end

A = flip(permute(im,[2 1 3]),2);

A_Final = [ padarray(A(:,:,1),[5, 8])   padarray(A(:,:,2),[5, 8])  ]; 
figure; imagesc(-A_Final,[ -250 250]); axis image off; colormap(jet);
F = getframe;

outfile = sprintf('%s/../figures/Figure_2f.tif', code_path);
imwrite(F.cdata,outfile,'tif');


