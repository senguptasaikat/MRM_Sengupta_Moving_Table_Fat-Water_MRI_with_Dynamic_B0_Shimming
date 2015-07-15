
%%% CODE FOR RECREATING CMT PAPER FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%%% This m file contains the script to recreate Figure 1 in the paper.
%%% The script reads reconstructed data from the ../data_output/ folder and
%%% puts the tif image in the  ../Figures folder.
%%% FIGURE 1 :  MULTIECHO CORONAL IMAGE

% clean slate
clear all; close all; clc;

code_path = fileparts(mfilename('fullpath'));
data_path = sprintf('%s/../data_output', code_path);

n = zeros(4);

%%%% 256 Recons : Dynamically Shimmed Python Recons 4 echoes

mat_file{1} = 'Sub1_DS_echo1_PYTHONREC';
mat_file{2} = 'Sub1_DS_echo2_PYTHONREC';
mat_file{3} = 'Sub1_DS_echo3_PYTHONREC';
mat_file{4} = 'Sub1_DS_echo4_PYTHONREC';

n =  160;
x = 1:720;
y = 55:265;

im = zeros(size(x,2),size(y,2),4);

for i = 1 : 4
    load( sprintf('%s/%s.mat', data_path, mat_file{i}) );    
    im(:,:,i) = img(x,y,n);
    clear img
    disp(sprintf('mat file %d loaded',i));  
end


im = flipdim(im,2);

A = abs(im);
A = A-min(A(:));
A = A/max(A(:));
A = reshape(A, [size(A,1) size(A,2) 1 size(A,3)]);

figure; montage(A); colormap gray; axis image off
F = getframe;

outfile = sprintf('%s/../figures/Figure_1.tif', code_path);
imwrite(F.cdata,outfile,'tif');


