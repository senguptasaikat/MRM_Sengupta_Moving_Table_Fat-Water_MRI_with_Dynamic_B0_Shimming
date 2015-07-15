
%%% CODE FOR RECREATING CMT PAPER FIGURES %%%%%%%%%%%%%%%%%%%%%%%%
%%% This m file contains the script to recreate Figure 5b in the paper.
%%% The script reads reconstructed data from the ../data_output/ folder and
%%% puts the tif image in the  ../Figures folder.
%%% FIGURE 5b :  Slicewise Fat/water plots for subject 1

% clean slate
clear all; close all; clc;

code_path = fileparts(mfilename('fullpath'));
data_path = sprintf('%s/../data_output', code_path);

mat_file{1} = 'Sub1_DS_outParamsQPBO_PYTHON';
load( sprintf('%s/%s.mat', data_path, mat_file{1}) );

    
water = squeeze(abs(outParams.species(1).amps)).* outParams.mask;
fat = squeeze(abs(outParams.species(2).amps)).* outParams.mask;

FSF = ones(size(water)) - water./(water+fat);
id = isnan(FSF); FSF(id) = 0;

water_max = max(water(:));
fat_max = max(fat(:));
    
water_fraction = zeros(size(water,3),1);
fat_fraction = zeros(size(water,3),1);
Total_signal = zeros(size(water,3),1);

for slice = 1:size(water,3) 
    FSF_slice = FSF(:,:,slice);
    water_slice = water(:,:,slice);
    fat_slice = fat(:,:,slice);
     
     id_w = find(water_slice );
    id_f = find(fat_slice );    
     
    water_fraction(slice)  = sum(water_slice(id_w));
    fat_fraction(slice)  = sum(fat_slice(id_f));
end
figure; plot(fat_fraction,'r'); hold on; plot(water_fraction,'b'); 
legend('fat','water');


F = getframe;
outfile = sprintf('%s/../figures/Figure_5b.tif', code_path);
imwrite(F.cdata,outfile,'tif');
