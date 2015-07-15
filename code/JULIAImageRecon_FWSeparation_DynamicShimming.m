
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % This file will read raw data from /data_input/ folder and call 
% % % % % the JULIA reconstruction script to reconstruct echo images 
% % % % %  The output images will be places in the /data_out/ folder.
% % % % % 
% % % % % Following that, tt will call FW separation algorithms
% % % % % FW separated data will be witten out
% % % % % in the /OutSlice folder at the OutParams structure.
% % % % % 
% % % % % In the Final step, the outParams structure will be read by the Dynamic
% % % % % Shimming tool to calculate shims. Dynamic Shims will be written 
%%%%%%%%% out per Phase encode or TR to the /data_out/ folder
 

tic;

fileno = 'Sub2_NS';

% uncomment individual line below to process other data sets
%fileno = 'Sub1_NS';
%fileno = 'Sub1_DS';
%fileno = 'Sub2_DS';

TR = num2str(6.6e-3);
npe_per_slice = num2str(128);
z_voxel = num2str(2.5);
nPy = 13732;
TE = [0.99 , 2.39, 3.79]* 1e-3;   % 3 echo 

pathname = '../data_output/';

 
%% KICK OFF RECON IN JULIA  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% recon.ji  'Echo mat file name'  TR  Voxel_Size
%%% No_of_profiles_per_slice

str1 = '/Applications/Julia-0.3.3.app/Contents/Resources/Julia/bin/julia';
str2 = ' -p 2 ';
str3 = 'Julia\ Recon\ Scripts/recon.jl ';


str4 = [fileno,'_echo1 ', TR,' ',z_voxel,' ',npe_per_slice];
command = [ str1 str2 str3 str4];
system(command) 


str4 = [fileno,'_echo2 ', TR,' ',z_voxel,' ',npe_per_slice];
command = [str1 str2 str3 str4];
system(command) 


str4 = [fileno,'_echo3 ', TR,' ',z_voxel,' ',npe_per_slice];
command = [str1 str2 str3 str4];
system(command) 


str4 = [fileno,'_echo4 ',TR,' ',z_voxel,' ',npe_per_slice]; 
command = [str1 str2 str3 str4];
system(command) 



%% PERFORMING FAT/WATER SEPARATION  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic
disp('Starting Fat/Water separation');

imDataParams.TE = TE;
imDataParams.FieldStrength = 3;
imDataParams.voxelSize  = [2.5 , 2.5 , str2double(z_voxel) ]; 
imDataParams.PrecessionIsClockwise = 1;

algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37];
algoParams.process_in_3D = true;


outParams = [];
outParams.species(1).name = 'water';
outParams.species(2).name = 'fat';

%%%% No need to open Echo 4 for MultiSeedRegionGrowing Algorithm
filename1  = [fileno,'_echo1', '_JULIA_RECON'];
load([pathname filename1]);
img1 = squeeze(Image); clear('Image');

   
filename2  = [fileno,'_echo2',  '_JULIA_RECON']; 
load([pathname filename2])  
img2 = squeeze(Image); clear('Image');
[ P, Q, slices] = size(img2);

filename3  = [fileno,'_echo3',  '_JULIA_RECON']; 
load([pathname filename3]);
img3 = squeeze(Image); clear('Image');

% filename4  = [fileno,'_echo4',  '_JULIA_RECON']; 
% load([pathname filename4]);
% img4 =squeeze(Image); clear('Image');




if(~isfield(outParams,'fieldmap'))
outParams.fieldmap = zeros(P ,Q, slices);
outParams.species(1).amps = zeros(P ,Q, slices);
outParams.species(2).amps = zeros(P ,Q, slices);
outParams.r2starmap = zeros(P ,Q, slices);
outParams.mask = zeros(P ,Q, slices);
end

mask = zeros(size(img1));
mag_im = abs(img1);     
max_im = max(mag_im(:));
mag_im = mag_im / max_im;
id  = mag_im > graythresh(mag_im) ;
mask(id) = 1;

step = 20;

 for n = 1 :step:slices
    
   start = n
   stop = n+step-1
   if stop > slices
       stop = slices;
   end
    
%     %%% QPBO Algorithm   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        im = cat(4,img1(:,:,start:stop),img2(:,:,start:stop),img3(:,:,start:stop),img4(:,:,start:stop)); 
%      im = cat(4,img1(:,:,start:stop),img2(:,:,start:stop),img3(:,:,start:stop)); 

%     %%% 3 POINT MultiSeedRegionGrowing Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       im = cat(4,img1(:,:,start:stop ),img2(:,:,start:stop ),img3(:,:,start:stop )); 
    
   imDataParams.images = reshape(im, [P,Q,stop-start+1,1,size(im,4)]);
   imDataParams.mask = mask(:,:,start:stop ); 
   clear im;

   
   if(~isempty(find(imDataParams.mask, 1)))
    
%     %%% QPBO Algorithm   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       out = fw_i3cm1i_3pluspoint_berglund_QPBO( imDataParams, algoParams );

     %%% 3 POINT MultiSeedRegionGrowing Algorithm    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      out = fw_i3cm0i_3point_berglund( imDataParams, algoParams );
    
     outParams.fieldmap(:,:,start:stop) = permute(flip(flip(out.fieldmap,1),2),[2,1,3]);
     outParams.species(1).amps(:,:,start:stop) = permute(flip(flip(out.species(1).amps,1),2),[2,1,3]);
     outParams.species(2).amps(:,:,start:stop) = permute(flip(flip(out.species(2).amps,1),2),[2,1,3]);
   
     
    %%% QPBO Algorithm   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % outParams.r2starmap = out.r2starmap;
    % outParams.mask(:,:,start:stop) = permute(flip(flip(out.mask,1),2),[2,1,3]);
   end
   
 end   
toc

outParams.fieldmap = -outParams.fieldmap;   %%% Negative sign added to correct 3Point field sign
outParams.mask = permute(flip(flip(mask,1),2),[2,1,3]);
outParams.nPy = nPy;
outParams.dz = z_voxel;
save([pathname, fileno,'_outParams_JULIA'],'outParams');

toc;


%%  CALLIING DYNAMIC SHIMMING TOOL FOR CALCULATING SHIMS

CMT_Dynamic_Shimming_GUI(outParams);






