%%%% This file will perform FW sepration with Berglund's QPBO   %%%%%%
%%%%% Fat Water Separation N slices at a time
%%%%% Change fileno to the required multiecho imaage files, ie
%%%%% either 'Sub1_DS' , 'Sub1_NS', 'Sub2_NS'  or 'Sub2_DS'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


tic;

pathname = '../data_output/';

fileno = 'Sub1_DS';

% uncomment individual line below to process other data sets
%fileno = 'Sub1_NS';
%fileno = 'Sub2_DS';
%fileno = 'Sub2_NS';

TE = [0.99 , 2.39, 3.79, 5.19]* 1e-3;   % 4 echo for QPBO

%%%% imDataParams %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

imDataParams.TE = TE;
imDataParams.FieldStrength = 3;
imDataParams.voxelSize  = [2.5 , 2.5 , 2.5 ]; 
imDataParams.PrecessionIsClockwise = 1;

%%%%% algoParams   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 4.70;
algoParams.species(1).relAmps = 1;
algoParams.species(2).name = 'fat';
algoParams.species(2).frequency = [0.90, 1.30, 1.60, 2.02, 2.24, 2.75, 4.20, 5.19, 5.29];
algoParams.species(2).relAmps   = [  88,  642,   58,   62,   58,    6,   39,   10,   37];
algoParams.process_in_3D = true;

%%%%%% outParams %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
outParams = [];   
outParams.species(1).name = 'water';
outParams.species(2).name = 'fat';


%%%%% Building imDataParams.images from PYTHON RECONS   %%%%%%%%%%%%%%%%%%%%%%

for echo = 1 : size(TE,2)
    
filename  = [fileno,'_echo',num2str(echo),'_PYTHONREC'];
load([pathname filename]);

[ slices, P, Q ] = size(img);


%%%%%%%% Dimension for splitting Fat/water separation %%%%%%%%%%

%%%% For Axial
    img = permute(img,[ 2 3 1]);
% %%%% For Coronal
%     img = permute(img,[ 2 1 3]);
%%%% For Sagittal
% data.img = permute(img,[ 3 1 2]);

im(:,:,:,echo) = img;

if echo == 1
    mask = zeros(size(img));
    mag_im = abs(img);
    max_im = max(mag_im(:));
    mag_im = mag_im / max_im;
    id  = mag_im > graythresh(mag_im);
    mask(id) = 1;
    clear mag_im id
end

clear img 
      
str = sprintf('echo number %d done',echo);
disp(str)

end


if(~isfield(outParams,'fieldmap'))
outParams.fieldmap = zeros(P ,Q, slices);
outParams.species(1).amps = zeros(P ,Q, slices);
outParams.species(2).amps = zeros(P ,Q, slices);
outParams.r2starmap = zeros(P ,Q, slices);
outParams.mask = zeros(P ,Q, slices);
end


imDataParams.images = reshape(im, [size(im,1),size(im,2),size(im,3),1,size(im,4)]);
imDataParams.mask = mask;
clear im;

%%%%%%%%%%%% QPBO FW Separation %%%%%%%%%%%%%%%%%%%%%%%%%

step = 20;
im_section = imDataParams;
cat_dim = 3;  %%% Dimension to break up the processing.
% cat_dim = 2;  %%% Dimension to break up the processing.
nstart = 1 ;
nstop  = size(imDataParams.images,cat_dim);


for n = nstart:step:nstop
    
    start = n
    if n+step-1 > size(imDataParams.images,cat_dim)
        stop = size(imDataParams.images,cat_dim)
    else
        stop = n+step-1
    end   
    
    if cat_dim == 3  
     im_section.images = imDataParams.images(:,:,start:stop,:,:);
     im_section.mask = imDataParams.mask(:,:,start:stop);
    elseif  cat_dim == 2
     im_section.images = imDataParams.images(:,start:stop,:,:,:);
     im_section.mask = imDataParams.mask(:,start:stop,:);
    end
    
    if(~isempty(find(im_section.mask, 1)))     
        out = fw_i3cm1i_3pluspoint_berglund_QPBO( im_section, algoParams );
           
        outParams.fieldmap(:,:,start:stop) = out.fieldmap;
        outParams.species(1).amps(:,:,start:stop) = out.species(1).amps;
        outParams.species(2).amps(:,:,start:stop) = out.species(2).amps;
        outParams.r2starmap(:,:,start:stop) = out.r2starmap;
        outParams.mask(:,:,start:stop) = out.mask;   
    end
    
        
    str = sprintf('section %d done',n);
    disp(str)
    
    clear out   
end

crop_range_LR = 1:P;
crop_range_AP = 1:P;

outParams.species(1).amps = outParams.species(1).amps(crop_range_LR,crop_range_AP,:);
outParams.species(2).amps = outParams.species(2).amps(crop_range_LR,crop_range_AP,:);
outParams.fieldmap = outParams.fieldmap(crop_range_LR,crop_range_AP,:);
outParams.mask = mask(crop_range_LR,crop_range_AP,:);
outParams.r2starmap = outParams.r2starmap(crop_range_LR,crop_range_AP,:);

save([pathname, fileno, '_', 'outParamsQPBO_PYTHON'], 'outParams', '-v7.3');

toc;