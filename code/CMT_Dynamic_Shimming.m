
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%     Dynamic 1st Order Shimming for COMBI Scan              

%%% Inputs : COMBI Fieldmaps, Patient Geometry
%%% Outputs : Slicewise X,Y,Z,Z0 shims, interpolated to 
%%% required slice thinkness 

%%% Writen By : Saikat Sengupta, Vanderbilt University, 2013

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global roi R Theta Phi;


roi = handles.fieldmap_Thresholded ~= 0;

[ Np, Nf, Ns ] = size(handles.fieldmap_Thresholded);
 
  
%%%%%%%%%%% Reading in Parameters and Creating the Coordinates %%%%%%%%%%%

if Nf > Np                    % Sometimes Nf < Np cos of less than 100% encoding.
    Np = Nf;
elseif Np > Nf
    Nf = Np;
end


block_ori = 1;                % axial

% FOV_AP = round(handles.info.recon_resolutions.vals(1) * handles.info.voxel_sizes.vals(1))/10;
% FOV_RL = round(handles.info.recon_resolutions.vals(2) * handles.info.voxel_sizes.vals(2))/10;
% FOV_FH = 3 * round(handles.info.recon_resolutions.vals(3) * handles.info.voxel_sizes.vals(3))/10;
% Slice_Thk = handles.info.voxel_sizes.vals(3)/10; 



% In R5 there might not be any SIN file, We have to use .Data/list ,
% therefore loadSIN does not work. Neither does loadPARREC as there is no
% PARREC file.FOV values are hardcoded in COMBI_Dynamic_Shimming for now

FOV_AP = 40;
FOV_RL = 40;
FOV_FH = 3.6;
Slice_Thk = 1.2;


% FOV_FH = 3* handles.dz;
% Slice_Thk = handles.dz;

%%%%%%%%%%%%%%%%%%%% Angulation and Orientation   %%%%%%%%%%%%%%%%%%%% 
 
ang_ap = 0;
ang_fh = 0;
ang_rl = 0;
 
Patient_Pos = 'Head First Supine' ;
 
[Tpom, Geo_T] = Rotation_Stack(block_ori,ang_ap,ang_rl,ang_fh,Patient_Pos);
 

%%%%%%   block_ori == 1 axial  

ip_degen_index =  [1;0;1;1 ];
tp_degen_index = 3;
              
      
%%%%%%%%%%%%%%%%%%%%%%% Processing %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% Coordinate System for Gradients Vanderbilt 3 Tesla  %%%%%%%%%%%%%       
% 
%         | +x
%         |
%         |
%         /--------> +y
%        /
%       /
%      /
%    +z
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%  block_ori == 1  axial      
      
FOV_ROW = FOV_AP; 
FOV_COL = FOV_RL; 
FOV_DEPTH = FOV_FH; 

row = linspace(-FOV_ROW/2,FOV_ROW/2,Np);
col = linspace(-FOV_COL/2,FOV_COL/2,Nf);
% depth = linspace(- FOV_DEPTH/2 + Slice_Thk/2, FOV_DEPTH/2-Slice_Thk/2,3);
% depth = [ 0.25 , 0 , -0.25];
depth = [FOV_AP/Np , 0 , -FOV_AP/Np];   %%% Images are reconstructed at isotropic resolution


[COL,ROW,DEPTH] = meshgrid(col,row,depth);

for row = 1 : Np
    for col = 1 : Nf
        for depth = 1 : 3
          Mag = Geo_T * [ROW(row,col,depth);COL(row,col,depth);DEPTH(row,col,depth)];
          Mag_X(row,col,depth) = Mag(1);
          Mag_Y(row,col,depth) = Mag(2);
          Mag_Z(row,col,depth) = Mag(3);
        end
    end
end


[Phi,Theta_complement,R] = cart2sph(Mag_X,Mag_Y,Mag_Z);
Theta = pi/2 - Theta_complement;

clear Mag_X;
clear Mag_Y;
clear Mag_Z;
clear ROW;
clear COL;
clear DEPTH;
 

%%%%%%%%%%%%%%%%%%%%%%%%%%%% Shimming  %%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
handles.A = zeros(Np,Nf,Ns);

for index=1:Ns
   
   index   
   leftSlice = index-1;
   rightSlice = index+1;
   
   if index == 1
      leftSlice = 1;
   elseif index == Ns
      rightSlice = Ns;
   end
    
   if index ~= 1
       if size(nonzeros(handles.fieldmap_Thresholded(:,:,index)),1) < 0.003*Np*Nf
           handles.coeffs(:,index) =  zeros(4,1);
           continue;
       end
   end
   
   index_RPT = 2;
   
   
    %%%%%% Shimming 1st Order only
    
%         if strcmp(handles.shim_type,'Dynamic_3Sl')
%         handles.coeffs(:,index) = shim04(leftSlice,rightSlice,1,fieldmap);
%         elseif strcmp(handles.shim_type,'Dynamic_IP')
   handles.coeffs(:,index) = shim04_Inplane(leftSlice,rightSlice,1,handles.fieldmap_Thresholded,index,ip_degen_index,tp_degen_index);   
%         end
            
   handles.A(:,:,index) = aftershim04(handles.coeffs(:,index),index,index_RPT,handles.fieldmap_Thresholded);
 
    
end

toc
handles.coeffs(2:end,:) = -handles.coeffs(2:end,:);    %%%% Corrections
handles.Freq_Offset  = handles.coeffs(1,:);

%%%%%%%%%%% Display Shimmed Fieldmaps %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


axes(handles.axes5);
imagesc(squeeze(handles.A(:,:,round(Ns/2))),[handles.fmap_min handles.fmap_max]);

axes(handles.axes6);
imagesc(squeeze(handles.A(32,:,:,1,1)),[handles.fmap_min handles.fmap_max]);

axes(handles.axes7);
imagesc(squeeze(handles.A(:,32,:,1,1)),[handles.fmap_min handles.fmap_max]);

 
axes(handles.axes8);

x = 1 :Ns;
y1 =  handles.coeffs(3,:);
y2 =  handles.coeffs(4,:);
y3 =  handles.coeffs(2,:);
y4 =  handles.Freq_Offset(1,:)/10;

plot(x,y1,'r'); hold on ;
plot(x,y2,'g'); hold on ;
plot(x,y3,'b'); hold on ;
plot(x,y4,'m'); ylim([ -50 50]);

ylabel(' Hz / cm  ;   Hz / 10  ');
legend(' X', 'Y','Z','Z0');
xlabel(' slice'); 



%%%%%%%%%%%%%%%%%%%%%%% Conversion to mT/m %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%To convert from Hz/cm to mT/m divide by 425.8

% First Order       
%   01 Z0           
%   02 Z            
%   03 X            
%   04 Y            

handles.coeffs_1st([1,2,3],:) = handles.coeffs([2,3,4],:)./425.8 ;  % Hz/cm - mT/m 


%%%%%% Filter Z direction Coeffs

handles.coeffs_1st(1,:) = medfilt1(handles.coeffs_1st(1,:),3);

id1 = find(handles.coeffs_1st > 0.1);
handles.coeffs_1st(id1) = 0.1; 
id2 = find(handles.coeffs_1st < -0.1);
handles.coeffs_1st(id2) = -0.1; 



%%%% Display in GUI
   
set(handles.text5,'String',num2str(handles.Freq_Offset(1))); % F0
set(handles.text6,'String',num2str(handles.coeffs_1st(1,1)));  %%% Z
set(handles.text7,'String',num2str(handles.coeffs_1st(2,1)));  %%% X
set(handles.text8,'String',num2str(handles.coeffs_1st(3,1)));  %%% Y


guidata(hObject, handles);




 
