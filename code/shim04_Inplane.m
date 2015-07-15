function coefficients = shim04_Inplane(firstSlice,lastSlice,dataSpace,fieldmap,index,degen_index,tp_shim_index)
% COEFFICIENTS = SHIM04(FIRSTSLICE,LASTSLICE,DATASPACE)
% returns coefficients of field expansion by multiple linear 
% regression ( 4 terms )
% input:
%   firstSlice: number of first slice whose data used in regression
%   lastSlice: lastslice number
%   dataSpace: data spacing in image plane
%   output: coefficients of field expansion

% Created by Yansong Zhao, Yale University June 01. 2001
% This is a function of Shimming Toolbox

%%%%% Modified By Saikat Sengupta, Vanderbilt University, Dec 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global  roi R Theta Phi;
                              
[Np,Nr,Ns]=size(fieldmap);

% number of data used in shimming
nSlice = size([firstSlice:lastSlice],2);
nPhase = size([1:dataSpace:Np],2);
nRead = size([1:dataSpace:Nr],2);

%%%%% Define slice limits for R,Theta,Phi since they are only Np x Nr x 3
 
   firstSlice_RPT = 1;
   lastSlice_RPT = 3;
   index_RPT = 2;
   
   if index == 1
      firstSlice_RPT = 2;        
   elseif index == Ns
      lastSlice_RPT = 2;
   end


%%%%%%%%%%%%%%%%%%% Data used in shimming inplane %%%%%%%%%%%
    
    shimField_ip = reshape(fieldmap(1:dataSpace:Np,1:dataSpace:Nr,...
        index),nRead*nPhase,1);

    shimRoi_ip = reshape(roi(1:dataSpace:Np,1:dataSpace:Nr,...
        index),nRead*nPhase,1);
    
    shimR_ip = reshape(R(1:dataSpace:Np,1:dataSpace:Nr,...
        index_RPT),nRead*nPhase,1);

    shimTheta_ip = reshape(Theta(1:dataSpace:Np,1:dataSpace:Nr,...
       index_RPT),nRead*nPhase,1);

    shimPhi_ip = reshape(Phi(1:dataSpace:Np,1:dataSpace:Nr,...
       index_RPT),nRead*nPhase,1);
    
    
  %%%%%%%%%%%%%%%%%%% Data used in shimming through plane %%%%%%%%%
    
    shimRoi_tp = reshape(roi(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice:lastSlice),nRead*nPhase*nSlice,1);
    
    shimR_tp = reshape(R(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice_RPT:lastSlice_RPT),nRead*nPhase*nSlice,1);

    shimTheta_tp = reshape(Theta(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice_RPT:lastSlice_RPT),nRead*nPhase*nSlice,1);

    shimPhi_tp = reshape(Phi(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice_RPT:lastSlice_RPT),nRead*nPhase*nSlice,1);
    
 coefficients_ip = zeros(4,1);
 coefficients_tp = zeros(4,1);


%%%%%%% Fittting In Plane Shims (Axial) %%%%%%%%%%%%%%%%%%%%%

index_ip = find(shimRoi_ip == 1);   %% Indices of masked areas where shimROI ==1 
ip_shim_index = find(degen_index);

if isempty(index_ip)
   coefficients = zeros(4,1);
    return
end


for k = 1 : size(index_ip,1)
    Temp = onepoint04(shimR_ip(index_ip(k)),shimTheta_ip(index_ip(k))...
        ,shimPhi_ip(index_ip(k)) );
     
    X_ip(k,:) = Temp(ip_shim_index);       
    Y_ip(k) = shimField_ip(index_ip(k));
end

coefficients_ip(ip_shim_index) = regress(Y_ip',X_ip); 
coeff_Z0 = coefficients_ip(1);




%%%%%%% Fittting Through Plane Shims, after Inplane shims 
%%%%%%% To compensate for through plane component of inplane shims
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i = 1 : size(tp_shim_index,2)
    if tp_shim_index(i) == 1
        tp_shim_index(i) = 3;           %%% x
    elseif tp_shim_index(i) == 2
        tp_shim_index(i) = 4;           %%% y
    elseif tp_shim_index(i) == 3
        tp_shim_index(i) = 2;           %%% z
    end
end

tp_shim_index = sort(tp_shim_index,2,'ascend');
tp_shim_index = [ 1 ,   tp_shim_index ];     
             
       
for i = 1:nSlice
    shimField_after_ip(:,:,i) = aftershim04(coefficients_ip,firstSlice+i-1,i,fieldmap);
end

Avg = mean(nonzeros(shimField_after_ip));
Stdev = std(nonzeros(shimField_after_ip));
id =  shimField_after_ip > (Avg + 2.5* Stdev) | shimField_after_ip < (Avg - 2.5*Stdev) ;
shimField_after_ip(id) = 0;

index_tp = find(shimRoi_tp == 1);   %% Indices of masked areas where shimROI ==1 

for k = 1 : size(index_tp,1)
   Temp = onepoint04(shimR_tp(index_tp(k)),shimTheta_tp(index_tp(k))...
        ,shimPhi_tp(index_tp(k)));
     
    X_tp(k,:) = Temp(tp_shim_index);      %% Axial
     
    Y_tp(k) = shimField_after_ip(index_tp(k));
end


coefficients_tp(tp_shim_index)= regress(Y_tp',X_tp); %% Axial


%%%%%%%%% Saving the Final Slicewise Field with Through Plane  %%%%%%%%%%%%%%%%

% fieldmap_copy = zeros(Np,Nr,Nsl);
% fieldmap_copy(:,:,firstSlice:lastSlice) = shimField_after_ip;
% shimField_after_tp = zeros(Np,Nr,3);
% for i = 1:nSlice
%     shimField_after_tp(:,:,i) = aftershim04(coefficients_tp,i+firstSlice-1,fieldmap_copy);
% end
% save( sprintf('GUI Signal Sim/ShimDS1_%d',index), 'shimField_after_tp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


coefficients =  coefficients_ip  +  coefficients_tp;

coefficients(1) = coeff_Z0 ;



% END of shim04.m
% This is a function of Shimming Toolbox
