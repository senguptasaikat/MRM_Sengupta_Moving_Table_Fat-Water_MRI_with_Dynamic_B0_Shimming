function coefficients = shim04(firstSlice,lastSlice,dataSpace,fieldmap)
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
                              
[Np,Nr,Nsl]=size(fieldmap);

% number of data used in shimming
nSlice = size([firstSlice:lastSlice],2);
nPhase = size([1:dataSpace:Np],2);
nRead = size([1:dataSpace:Nr],2);

%%%%%%%%%%%%%%%%%%% Data used in shimming
    
    shimField = reshape(fieldmap(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice:lastSlice),nRead*nPhase*nSlice,1);

    shimRoi = reshape(roi(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice:lastSlice),nRead*nPhase*nSlice,1);
    
    shimR = reshape(R(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice:lastSlice),nRead*nPhase*nSlice,1);

    shimTheta = reshape(Theta(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice:lastSlice),nRead*nPhase*nSlice,1);

    shimPhi = reshape(Phi(1:dataSpace:Np,1:dataSpace:Nr,...
        firstSlice:lastSlice),nRead*nPhase*nSlice,1);


index = find(shimRoi == 1);   %% Indices of masked areas where shimROI ==1 

if isempty(index)
   coefficients = zeros(4,1);
    return
end


for k = 1 : size(index,1)
    X(k,:) = onepoint04(shimR(index(k)),shimTheta(index(k))...
        ,shimPhi(index(k)));
    Y(k) = shimField(index(k));
end

coefficients= regress(Y',X);

if nSlice <= 3
     if firstSlice == 1 && lastSlice == 2
         slice = 1;
     else
         slice = firstSlice + 1;
     end
     
     %%% Since Z0 is calculated separately later on in DynShim_3sl,
     %%% coefficients(1) is changed. It doesnt affect Global Shim.
     coefficients(1) = median(nonzeros(fieldmap(:,:,slice)));
 end






% END of shim04.m
% This is a function of Shimming Toolbox