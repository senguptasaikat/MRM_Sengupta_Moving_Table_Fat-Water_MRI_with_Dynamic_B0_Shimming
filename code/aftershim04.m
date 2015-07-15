function Y = aftershim04(coefficients,slice,slice_RPT,fieldmap)
% Y = AFTERSHIM04(COEFFICIENTS,SLICE)
% returns field strength ( in Hz ) of one slice after 4-coil-shim
% shim terms : Z0, Z, X, Y
% input:
%   coefficients: Coefficients of field expansion
%   slice : the number of slice which is being shimmed
% output:
%   Y : field strength after shim

%   Created by Yansong Zhao, Yale University June 01. 2001
%   This is a function of Shimming Toolbox

%%%%% Modified By Saikat Sengupta, Vanderbilt University, Dec 2005

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global roi R Theta Phi;

[Np,Nr,Nsl]=size(fieldmap);
firstSlice = 1;
lastSlice = 1;
dataSpace = 1;

nSlice = size([firstSlice:lastSlice],2);
nPhase = size([1:dataSpace:Np],2);
nRead = size([1:dataSpace:Nr],2);


    
    shimField = reshape(fieldmap(1:dataSpace:Np,1:dataSpace:Nr,...
       slice),nRead*nPhase,1);

    shimRoi = reshape(roi(1:dataSpace:Np,1:dataSpace:Nr,...
       slice),nRead*nPhase,1);
    
    shimR = reshape(R(1:dataSpace:Np,1:dataSpace:Nr,...
        slice_RPT),nRead*nPhase,1);

    shimTheta = reshape(Theta(1:dataSpace:Np,1:dataSpace:Nr,...
       slice_RPT),nRead*nPhase,1);

    shimPhi = reshape(Phi(1:dataSpace:Np,1:dataSpace:Nr,...
        slice_RPT),nRead*nPhase,1);


index = find(shimRoi == 1);

Y = zeros(1,Nr*Np);

for k = 1: size(index,1)
    
    X = onepoint04(shimR(index(k)),shimTheta(index(k))...
        ,shimPhi(index(k)));
    
    Y(index(k)) = shimField(index(k)) - X*coefficients;   %%% Corrections for first order shims ?
    
end

Y = reshape(Y,Np,Nr);

% END of aftershim09.m
% This is a function of Shimming Toolbox