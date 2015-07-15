function x_line=onepoint04(r,theta,phi)
% X_LINE = ONEPOINT04(R,THETA,PHI) returns vector
% r^n*Pnm(cos(theta))*cos(m*phi)|sin(m*phi) of one point in spherical coordinate.
%         r: radius; theta: angle displacement ; phi: azimuth
%  Z0 and all 1st order spherical harmonics
%  ------------------------
%  #    n     m   shim
%  01   0     0    Z0
%  02   1     0    Z
%  03   1     1c   X
%  04   1     1s   Y
%  ------------------------
%  c : cos term , s: sin term

%   Created by Yansong Zhao, Yale University Oct. 1999
%   This is a function of Shimming Toolbox

s=sin(theta);
c=cos(theta);

x_line(1)=1;                      %z0
x_line(2)=r*c;                    %z
x_line(3)=r*s*cos(phi);           %x
x_line(4)=r*s*sin(phi);           %y

% END of onepoint04.m  
% This is a function of Shimming Toolbox
