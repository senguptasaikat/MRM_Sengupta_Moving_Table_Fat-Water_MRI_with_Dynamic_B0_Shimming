function [Tpom T ] = Rotation_Stack(orient,ang_ap,ang_rl,ang_fh,Patient_Pos)

A = strread(Patient_Pos,'%s');

Pat_Ori = A{3};
Pat_Pos = strcat(A{1},A{2});


%%%%%%%%%  Rotation Matrices 

%%% Slice Orientation Matrix (SOM)
%%% Image/Matlab System (SEV , South , East ,View  or Row, Col, Depth) to 
%%% Angulated Patient Axis System (L'P'H', Left, Posterior, Head)

if orient == 1              %%% Axial
    
    Tsom = [ 0 1 0;
             1 0 0; 
             0 0 1];
         
elseif  orient == 3        %%%% Coronal
    
    Tsom = [ 0 1 0;
             0 0 1; 
            -1 0 0];
         
elseif  orient == 2         %%% Sagittal
    
    Tsom = [ 0 0 -1;
             0 1  0; 
            -1 0  0];
end


%%%%% Angulation Matrix (ANG)
%%%%  Angulated Patient Coordinate System ( L'P'H') to
%%%%  Non Angulated Patient Coordinate System (L,P,H)

Trl = [ 1      0            0 ;
        0 cos(ang_rl) -sin(ang_rl);
        0 sin(ang_rl)  cos(ang_rl)];
    
Tap = [ cos(ang_ap)   0   sin(ang_ap) ; 
           0          1        0     ; 
       -sin(ang_ap)   0   cos(ang_ap)];
   
Tfh = [ cos(ang_fh)   -sin(ang_fh)  0 ;
        sin(ang_fh)   cos(ang_fh)   0 ;
            0             0         1];
        
Tang = Trl * Tap * Tfh;
   
   
 %%%%% Patient Orientation Matrix %%%%%%%%%%%%%%%%%%%%%%%5
 %%%%% Non Angulated Patient Coordinate System (L,P,H) to
 %%%%% Magnet/Shim Coordinate system (X,Y,Z)
 
 %%% Patient Orientation :Supine 

 
 if strcmp(Pat_Ori,'Supine')== 1        %%% Patient Orientation :Supine 
     
     Tpo = [ 1  0  0;
             0  1  0;
             0  0  1];  

 elseif strcmp(Pat_Ori,'Prone')== 1     %%% Patient Orientation : Prone
 
     Tpo = [ -1  0 0;
              0 -1 0;
              0  0 1];
          
 end
%  
%   %%% Patient Orientation : Right Decubitus
%   
%  Tpo = [ 0 -1 0;
%          1  0 0;
%          0  0 1];
%    
%   %%% Patient Orientation :Left Decubitus
%   
%  Tpo = [ 0 1 0;
%         -1 0 0;
%          0 0 1];
%    
if strcmp(Pat_Pos,'FeetFirst')== 1       %   %%% Patient Position :Feet First  
  
      Tpp= [ 0 -1 0;
            -1  0 0;
             0  0 1];
  
elseif  strcmp(Pat_Pos,'HeadFirst')== 1       %   %%% Patient Position :Head First  
   
  Tpp = [ 0 1  0;
         -1 0  0;
          0 0 -1];
      
end
                 
         
   Tpom = Tpo * Tpp;          
   Tpom = inv(Tpom);      
        
                
                          
  T =  Tpom * Tang *  Tsom ;
    