
function [area]=amar_area(x,y)
%#codegen
%% Code to calculate the area of contours


%% Discription
% This function takes care of the fact that area is positive when curve moves in clock wise
% direction -- for outer contour and
% area is negative for curve in anti clock wise direction -- for voids.

% Algorithm is given by
% $$ area=\sum{\frac{1}{2}(x_i+x_{i+1})(y_i-y_{i+1})}  $$
% X is the input abscissa,
% Y is the input ordinates and

%% Code
n=length(x); % number of nodes in the current contour
 x=reshape(x,[1,n]);
y=reshape(y,[1,n]);           
xf=[x(2:n),x(2)];yf=[y(2:n),y(2)];
area=-0.5*(x(1:end-1)+xf(1:end-1))*(y(1:end-1)-yf(1:end-1))';    
                              
                
               
                
               
            
   

