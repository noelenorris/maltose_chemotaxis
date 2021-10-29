
% Goodness of fit measure J_x
%
% Noele Norris
%

function [ distance ] = square_measure(dist1, dist2)
distance = sum((dist1-dist2).^2);
   if(isnan(distance))
       distance = 100;
   end
end