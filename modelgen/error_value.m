function [MeanSE, MaxSE]= error_value(Points, P_Referenz)
% calculate error between PWLModel and referenceModel
%
% Author: Gerald Koroa
% Written: Oktober 2016
% Last update: 28.10.2016
% 
% refine PWLModel points using interpolation
F = scatteredInterpolant(Points(:,1),Points(:,2),Points(:,3));
P_new = F(P_Referenz(:,1),P_Referenz(:,2));

% tx=-5:0.5:20;
% ty=-5:0.5:20;
% [x,y]=meshgrid(tx,ty);
% 
% P_Grid = griddata(Points(:,1), Points(:,2), Points(:,3), x, y);
% P_Grid(find(isnan(P_Grid)))=[];
% P_new = P_Grid';

% using Mean Square Error(MSE)
dP=(P_Referenz(:,3)-P_new).^2;
MeanSE=sum(dP)/size(dP,1);
[Value,Ind]=max(dP);
%  Point = [P_Referenz(Ind,1),P_Referenz(Ind,2)];
%  MaxSE=struct('Value', Value, 'Point', Point);
MaxSE=[P_Referenz(Ind,1),P_Referenz(Ind,2),Value];
end