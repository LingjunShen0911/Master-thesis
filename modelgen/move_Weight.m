function weight = move_Weight(z1_grid, z2_grid, PWLModel, Bias, vertexindex, Modellgrenzen)
%MOVE_WEIGHT calculating directional weighting of point
% 
% Inputs:
%Bias - Higher the number, the less effect weighting have
% 
% Author: Gerald Koroa
% Last Edit: 21.10.16

    % x1-richtung
    m1 = floor((PWLModel(vertexindex,2)-Modellgrenzen(5))/Modellgrenzen(6)+1e-8)+1;
    n1 = floor((PWLModel(vertexindex,1)-Modellgrenzen(2))/Modellgrenzen(3)+1e-8)+1;
%     A1 = sum(z1_grid(m1:m1+1,:),1);
%     Left = Bias+(sum(A1(1,1:n1))/n1);
%     Right = Bias+(sum(A1(1,n1+1:end))/(size(A1,2)-n1));
    if m1 > size(z1_grid,1)
        m1 = size(z1_grid,1);
    end
    if n1 > size(z1_grid,2)
        n1 = size(z1_grid,2);
    end
    
    Left = Bias+(sum(z1_grid(m1,1:n1))/n1);
    Right = Bias+(sum(z1_grid(m1,n1:end))/(size(z1_grid(1,n1:end),2)));

    % x2-richtung
    m2 = floor((PWLModel(vertexindex,1)-Modellgrenzen(2))/Modellgrenzen(3)+1e-8)+1;
    n2 = floor((PWLModel(vertexindex,2)-Modellgrenzen(5))/Modellgrenzen(6)+1e-8)+1;
%     A2 = sum(z2_grid(:,m2:m2+1),2);
%     Down = Bias+(sum(A2(1:n2,1))/n2);
%     Up = Bias+(sum(A2(n2+1:end,1))/(size(A2,1)-n2));
    if m2 > size(z2_grid,1)
        m2 = size(z2_grid,1);
    end
    if n2 > size(z2_grid,2)
        n2 = size(z2_grid,2);
    end
    
    Down = Bias+(sum(z2_grid(1:n2,m2))/n2);
    Up = Bias+(sum(z2_grid(n2:end,m2))/(size(z2_grid(n2:end,m2),1)));

    weight(1) = Left/(Left+Right);
    weight(2) = Right/(Left+Right);
    weight(3) = Up/(Up+Down);
    weight(4) = Down/(Up+Down);
end