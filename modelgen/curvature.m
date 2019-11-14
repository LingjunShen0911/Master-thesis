function [z1_grid, z2_grid] = curvature(Points, Modellgrenzen, options)
%Curvature finding the curvature of PWLModel in x and y direction
% 
% Author: Gerald Koroa
% Last Edit: 21.10.16

    % --- Finding derivation 1 and 2 of x2-y plane ----------------------------
    Pdiff_x2 = diff(Points,1);
    Pdiff_x2(:,4) = Pdiff_x2(:,3)./Pdiff_x2(:,2);
    Pdiff2_x2 = diff(Pdiff_x2,1);

    % --- Processing matrix ---------------------------------------------------
    Rem = find(Pdiff2_x2(:,1));
    Pdiff2_x2(:,3) = [];
    Pdiff2_x2(Rem,:) = NaN;
    Pdiff2_x2 = [NaN(1,3);Pdiff2_x2;NaN(1,3)];
    Pdiff2_x2(:,1) = Points(:,1);
    Pdiff2_x2(:,2) = Points(:,2);

    Pdiff2_x2(:,3) = Pdiff2_x2(:,3)/Modellgrenzen(6);

    Rem = find(isnan(Pdiff2_x2(:,3)));
    Pdiff2_x2(Rem,:) = [];

    % --- Sorting data to x2 --------------------------------------------------
    [~,ind] = sort(Points(:,2));
    Points_x2= Points(ind,:);

    % --- Finding derivation 1 and 2 of x1-y plane ----------------------------
    Pdiff_x1 = diff(Points_x2,1);
    Pdiff_x1(:,4) = Pdiff_x1(:,3)./Pdiff_x1(:,1);
    Pdiff2_x1 = diff(Pdiff_x1,1);

    % --- Processing matrix ---------------------------------------------------
    Rem = find(Pdiff2_x1(:,1));
    Pdiff2_x1(:,3) = [];
    Pdiff2_x1(Rem,:) = NaN;
    Pdiff2_x1 = [NaN(1,3);Pdiff2_x1;NaN(1,3)];
    Pdiff2_x1(:,1) = Points_x2(:,1);
    Pdiff2_x1(:,2) = Points_x2(:,2);

    Pdiff2_x1(:,3) = Pdiff2_x1(:,3)/Modellgrenzen(3);

    Rem = find(isnan(Pdiff2_x1(:,3)));
    Pdiff2_x1(Rem,:) = [];

    % --- Reshape to Griddata -------------------------------------------------
    tx1 = Modellgrenzen(2):Modellgrenzen(3):Modellgrenzen(1);
    tx2 = Modellgrenzen(5):Modellgrenzen(6):Modellgrenzen(4);
    [x1q,x2q] = meshgrid(tx1,tx2);
    % Using extrapolation to predict the value of edge-points.
    F1 = scatteredInterpolant(Pdiff2_x1(:,1),Pdiff2_x1(:,2),Pdiff2_x1(:,3));
    F2 = scatteredInterpolant(Pdiff2_x2(:,1),Pdiff2_x2(:,2),Pdiff2_x2(:,3));

    z1_grid = abs(F1(x1q,x2q));
    z2_grid = abs(F2(x1q,x2q));

    % --- Plot ----------------------------------------------------------------
    if (options.plot)
        
        scrsz=get(groot,'ScreenSize');
        
        fh = figure(2);
        set(fh,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
        mesh(x1q,x2q,z1_grid,'LineWidth',2);
        view(0,90);
        xlabel('U_{DS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
        ylabel('U_{GS}(V)', 'FontSize', 18, 'FontName', 'Helvetica') ;   
        title('2nd Derivative of x-z plane');
        colormap(jet);
%         hold on;
%         plot3(Points(:,1),Points(:,2),Points(:,3));
%         hold off;

        fh = figure(3);
        set(fh,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2])
        mesh(x1q,x2q,z2_grid,'LineWidth',2);
        view(0,90);    
        xlabel('U_{DS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
        ylabel('U_{GS}(V)', 'FontSize', 18, 'FontName', 'Helvetica');
        title('2nd Derivative of y-z plane');
        colormap(jet);
%         hold on;
%         plot3(Points(:,1),Points(:,2),Points(:,3),'.','15');
%         hold off;
    end
end

