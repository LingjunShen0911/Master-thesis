function [ PWLModel, PWLCharacteristics, variantName] = simulatedannealing(paths)
% SIMULATEDANNEALING Optimization of a PWL Model by randomly moving vertices
% 
% Literature: Modellierung komplexer Prozesse durch naturanaloge Verfahren
%             2012 (2. Edition), page 89, Authors: Klüver, Klüver, Schmidt
%
% Inputs:
% - filename_referenceModel (XML file): e.g. 'models/reference.xml'
% - filename_initialModel (XML file): e.g. 'models/model.xml'
%
% Requirements:
% - MATLAB version: R2014b (tested)
% 
% - error_value.m
% - move_Point.m
% - plotpoints.m
% - read_xml.m
% - Starttemperatur.m
% - curvature.m
% - move_Weight.m
%
% Example:
% - simulatedann50
% Author: Lukas Lee (Lukas.Lee@ims.uni-hannover.de), Gerald Koroa
% Last Edit: 28.10.16



% Load Settings
if (exist(paths.settingsPath, 'file') ~= 2)
    error('No settings file available!');
else
    try
        clear(paths.settingsPath); %Clear old cache
        run(paths.settingsPath);
    catch
        error('Invalid settings!');
    end
end

% Load Reference Model
if exist([paths.referenceModelPath,'.mat']) == 2
    load([paths.referenceModelPath,'.mat']);
else
    referenceModel.Points=read_xml([paths.referenceModelPath,'.xml']);
    
    %Sort data order(important for Modified SA)
    [~,I]=sort(referenceModel.Points(:,2));
    referenceModel.Points=referenceModel.Points(I,:);
    [~,I]=sort(referenceModel.Points(:,1));
    referenceModel.Points=referenceModel.Points(I,:);
    
    referenceModel.Triangles =delaunay(referenceModel.Points(:,1),referenceModel.Points(:,2));
    referenceModel.xGrid = length(unique(referenceModel.Points(:,1)));
    referenceModel.yGrid = length(unique(referenceModel.Points(:,2)));
    save([paths.referenceModelPath,'.mat'],'referenceModel');
end

% Calculate Reference Model Grid Size
P_Referenz=referenceModel.Points;
Xmax=max(P_Referenz(:,1));
Ymax=max(P_Referenz(:,2));
Zmax=max(P_Referenz(:,3));
Xmin=min(P_Referenz(:,1));
Ymin=min(P_Referenz(:,2));
Zmin=min(P_Referenz(:,3));
X_List_unique=unique(P_Referenz(:,1));
XstepSize=X_List_unique(2)-X_List_unique(1);
Y_List_unique=unique(P_Referenz(:,2));
YstepSize=Y_List_unique(2)-Y_List_unique(1);

Modellgrenzen(1)=Xmax;
Modellgrenzen(2)=Xmin;
Modellgrenzen(3)=XstepSize;
Modellgrenzen(4)=Ymax;
Modellgrenzen(5)=Ymin;
Modellgrenzen(6)=YstepSize;
Modellgrenzen(7)=Zmax;
Modellgrenzen(8)=Zmin;

% - parse initial PWL model from XML file ---------------------------------
[initialPWLModel.Points, ~, initialPWLModel.Triangles,variantName]=read_xml(paths.initialModelPath);

% % - determine grid size of initial PWL model (e.g. 6x5) -------------------
% Xpoint=unique(initialPWLModel.Points(:,1));
% Ypoint=unique(initialPWLModel.Points(:,2));
% isgrid=1;
% for i = 1:length(Xpoint)-2
%    if Xpoint(2)-Xpoint(1) ~= Xpoint(i+2)-Xpoint(i+1)
%        isgrid = 0;
%        break
%    end
% end
% for i = 1:length(Ypoint)-2
%    if Ypoint(2)-Ypoint(1) ~= Ypoint(i+2)-Ypoint(i+1)
%        isgrid = 0;
%        break
%    end
% end
% if isgrid == 1
%     initialPWLModel.xGrid = length(Xpoint);
%     initialPWLModel.yGrid = length(Ypoint);
% else
%     initialPWLModel.xGrid = 0;
%     initialPWLModel.xGrid = 0;
% end

% - indexing points -------------------------------------------------------
% 0 = inside ;1 = Xedge; 2 = Yedge; 3 = corner
initialPWLModel.Index = or(initialPWLModel.Points(:,1)==Xmax,initialPWLModel.Points(:,1)==Xmin);
initialPWLModel.Index = initialPWLModel.Index + 2*or(initialPWLModel.Points(:,2)==Ymax,initialPWLModel.Points(:,2)==Ymin);
%add revert counter
initialPWLModel.RevertCount=0;
% - assume initial PWL model ----------------------------------------------
PWLModel=initialPWLModel;
else
   error('Device Type nonexistent! In device folder please create an empty file with the name "Is___", where ___ is "Diode" or "Transistor".');

if (options.plot)
    
    scrsz=get(groot,'ScreenSize');

    fh=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    plotpoints(PWLModel.Points, PWLModel.Triangles,1)
    grid on;
    title('Initial PWL model');
    
    % save figure as png file
    saveas(fh, fullfile(paths.resultPath,'initialModel.png'), 'png')
    saveas(fh, fullfile(paths.resultPath,'initialModel.eps'), 'epsc2')
    clear fh;
end
if (strcmp(Typ,'MOD'))
    % - determine curvature of reference model --------------------------------
    [xz_grid, yz_grid] = curvature(referenceModel.Points, Modellgrenzen, options);

    % - create interpolation grid of reference model --------------------------
    Interpolated_Ref = scatteredInterpolant(P_Referenz(:,1),P_Referenz(:,2),P_Referenz(:,3));
else
  xz_grid = [];
  yz_grid = [];
  Interpolated_Ref = [];
end

% - start SA run for determine start temperature --------------------------
if (options.debugMode)
    disp('INFO: Determine Start Temperature T0 ...');
    tic;
end
[T0, currentPWLModelPoints, Var_Bewertung, Var_Bewertung_Max, newRevertCount]=Starttemperatur(PWLModel,P_Referenz,Typ,options,Modellgrenzen,chi,relativeSchrittweite, xz_grid, yz_grid, Bias, Interpolated_Ref);
if (options.debugMode)
    toc;
    disp('... DONE: T0 determined.')
end
% - assume optimized PWL model from previous SA run -----------------------
PWLModel.Points=currentPWLModelPoints;
PWLModel.RevertCount=newRevertCount;
%update Triangulation
PWLModel.Triangles=delaunay(PWLModel.Points(:,1),PWLModel.Points(:,2));

if (options.plot)
    plotpoints(PWLModel.Points,PWLModel.Triangles,5)
    grid on;
    title('1. PWL model (after start temperature)');
end

[Var_Bewertung_Initial, Var_Bewertung_Initial_Max]=error_value(initialPWLModel.Points, P_Referenz);
Var_Bewertung_T0=Var_Bewertung;
Var_Bewertung_T0_Max=Var_Bewertung_Max;
Var_Bewertung_Best=Var_Bewertung;
Var_Bewertung_Best_Max=Var_Bewertung_Max;
% Var_Bewertung_Worst=Var_Bewertung;
% P_Init=PWLModel.Points;
P_Best=PWLModel.Points;
% P_Worst=PWLModel.Points;


% Start SA

% initialize some variables
Vektor_T=T0;
Vektor_C=Var_Bewertung_T0;
Vektor_C_Max=Var_Bewertung_T0_Max;
Vektor_C_accept=Var_Bewertung_T0;
Vektor_C_accept_Max=Var_Bewertung_T0_Max;
Vektor_C_best=Var_Bewertung_T0;
Vektor_C_best_Max=Var_Bewertung_T0_Max;
Vektor_dC=0;
xSteps=[0];

if (options.debugMode)
    disp('INFO: Execute Simulated Annealing ...');
    tic;
end
%else
   error('Device Type nonexistent! In device folder please create an empty file with the name "Is___", where ___ is "Diode" or "Transistor".'); open waitbar for progress status
progressBar = waitbar(0,'','Name','Optimize PWL model...');

for i=1:T_Steps
    for j=1:T_Hold
        xSteps=[xSteps (((i-1)*T_Hold)+j)];
        % determine progress status
        progressStatus=((((i-1)*T_Hold)+j)/(T_Steps*T_Hold))*100;
        % update waitbar progress status
        waitbar(progressStatus/100,progressBar,[num2str((progressStatus)) '%']);
        
        % random selection of vertex index
        k=randi(length(PWLModel.Points));

        % save old vertices of PWL model
        oldPWLModel.Points=PWLModel.Points;
        
        % find probability weighting
        if (strcmp(Typ,'MOD'))
            weight = move_Weight(xz_grid, yz_grid, PWLModel.Points, Bias, k, Modellgrenzen);
        else
            weight=[];
        end

        % move selected vertex of PWL model
        [vertex, newRevertCount] = move_Point(PWLModel, relativeSchrittweite, weight, k, Interpolated_Ref, Typ, Modellgrenzen);
        PWLModel.Points(k,1)=vertex(1);
        PWLModel.Points(k,2)=vertex(2);
        PWLModel.Points(k,3)=vertex(3);
        PWLModel.RevertCount=newRevertCount;
        assert(isempty(find(isnan(PWLModel.Points) == 1)))

        % evaluate changed PWL model
        [Var_Bewertung_Current_Mean, Var_Bewertung_Current_Max]=error_value(PWLModel.Points, P_Referenz);

        %Bewertung
        deltaCost=Var_Bewertung_Current_Mean-Vektor_C_accept(end);

        if (options.debugMode)
            Vektor_T=[Vektor_T Funktion(T0,i,options.Temperaturfaktor)];
            Vektor_C=[Vektor_C Var_Bewertung_Current_Mean];
            Vektor_C_Max=[Vektor_C_Max Var_Bewertung_Current_Max];
            Vektor_dC=[Vektor_dC deltaCost];
        end
        
        if deltaCost<0; % Improvement of PWL model
            Vektor_C_accept=[Vektor_C_accept Var_Bewertung_Current_Mean];
            % Check for best model
            if Var_Bewertung_Current_Mean<Var_Bewertung_Best
                Var_Bewertung_Best=Var_Bewertung_Current_Mean;
                Var_Bewertung_Best_Max=Var_Bewertung_Current_Max;
                Vektor_C_best=[Vektor_C_best Var_Bewertung_Best];
                Vektor_C_best_Max=[Vektor_C_best_Max Var_Bewertung_Best_Max];
                P_Best=PWLModel.Points;
            else
                Vektor_C_best=[Vektor_C_best Vektor_C_best(end)];
                Vektor_C_best_Max=[Vektor_C_best_Max Vektor_C_best_Max(end-2:end)];
            end
        else % deltaCost >= 0: Deterioration of PWL model
            % some Debug information
%             if (options.debugMode)
                Vektor_C_best=[Vektor_C_best Vektor_C_best(end)];
                Vektor_C_best_Max=[Vektor_C_best_Max Vektor_C_best_Max(end-2:end)];
%             end
%             if Var_Bewertung_Current>Var_Bewertung_Worst
%                 Var_Bewertung_Worst=Var_Bewertung_Current;
%                 P_Worst=PWLModel.Points;
%                 T_Worst=T;
%                 P_Worst_fine=P_prev_fine;
%             end

            %f(i) Temperaturfunktion
            %Boltzmann distribution: p=exp(-dE/kT) with k=1
            p=exp(-deltaCost/Funktion(T0,i,options.Temperaturfaktor));
            if rand(1)<=p;
                Vektor_C_accept=[Vektor_C_accept Var_Bewertung_Current_Mean];
                Vektor_C_accept_Max=[Vektor_C_accept_Max Var_Bewertung_Current_Max];
            else
                % revert vertex move, restore old PWL model
                PWLModel.Points=oldPWLModel.Points;
                Vektor_C_accept=[Vektor_C_accept Vektor_C_accept(end)];
                Vektor_C_accept_Max=[Vektor_C_accept_Max Vektor_C_accept_Max(end-2:end)];
            end
        end
        assert( length(xSteps) == length(Vektor_C_best));
    end
end

% update PWL model
PWLModel.Points=P_Best;
%update Triangulation
PWLModel.Triangles=delaunay(PWLModel.Points(:,1),PWLModel.Points(:,2));
% count number move reverted
PercentReverted=PWLModel.RevertCount/(options.T_Start_Steps+T_Steps*T_Hold+PWLModel.RevertCount)*100;
% calculate percentage improvement
Improvement=(Var_Bewertung_Initial-Vektor_C_best(end))/Var_Bewertung_Initial*100;

if (options.plot)
    if (strcmp(Typ,'MOD'))
        % add points trajectory
        figure(2)
        hold on;
        % t = plot3(PWLModel.Points(:,1),PWLModel.Points(:,2),20.*ones(size(PWLModel.Points,1),1));
        % set(t,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',8);
        % hold off;
        for j = 1:size(PWLModel.Points,1)
            t =plot3([initialPWLModel.Points(j,1),PWLModel.Points(j,1)],[initialPWLModel.Points(j,2),PWLModel.Points(j,2)],[20,20]);
            set(t,'LineStyle','-','LineWidth',3,'Color','k','Marker','.','MarkerSize',20);
            text(initialPWLModel.Points(j,1),initialPWLModel.Points(j,2),20,int2str(j),'HorizontalAlignment','left','Color','k','FontSize',15,'FontWeight','bold');
        end
        hold off;
        % save figure as png file
        fh = figure(2);
        saveas(fh, fullfile(paths.resultPath,'curvature_x-z_plane.png'), 'png')
        saveas(fh, fullfile(paths.resultPath,'curvature_x-z_plane.eps'), 'eps2c')
        clear fh;

        % add points trajectory
        figure(3)
        hold on;
        % t = plot3(PWLModel.Points(:,1),PWLModel.Points(:,2),20.*ones(size(PWLModel.Points,1),1));
        % set(t,'LineStyle','none','Marker','o','MarkerFaceColor','k','MarkerSize',8);
        % hold off;
        for j = 1:size(PWLModel.Points,1)
            t =plot3([initialPWLModel.Points(j,1),PWLModel.Points(j,1)],[initialPWLModel.Points(j,2),PWLModel.Points(j,2)],[20,20]);
            set(t,'LineStyle','-','LineWidth',3,'Color','k','Marker','.','MarkerSize',20);
            text(initialPWLModel.Points(j,1),initialPWLModel.Points(j,2),20,int2str(j),'HorizontalAlignment','left','Color','k','FontSize',15,'FontWeight','bold');
        end
        hold off;
        % save figure as png file
        fh = figure(2);
        saveas(fh, fullfile(paths.resultPath,'curvature_y-z_plane.png'), 'png')
        saveas(fh, fullfile(paths.resultPath,'curvature_y-z_plane.eps'), 'eps2c')
        clear fh;
    end
end


if (options.debugMode)
    toc;
    disp('... DONE: SA');

    PWLCharacteristics.Vektor_C                 = Vektor_C;
    PWLCharacteristics.Vektor_C                 = Vektor_C_Max;
    PWLCharacteristics.Vektor_C_accept          = Vektor_C_accept;
    PWLCharacteristics.Vektor_C_accept_Max      = Vektor_C_accept_Max;
    PWLCharacteristics.Vektor_C_best            = Vektor_C_best;
    PWLCharacteristics.Vektor_C_best_Max        = Vektor_C_best_Max;
    PWLCharacteristics.Vektor_dC                = Vektor_dC;
    PWLCharacteristics.Vektor_T                 = Vektor_T;
    PWLCharacteristics.Var_Bewertung_Best       = Var_Bewertung_Best;
    PWLCharacteristics.Var_Bewertung_Best_Max	= Var_Bewertung_Best_Max;
%     PWLCharacteristics.Var_Bewertung_Worst    = Var_Bewertung_Worst;
%     PWLCharacteristics.T_Best                 = T_Best;
%     PWLCharacteristics.T_Worst                = T_Worst;
else
    PWLCharacteristics.Vektor_C_best            = Vektor_C_best;
    PWLCharacteristics.Vektor_C_best_Max    	= Vektor_C_best_Max;
    PWLCharacteristics.PercentReverted          = PercentReverted;
    PWLCharacteristics.Improvement              = Improvement;
end

% write PWL model characteristics
fid = fopen(fullfile(paths.resultPath,'metric_pwl.csv'), 'w');
fprintf(fid, '%s,', 'MSE(initial):');
fprintf(fid, '%f\n',Var_Bewertung_Initial);
fprintf(fid, '%s,', 'MaxSE(initial):');
fprintf(fid, '%f\n',Var_Bewertung_Initial_Max(3));
fprintf(fid, '%s,', 'MaxSE(initial)X.Y:');
fprintf(fid, '%f;%f\n',Var_Bewertung_Initial_Max(1),Var_Bewertung_Initial_Max(2));
fprintf(fid, '%s,', 'MSE(T0):');
fprintf(fid, '%f\n',PWLCharacteristics.Vektor_C_best(1));
fprintf(fid, '%s,', 'MaxSE(T0):');
fprintf(fid, '%f\n',PWLCharacteristics.Vektor_C_best_Max(3));
fprintf(fid, '%s,', 'MaxSE(T0)X.Y:');
fprintf(fid, '%f;%f\n',PWLCharacteristics.Vektor_C_best_Max(1),PWLCharacteristics.Vektor_C_best_Max(2));
fprintf(fid, '%s,', 'MSE(best):');
fprintf(fid, '%f\n',PWLCharacteristics.Vektor_C_best(end));
fprintf(fid, '%s,', 'MaxSE(best):');
fprintf(fid, '%f\n',PWLCharacteristics.Vektor_C_best_Max(end));
fprintf(fid, '%s,', 'MaxSE(best)X.Y:');
fprintf(fid, '%f;%f\n',PWLCharacteristics.Vektor_C_best_Max(end-2),PWLCharacteristics.Vektor_C_best_Max(end-1));
% fprintf(fid, '%s,', 'X-Grid:');
% fprintf(fid, '%d\n',PWLModel.xGrid-1);
% fprintf(fid, '%s,', 'Y-Grid:');
% fprintf(fid, '%d\n',PWLModel.yGrid-1);
fprintf(fid, '%s,', 'Number of Move reverted(%):');
fprintf(fid, '%d\n',PercentReverted);
fprintf(fid, '%s,', 'Improvement(%):');
fprintf(fid, '%d\n',Improvement);
fclose(fid);

% For DEBUG
if (options.plot)
    fh=figure(6);
    set(fh,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    plotpoints(PWLModel.Points,PWLModel.Triangles,6)
    grid on;
    title('Optimized PWL model');

    % save figure as png file
    saveas(fh, fullfile(paths.resultPath,'optimizedModel.png'), 'png')
    saveas(fh, fullfile(paths.resultPath,'optimizedModel.eps'), 'eps2c')
    clear fh;
    
    plotpoints(PWLModel.Points,PWLModel.Triangles,7)
    hold on;
    plotpoints(referenceModel.Points,referenceModel.Triangles,7)
    title('Optimized PWL model with reference PWL model');
    hold off;

end

if (options.debugMode)
    fh=figure(8);
    set(fh,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
    subplot(4,1,1);
    plot(xSteps, Vektor_C_accept);
    grid on;
    title('accepted cost profile');
    xlabel('N');
    ylabel('accepted cost');

    subplot(4,1,2);
    plot(xSteps, Vektor_dC);
    range =sqrt(sum(Vektor_dC.^2)/size(Vektor_dC,2));
    grid on;
    title('delta cost profile');
    axis([-inf inf -(range/2) range]);
    xlabel('N');
    ylabel('delta cost');

    subplot(4,1,3);
    plot(xSteps, Vektor_C_best);
    grid on;
    title('best cost profile');
    xlabel('N');
    ylabel('best cost');
    
    subplot(4,1,4);
    plot(xSteps, Vektor_T);
    grid on;
    title('temperature profile');
    xlabel('N');
    ylabel('temperature');
    
    % save figure as eps file
    saveas(fh, fullfile(paths.resultPath,'metric_pwl.fig'), 'fig')
    saveas(fh, fullfile(paths.resultPath,'metric_pwl.eps'), 'eps2c')
    clear fh;
end

try
    % close waitbar for progress status
    close(progressBar);
end
end
