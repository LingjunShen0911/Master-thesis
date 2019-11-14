function [T0, P_Best, Var_Bewertung_Best, Var_Bewertung_Best_Max, newRevertCount]=Starttemperatur(PWLModel,P_Referenz,Typ,options,Modellgrenzen,chi,relativeSchrittweite, xz_grid, yz_grid, Bias, Interpolated_Ref)
% STARTTEMPERTATUR
% 
% Quelle: Modellierung komplexer Prozesse durch naturanaloge Verfahren
% Autoren: Kl�ver, Kl�ver, Schmidt
% Erschienen: 2012 (2. Auflage)
% 
%   Inputs:
%       TODO
%   Output: Starttemperatur für den Haupt SA-Lauf
%       TODO
%
% Arbeitsweise: Das SA wird eine angegebene Anzahl oft ausgeführt.
% Verschlechterungen werden aufsummiert. Der "Mittelwert aller
% Verschlechterungen / log(chi) " bildet die Starttemperatur.
%
% Authors: Axel Mensing, Lukas Lee (Lukas.Lee@ims.uni-hannover.de), Gerald Koroa
% Last Edit: 28.10.16

    if (options.debugMode)
        Vektor_C=[];
        Vektor_C_Max=[];
        Vektor_dC=[];
    end

    % initialize some variables
    T_Start_Steps=options.T_Start_Steps;
    dCSum=0;
    NoOfDeterioration=0;
    xSteps=[0];
    
    % evaluate intial PWL model
    P_Best=PWLModel.Points;
    [Var_Bewertung_New_Mean, Var_Bewertung_New_Max]=error_value(PWLModel.Points, P_Referenz);
    Var_Bewertung_Current=Var_Bewertung_New_Mean;
    Var_Bewertung_Current_Max=Var_Bewertung_New_Max;
    Var_Bewertung_Best=Var_Bewertung_Current;
    Var_Bewertung_Best_Max=Var_Bewertung_Current_Max;
 
    if (options.debugMode)
        Vektor_C=[Vektor_C Var_Bewertung_Current];
        Vektor_C_Max=[Vektor_C_Max Var_Bewertung_Current_Max];
        Vektor_dC=[Vektor_dC 0];
    end

    % start SA for compute T0
    
    % open waitbar for progress status
    progressBar = waitbar(0,'','Name','Compute start temperature ...');

    for i=1:T_Start_Steps
        xSteps=[xSteps i];
        % determine progress status
        progressStatus=i/T_Start_Steps*100;
        % update waitbar progress status
        waitbar(progressStatus/100,progressBar,[num2str((progressStatus)) '%']);
        
        % random selection of vertex index
        k=randi(length(PWLModel.Points));
        
        % find probability weighting
        if (strcmp(Typ,'MOD'))
            weight = move_Weight(xz_grid, yz_grid, PWLModel.Points, Bias, k, Modellgrenzen);
        else
            weight=[];
        end
        
        % move selected vertex of PWL model
        [vertex, newRevertCount] = move_Point(PWLModel, relativeSchrittweite, weight, k, Interpolated_Ref, Typ, Modellgrenzen);

        % update PWL model
        PWLModel.Points(k,1)=vertex(1);
        PWLModel.Points(k,2)=vertex(2);
        PWLModel.Points(k,3)=vertex(3);
        PWLModel.RevertCount=newRevertCount;
        assert(isempty(find(isnan(PWLModel.Points) == 1)))

        % evaluate changed PWL model
        [Var_Bewertung_New_Mean,Var_Bewertung_New_Max] =error_value(PWLModel.Points, P_Referenz);

        %evaluate cost
        deltaCost=Var_Bewertung_New_Mean-Var_Bewertung_Current;

        if (options.debugMode)
            Vektor_dC=[Vektor_dC deltaCost];   
            Vektor_C=[Vektor_C Var_Bewertung_New_Mean];
            Vektor_C_Max=[Vektor_C_Max Var_Bewertung_New_Max];
        end
        
        % IF: deltaCost < 0: Improvement of PWL model
        if (deltaCost < 0);
            Var_Bewertung_Current=Var_Bewertung_New_Mean;
            Var_Bewertung_Current_Max=Var_Bewertung_New_Max;

            % search for best PWL model
            if (abs(Var_Bewertung_New_Mean) < abs(Var_Bewertung_Best))
                P_Best=PWLModel.Points;
                Var_Bewertung_Best=Var_Bewertung_New_Mean;
                Var_Bewertung_Best_Max=Var_Bewertung_New_Max;
            end
            
        % ELSE: deltaCost > 0: Deterioration of PWL model
        elseif (deltaCost>=0)
            NoOfDeterioration=NoOfDeterioration+1;
            dCSum=dCSum+deltaCost;
            Var_Bewertung_Current=Var_Bewertung_New_Mean;
            Var_Bewertung_Current_Max=Var_Bewertung_New_Max;
        else
            assert(false, 'ERROR: deltaCost'); % Should never reached
        end
    end
    
    try
        % close waitbar for progress status
        close(progressBar);
    end

    % compute start temperature from dCSum and NoOfDeterioration
    T0=abs(dCSum/(NoOfDeterioration*log(chi)));

    if (options.debugMode)
        fh = figure;
        subplot(2,1,1); % top subplot
        plot(xSteps, Vektor_C);
        grid on;
        title('Start temp.: cost profile');
        xlabel('N');
        ylabel('cost');

        subplot(2,1,2); % bottom subplot
        plot(xSteps, Vektor_dC);
        range =sqrt(sum(Vektor_dC.^2)/size(Vektor_dC,2));
        grid on;
        title('Start temp.: delta cost profile');
        axis([-inf inf -(range/2) range]);
        xlabel('N');
        ylabel('delta cost');
        
        clear fh;

    end
    
end
