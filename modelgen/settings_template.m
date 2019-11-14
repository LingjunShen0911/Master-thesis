% parse reference PWL model 
xmlName = '###referenceModel###';
filename_referenceModel= ['models/', xmlName, '.xml'];
disp(['INFO: Reference model: ' filename_referenceModel]);
filename_initialModel='###initialModel###';
disp(['INFO: Initial model: ' filename_initialModel]);

if exist(['models/', xmlName, '.mat']) == 2
    load(['models/', xmlName, '.mat']);
else
    referenceModel.Points=read_xml(filename_referenceModel);
    
    %Sort data order(important for Modified SA)
    [~,I]=sort(referenceModel.Points(:,2));
    referenceModel.Points=referenceModel.Points(I,:);
    [~,I]=sort(referenceModel.Points(:,1));
    referenceModel.Points=referenceModel.Points(I,:);
    
    referenceModel.Triangles =delaunay(referenceModel.Points(:,1),referenceModel.Points(:,2));
    referenceModel.xGrid = length(unique(referenceModel.Points(:,1)));
    referenceModel.yGrid = length(unique(referenceModel.Points(:,2)));
    save(['models/', xmlName, '.mat'],'referenceModel');
end

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

options.debugMode=true;
options.plot=true;

%function must be >=0
Funktion=@(T0,x,Temperaturfaktor) T0*10.^(x*log(Temperaturfaktor));
options.Temperaturfaktor=0.8;
chi=0.8;
Bias = 0.3;

%%% variables set by run script %%% 
% relativeSchrittweite.X = 1;
% relativeSchrittweite.Y = 1;
% relativeSchrittweite.Z = 1;
% options.T_Start_Steps = 2;
% T_Steps = 2;
% T_Hold = 2;
% Typ = 'MOD';  
