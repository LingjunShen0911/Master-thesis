% Generate Constructive Model
%   TODO
%
% Outputs :
%   -
%
% Example:
%   TODO
%
% TODO List :
%
% Dependency : none
%
% Source : Dynamic Programming, Douglas-something
% 
% TODO :
%   > Exit(1) when both signal have different number of dimension
%
% Documentation :
%   TODO
% 
% Author: Yeremia Gunawan A. (yeremia.gunawan@ims.uni-hannover.de)
% Written: 24-Oktober-2016
% Last update: 26-Oktober-2016
% Version : 1.0 
%% Create new structure
filename_initialModel='models/irff130_ref_57800_t.xml';
if exist('models/reference.mat') == 2
    load('models/reference.mat');
else
    [referencePWLModel.Points, ~, referencePWLModel.Triangles]=read_xml(filename_initialModel);
    referencePWLModel.xGrid = length(unique(referencePWLModel.Points(:,1)));
    referencePWLModel.yGrid = length(unique(referencePWLModel.Points(:,2)));
    save('models/reference.mat','referencePWLModel');
end
PWLModel=referencePWLModel;

%% Initialize
RefPts = PWLModel.Points;
Corner = and( ...
        or(RefPts(:,1) == min(RefPts(:,1)),  RefPts(:,1) == max(RefPts(:,1))) ...
        ,or(RefPts(:,2) == min(RefPts(:,2)),  RefPts(:,2) == max(RefPts(:,2))) ...
        );

PtsTaken = RefPts(Corner,:);
PtsAvail = RefPts(~Corner,:);
numOfTri = 0;

F = scatteredInterpolant(PtsTaken(:,1),PtsTaken(:,2),PtsTaken(:,3));
zInterpolated = F(RefPts(:,1),RefPts(:,2));
squaredError = (zInterpolated-RefPts(:,3)).^2;
oldMSE = 1/size(zInterpolated,1)*sum(squaredError);
saturationTreshold = 0.08;
%%
% numOfTri = 0;
while numOfTri < 200 
    MSEs = zeros(size(PtsAvail,1),1);
    squaredErrors = zeros(size(PtsAvail,1),1);
%     MSEs = zeros(size(PtsAvail,1),1);
    squaredErrors = zeros(size(PtsAvail,1),1);
%     tic
%     h = waitbar(0);
    for i=1:size(PtsAvail,1)
        PtsTakenWithNewPts = cat(1,PtsTaken,PtsAvail(i,:));
        F = scatteredInterpolant(PtsTakenWithNewPts(:,1),PtsTakenWithNewPts(:,2),PtsTakenWithNewPts(:,3));
        zInterpolated = F(RefPts(:,1),RefPts(:,2));
        squaredError = (zInterpolated-RefPts(:,3)).^2;
        % First column of errorMatrix is MSE and second column is maximum
        % squared error
        MSEs(i) = 1/size(zInterpolated,1)*sum(squaredError);
        squaredErrors(i) = max(squaredError);
%         waitbar(i/size(PtsAvail,1), h, [num2str(i),'/',num2str(size(PtsAvail,1))]);
%         waitbar(i/size(PtsAvail,1), h, sum(MSEs(:) ~= 0));
    end
%     close(h)
    toc
    [ minimumNewMSE, ind ] =  min(MSEs);
    if (oldMSE - minimumNewMSE)/oldMSE < saturationTreshold
        [maxSquaredError , ind] = max(squaredErrors);
    end
    PtsTaken = cat(1,PtsTaken, PtsAvail(ind,:));
    PtsAvail(ind,:) = []; 
    tri = delaunay(PtsTaken(:,1),PtsTaken(:,2));
    numOfTri = size(tri,1);
    oldMSE = minimumNewMSE;
    disp(['MSE : ',num2str(oldMSE), ', Max Squared Error : ', num2str(maxError)]);    
    save(['vars\','tri_', num2str(numOfTri)]);
end