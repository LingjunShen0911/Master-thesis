%% Simple op-amp curve creation
P = dlmread("examples/lt1363/LT1363.csv",'\t',1,0);

ind = findchangepts(P(:,2),'MaxNumChanges',2,'Statistic','linear');

%% Inputvalues
u = [P(1,1) P(ind(1),1) 0 P(ind(2),1) P(size(P,1),1)];

%% Outputvalues
i = [P(1,2) P(ind(1),2) 0 P(ind(2),2) P(size(P,1),2)];

disp("Input values");
disp(u);
disp("Output values");
disp(i);

u(1) = -3;
u(5) = 3;

Points = [];
Tris = [];
disp(Points);

%                 P1                     P2                P3
Points = [Points; [u(1) u(1)-u(2) i(2)]; [u(1) -u(1) i(1)]; [-u(1)+u(2) -u(1) i(2)]];
%                 Pa                    P5                     P6
Points = [Points; [-u(5)+u(4) -u(5) i(4)];[u(5) u(5)-u(4) i(4)]; [u(5) -u(5) i(5)]];
Tris = [[1 2 3]; [1 3 4]; [3 4 5]; [4 5 6]];
scrsz=get(groot,'ScreenSize');

fh=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
plotpoints(Points, Tris,1)
grid on;
title('Initial PWL model');

clear fh;