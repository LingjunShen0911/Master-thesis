% Settings
options.debugMode=true;
options.plot=true;

%Temperature Function
%function must be >=0
Funktion=@(T0,x,Temperaturfaktor) T0*10.^(x*log(Temperaturfaktor));
options.Temperaturfaktor=0.8;
chi=0.8;
Bias = 0.3;

%Variables
relativeSchrittweite.X = .01666;
relativeSchrittweite.Y = .01666;
relativeSchrittweite.Z = 1e-4;
options.T_Start_Steps = 1200;
T_Steps = 30;
T_Hold = 100;

%" available modes: 'NORM' and 'MOD'
Typ = 'NORM';
