%Dioden Modell

clear;

if exist(fullfile('data', 'settings.m'), 'file')
	run(fullfile('data', 'settings.m'));
end

if ~(exist('element_name', 'var'))
	error('Variable "element_name" is not defined!')
end

if ~(ischar(element_name))
    error('Variable "element_name" is not a string!')
end

if ~(exist('sigma', 'var'))
        error('Variable "sigma" is not defined!')
end
if ~(isfloat(sigma))
    error('Variable "sigma" is not a number!')
end

data = dlmread(fullfile('data', 'saber_output.txt'), '\t', 1, 0);
Vd = data(:, 1);
Id = data(:, 2);

plot(Vd, Id);



% %set number of breakpoints
% %sigma=...;
% %sigma=6
% extra_breakpoint=1;
% sorted=0;
% 
% %%generate pwl-model
% %[Xp,Yp,Z2,P,Ypwl,sigma] = pwl_optimal_diode(Vd,Id,sigma,extra_breakpoint,sorted);
% %
% %%parameter
% %P
% 
% S = sigma;
% 
% [er, a, b, xb] = optimalLinearApproximation(Vd, Id, S, 0);
% 
% for i = 1:S
%     Xp(1,i) = xb(S,i);
%     Yp(1,i) = Xp(1,i) * a(1,S,i) + b(1,S,i);
% end
% 
% Xp(1,S+1) = xb(S,S+1);
% Yp(1,S+1) = Xp(1,S+1) * a(1,S,S) + b(1,S,S);
% 
% %plot
% plot(Vd,Id,Xp,Yp,'-*')

%%create outputfile
%fid = fopen('model.dat', 'w');
%    fprintf(fid, 'MODULENAME:\n');
%    fprintf(fid, '%s\n', element_name);
%    fprintf(fid, 'NODES:\n');
%    fprintf(fid, 'p n\n');
%    fprintf(fid, 'ELEMENTS:2\n');
%    fprintf(fid, 'r.ID p:p m:n = rnom=ID\n');
%    fprintf(fid, 'i.ID p:p m:n = dc=ID, ac=ID\n');
%    fprintf(fid, 'INPUT1:\n');
%    fprintf(fid, 'voltage p n %d [',sigma+1);
%    fprintf(fid, '%g ', Vd(1,1));
%    for i=1:sigma
%        fprintf(fid, '%g ', Z2(i));
%    end
%    fprintf(fid, '%g', Vd(l_d));
%    fprintf(fid, ']\n');
%    fprintf(fid, 'PWLMODEL:%d\n', sigma+1);
%    fprintf(fid, '[');
%    for i=1:sigma+1
%%      fprintf(fid, '(%g %g)', P(i,1), P(i,2)); 
%       fprintf(fid, '(%g %g)', inv(P(i,1)), P(i,2)); 
%       if i<sigma+1
%       fprintf(fid,',');
%       end      
%    end
%    fprintf(fid, ']\n');
%fclose(fid);

fxml = fopen('model.xml', 'w');
fprintf(fxml, '<praise>\n');
fprintf(fxml, '<pwlmodels>\n');
fprintf(fxml, '<element name="%s" variant="VARIANT_NAME">\n', element_name);
fprintf(fxml, '\t<model1dim name="id">\n');
fprintf(fxml, '\t\t<input type="voltage" node1="p" node2="n"/>\n');
fprintf(fxml, '\t\t<output type="current" node1="p" node2="n"/>\n');
fprintf(fxml, '\t\t<datapoints>\n');
for i=1:size(Vd,1)
    fprintf(fxml, '\t\t\t<datapoint x="%g" y="%g"/>\n', Vd(i,1), Id(i,1));
end
fprintf(fxml, '\t\t</datapoints>\n');
fprintf(fxml, '\t</model1dim>\n');
fprintf(fxml, '</element>\n');
fprintf(fxml, '</pwlmodels>\n');
fprintf(fxml, '</praise>\n');
fclose(fxml);

dlmwrite('model.csv', [Vd, Id]);

save as .png image 
 print -dpng Dioden-Kennlinie.png;
