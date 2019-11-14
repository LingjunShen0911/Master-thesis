function generate_diode_xml(elementName, variantName, paths)

%Parse simulation Data
P = dlmread(paths.referenceSimDataPath,'	',1,1);

%Calculate segment with least error
nseg=str2num(variantName);
K=nseg-1;

[ipt,Error] = findchangepts(P(:,2),'MaxNumChanges',K,'Statistic','linear');

y = P(:,2)';
icp = ipt';

m = size(y,1);
n = size(y,2);

% algorithm taken from function cpplot (part of findchangepts) ----------------
% create an approximated line that fit the curve best on each segment
for r=1:m
    xData = [1 icp-0.5; icp-0.5 n; nan(1, nseg)];
    yData = NaN(3,nseg);
    istart = [1 icp];
    istop = [icp-1 n];

    for s=1:nseg
        ix = (istart(s):istop(s));
        b1 = 0.5*(s>1);
        b2 = 0.5*(s<nseg);
        yData(1:2,s) = polyval(polyfit(ix,y(r,ix),1),[ix(1)-b1 ix(end)+b2]);
    end
end
xData = reshape(xData, 3*nseg, 1);
yData = reshape(yData, 3*nseg, 1);
% -----------------------------------------------------------------------------

%convert line Data into points Data
yPoint(1) = yData(1);
xPoint(1) = P(1,1);
for i = 1:K
    yPoint(i+1) = (yData(2+3*(i-1))+yData(4+3*(i-1)))/2;
    xPoint(i+1) = P(ipt(i),1);
end
yPoint(nseg+1) = yData(end-1);
xPoint(nseg+1) = P(end,1);

%Make .xml
fxml = fopen(paths.resultPath, 'w');
	fprintf(fxml, '<praise>\n');
	fprintf(fxml, '<pwlmodels>\n');
	fprintf(fxml, ['<element name="%s" variant="' variantName '">\n'], elementName);
	fprintf(fxml, '\t<model1dim name="id">\n');
	fprintf(fxml, '\t\t<input type="voltage" node1="p" node2="n"/>\n');
	fprintf(fxml, '\t\t<output type="current" node1="p" node2="n"/>\n');
	fprintf(fxml, '\t\t<datapoints>\n');
	for i=1:size(xPoint,2)
		x = xPoint(1,i);
		y = yPoint(1,i);
		fprintf(fxml, '\t\t\t<datapoint x="%g" y="%g"/>\n', x, y);
	end;
	fprintf(fxml, '\t\t</datapoints>\n');
	fprintf(fxml, '\t</model1dim>\n');
    fprintf(fxml, '</element>\n');
    fprintf(fxml, '</pwlmodels>\n');
	fprintf(fxml, '</praise>\n');
fclose(fxml);