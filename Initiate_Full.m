clear, close all, clc;

load 5Cases.mat

GeoVariableOpt = GeoVariableOpt';
VinMaxOpt_L = VinMaxOpt;
VoutMaxOpt_L = VoutMaxOpt;
VinMaxOpt_L(5) = 15;
VoutMaxOpt_L(5) = 15;

for geonumber = 2:2 % Which case to go to

GeoVariables = GeoVariableOpt(geonumber,:);
VinMaxOpt = VinMaxOpt_L(geonumber);
VoutMaxOpt = VoutMaxOpt_L(geonumber);
    
[Inequality, Equality] = GeoConstraint(GeoVariables);

[OIB,IIB, OOB,IOB, ORBr,ORBl,IRB, OHLt,OHLb,OHLl,OHLr,IHL, OMLtr,OMLtl,OMLbr,OMLbl,IML, Case] = All_STL_Gen( GeoVariables );

run('Thermal.m');
run('Mechanical.m');

Final(1,:) = Disp_Out(1);
Final(2,:) = Disp_Out(2);
Final(3,1) = Stroke;
Final(4,1) = Mech_stiffness;
Final(5,1) = F_Output;
end

figure

scatter(Coordinates(:,1),Coordinates(:,2));
a = [1:num_nodes]'; b = num2str(a); c = cellstr(b);
dx = 0.00001; dy = 0.00001; % displacement so the text does not overlay the data points
text(Coordinates(:,1)+dx, Coordinates(:,2)+dy, c);

hold on

x1 = zeros(num_links,1);
x2 = zeros(num_links,1);
y1 = zeros(num_links,1);
y2 = zeros(num_links,1);
for i = 1:num_links
    x1(i,1) = Coordinates(Links(i,1),1);
    x2(i,1) = Coordinates(Links(i,1),2);
    y1(i,1) = Coordinates(Links(i,2),1);
    y2(i,1) = Coordinates(Links(i,2),2);
end
x = [x1 y1];
y = [x2 y2];
plot(x',y')