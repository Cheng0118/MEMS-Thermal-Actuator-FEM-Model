 %% Set Parameters %%
% Vin
Vin_E = VoutMaxOpt;
Vin_C = VinMaxOpt;

% Material Properties
k = 148;
k_O = 1.4;
row = 0.0002; % 0.005 - 0.020 ohm-cm

% Structure Parameters
TD = 0.4;
TO = 0.004;
TH = 0.05;
Gtb = 0.01;
Gcb = 0.01;
Gbb = 0.01;
Gmc = 0.05;

% Number of beams
N1 = GeoVariables(9); N2 = GeoVariables(1); N3 = GeoVariables(5);

%% Parameters to Input
% Number of nodes in beam
n_n = 5;

% Nodes
num_nodes_outer = 2+7*2;
num_nodes_inner = 4*2;

num_nodes_anchor = 5*2;

num_nodes = num_nodes_outer + num_nodes_inner + num_nodes_anchor;

% Links
num_links_outer = 2*(1+N1+2+N3+2)+1;
num_links_inner = N2+3*2;

num_links_anchor = 6*2;

num_links = num_links_outer+num_links_inner+num_links_anchor;

% Conductivity
Cond = zeros(num_links,1);

%% Connections %%
Links = zeros(num_links,2);

%% Outer
Links(1,1) = 1;
Links(1,2) = 2;
for i = 1:N1
    Links(1+i,1) = 2;
    Links(1+i,2) = 3;
end
for i = 1:2
    Links(1+N1+i,1) = 3+(i-1);
    Links(1+N1+i,2) = 4+(i-1);
end
for i = 1:N3
    Links(1+N1+2+i,1) = 5;
    Links(1+N1+2+i,2) = 6;
end
for i = 1:2
    Links(1+N1+2+N3+i,1) = 6+(i-1);
    Links(1+N1+2+N3+i,2) = 7+(i-1);
end
%--------------------------------------------------------------------------
Links(1+N1+2+N3+2+1,1) = 1;
Links(1+N1+2+N3+2+1,2) = 9;
for i = 1:N1+2+N3+2
    Links(1+N1+2+N3+2+1+i,1) = Links(1+i,1)+7;
    Links(1+N1+2+N3+2+1+i,2) = Links(1+i,2)+7;
end
%--------------------------------------------------------------------------
Links(num_links_outer,1) = 1;
Links(num_links_outer,2) = 16;

%% Inner
for i = 1:3
    Links(num_links_outer+i,1) = num_nodes_outer+i;
    Links(num_links_outer+i,2) = num_nodes_outer+i+1;
end
for i = 1:3
    Links(num_links_outer+3+i,1) = num_nodes_outer+4+i;
    Links(num_links_outer+3+i,2) = num_nodes_outer+4+i+1;
end
for i = 1:N2
    Links(num_links_outer+6+i,1) = num_nodes_outer+1;
    Links(num_links_outer+6+i,2) = num_nodes_outer+5;
end

%% Anchor
Links(num_links_outer+num_links_inner+1,1) = 4;
Links(num_links_outer+num_links_inner+1,2) = num_nodes_outer+num_nodes_inner+1;
for i = 1:4
    Links(num_links_outer+num_links_inner+1+i,1) = num_nodes_outer+num_nodes_inner+1+(i-1);
    Links(num_links_outer+num_links_inner+1+i,2) = num_nodes_outer+num_nodes_inner+1+i;    
end
Links(num_links_outer+num_links_inner+6,1) = num_nodes_outer+num_nodes_inner+5;
Links(num_links_outer+num_links_inner+6,2) = num_nodes_outer+1;
%--------------------------------------------------------------------------
Links(num_links_outer+num_links_inner+7,1) = 11;
Links(num_links_outer+num_links_inner+7,2) = num_nodes_outer+num_nodes_inner+6;
for i = 1:4
    Links(num_links_outer+num_links_inner+7+i,1) = num_nodes_outer+num_nodes_inner+6+(i-1);
    Links(num_links_outer+num_links_inner+7+i,2) = num_nodes_outer+num_nodes_inner+6+i;    
end
Links(num_links_outer+num_links_inner+12,1) = num_nodes_outer+num_nodes_inner+10;
Links(num_links_outer+num_links_inner+12,2) = num_nodes_outer+5;

%% Geometry %%
Geometry = zeros(num_links,3);
Geometry(:,3) = TD;

%% Outer
Geometry(1,1) = abs(OOB(2,1));
Geometry(1,2) = OOB(3,2)-OOB(17,2);

for i = 1:N1
    Geometry(1+i,1) = GeoVariables(10);
    Geometry(1+i,2) = GeoVariables(11);
end

Geometry(1+N1+1,1) = abs(OOB(16,2));
Geometry(1+N1+1,2) = abs(OOB(4,1)-OOB(16,1));

Geometry(1+N1+2,1) = abs(OOB(15,2));
Geometry(1+N1+2,2) = abs(OOB(4,1)-OOB(16,1));

for i = 1:N3
    Geometry(1+N1+2+i,1) = GeoVariables(6);
    Geometry(1+N1+2+i,2) = GeoVariables(7);
end

Geometry(1+N1+2+N3+1,1) = abs(OOB(12,1)-OOB(11,1));
Geometry(1+N1+2+N3+1,2) = abs(OOB(8,2)-OOB(11,2));

Geometry(1+N1+2+N3+2,1) = (OOB(8,2)-OOB(9,2));
Geometry(1+N1+2+N3+2,2) = abs(OOB(8,1)-OOB(11,1));

%--------------------------------------------------------------------------
Geometry(1+N1+2+N3+3:num_links_outer-1,1:3) = Geometry(1:1+N1+2+N3+2,1:3);
%--------------------------------------------------------------------------
Geometry(num_links_outer,1) = OOB(2,2) - OOB(3,2);
Geometry(num_links_outer,2) = 2*abs(OOB(2,1));

%% Inner
Geometry(num_links_outer+1,1) = abs(OIB(5,1)-OIB(6,1));
Geometry(num_links_outer+1,2) = OIB(12,2)-OIB(5,2);
Geometry(num_links_outer+2,1) = OIB(6,2)-OIB(7,2);
Geometry(num_links_outer+2,2) = OIB(12,2)-OIB(5,2);
Geometry(num_links_outer+3,1) = OIB(8,2)-OIB(9,2);
Geometry(num_links_outer+3,2) = abs(OIB(8,1)-OIB(10,1));
%--------------------------------------------------------------------------
Geometry(num_links_outer+4:num_links_outer+6,1:3) = Geometry(num_links_outer+1:num_links_outer+3,1:3);
%--------------------------------------------------------------------------
for i = 1:N2
    Geometry(num_links_outer+6+i,1) = 2*GeoVariables(2);
    Geometry(num_links_outer+6+i,2) = GeoVariables(3);
end

%% Anchor
Geometry(num_links_outer+num_links_inner+1,1) = TO;
Geometry(num_links_outer+num_links_inner+1,2) = abs(OHLl(1,1)-OHLl(9,1));
Geometry(num_links_outer+num_links_inner+1,3) = abs(OHLl(1,2)-OHLl(2,2));

Geometry(num_links_outer+num_links_inner+2,1) = TH;
Geometry(num_links_outer+num_links_inner+2,2) = abs(OHLl(1,1)-OHLl(9,1));
Geometry(num_links_outer+num_links_inner+2,3) = abs(OHLl(1,2)-OHLl(2,2));

Geometry(num_links_outer+num_links_inner+3,1) = abs(OHLl(8,1)-OHLl(9,1));
Geometry(num_links_outer+num_links_inner+3,2) = abs(OHLl(8,2)-OHLl(5,2));
Geometry(num_links_outer+num_links_inner+3,3) = TH;

Geometry(num_links_outer+num_links_inner+4,1) = TH;
Geometry(num_links_outer+num_links_inner+4,2) = abs(OHLl(5,1)-OHLl(6,1));
Geometry(num_links_outer+num_links_inner+4,3) = 1/2*(abs(OHLl(8,2)-OOB(5,2))+abs(OOB(7,2)-OOB(6,2)));

Geometry(num_links_outer+num_links_inner+5,1) = TO;
Geometry(num_links_outer+num_links_inner+5,2) = abs(OHLl(5,1)-OHLl(6,1));
Geometry(num_links_outer+num_links_inner+5,3) = 1/2*abs(OHLl(8,2)-OOB(5,2))+abs(OOB(7,2)-OOB(6,2));

Geometry(num_links_outer+num_links_inner+6,1) = 1/2*TD;
Geometry(num_links_outer+num_links_inner+6,2) = abs(OHLl(5,1)-OHLl(6,1));
Geometry(num_links_outer+num_links_inner+6,3) = 1/2*abs(OHLl(8,2)-OOB(5,2))+abs(OOB(7,2)-OOB(6,2));
%--------------------------------------------------------------------------
Geometry(num_links_outer+num_links_inner+7:num_links_outer+num_links_inner+12,1:3) = Geometry(num_links_outer+num_links_inner+1:num_links_outer+num_links_inner+6,1:3);

%% Convert
% Units mm to m
Geometry(:,1:3) = Geometry(:,1:3).*(10^(-3));

% Lengths
L = Geometry(:,1);

% Cross-Section Areas
A = zeros(num_links,1);
for i = 1:num_links
    A(i,1) = Geometry(i,2)*Geometry(i,3);
end

%% Conductivity %%
Cond(:,1) = k;

%% Anchor
Cond(num_links_outer+num_links_inner+1,1) = k_O;
Cond(num_links_outer+num_links_inner+5,1) = k_O;
Cond(num_links_outer+num_links_inner+7,1) = k_O;
Cond(num_links_outer+num_links_inner+num_links_anchor-1,1) = k_O;

%% Solve %%
run('Solve_E.m');
run('Solve_C.m');