%% Set Parameters
TD = 0.4;
TO = 0.004;
TH = 0.05;
Gtb = 0.01;
Gcb = 0.01;
Gbb = 0.01;
Gmc = 0.05;

E = 1.60e11;
v = 0.22;
G = E/(2*(1+v));

%% Parameters to Input
% Number of beams
N1 = GeoVariables(9); N2 = GeoVariables(1); N3 = GeoVariables(5);

% Angles
theta1 = deg2rad(GeoVariables(8)); theta2 = -deg2rad(GeoVariables(4));

% Beam Widths
h(1) = GeoVariables(11); h(2) = GeoVariables(3); h(3) = GeoVariables(7);

% Nodes
num_nodes_outer = N1+1+N3+N1+N3+2;

if mod(N2,2) == 0
    num_nodes_inner = 1+N2+N2+3;
else
    num_nodes_inner = N2+N2+3;
end

num_nodes = num_nodes_outer+num_nodes_inner + 5;

% Links
num_links_outer = (N1-1)+2+(N3-1)+N1+N3+(N1-1)+2;

if mod(N2,2) == 0
    num_links_inner = (N2-1)+1+N2+3;
else
    num_links_inner = (N2-1)+N2+3;
end

num_links = num_links_outer+num_links_inner + 6;

%% Coordinates
Coordinates = zeros(num_nodes,3);

% Outer
for i = 1:N1
    Coordinates(i,1) = OOB(16,1) - 1/2*h(1)*sin(theta1) - (N1-i)*(h(1)+Gtb)*sin(theta1);
    Coordinates(i,2) = OOB(16,2) + 1/2*h(1)*cos(theta1) + (N1-i)*(h(1)+Gtb)*cos(theta1);
end
for i = 1:1
    Coordinates(N1+i,1) = Coordinates(N1,1);
    Coordinates(N1+i,2) = 0;
end
for i = 1:N3
    Coordinates(N1+1+i,1) = Coordinates(N1,1) + (i-1)*(h(3)+Gbb)*sin(theta2);
    Coordinates(N1+1+i,2) = OOB(15,2) - 1/2*h(3)*cos(theta2) - (i-1)*(h(3)+Gbb)*cos(theta2);
end
for i = 1:N1
    Coordinates(N1+1+N3+i,1) = OOB(17,1) - 1/2*h(1)*sin(theta1) - (N1-i)*(h(1)+Gtb)*sin(theta1);
    Coordinates(N1+1+N3+i,2) = OOB(17,2) + 1/2*h(1)*cos(theta1) + (N1-i)*(h(1)+Gtb)*cos(theta1);
end
for i = 1:N3
    Coordinates(N1+1+N3+N1+i,1) = OOB(14,1) + 1/2*h(3)*sin(theta2) + (i-1)*(h(3)+Gbb)*sin(theta2);
    Coordinates(N1+1+N3+N1+i,2) = OOB(14,2) - 1/2*h(3)*cos(theta2) - (i-1)*(h(3)+Gbb)*cos(theta2);
end
for i = 1:1
    Coordinates(N1+1+N3+N1+N3+i,1) = 0;
    Coordinates(N1+1+N3+N1+N3+i,2) = Coordinates(N1+1+N3+N1,2);
end
for i = 1:1
    Coordinates(N1+1+N3+N1+N3+1+i,1) = Coordinates(N1+1+N3+1,1);
    Coordinates(N1+1+N3+N1+N3+1+i,2) = OOB(2,2);
end

% Inner
if mod(N2,2) == 0
    for i = 1:1
        Coordinates(num_nodes_outer+i,1) = OIB(2,1);
        Coordinates(num_nodes_outer+i,2) = 0;
    end
    for i = 1:N2/2
        Coordinates(num_nodes_outer+1+i,1) = OIB(2,1);
        Coordinates(num_nodes_outer+1+i,2) = 1/2*(Gcb+h(2))+(N2/2-i)*(Gcb+h(2));
    end
    for i = 1:N2/2
        Coordinates(num_nodes_outer+1+N2/2+i,1) = OIB(2,1);
        Coordinates(num_nodes_outer+1+N2/2+i,2) = -1/2*(Gcb+h(2))-(i-1)*(Gcb+h(2));
    end
    for i = 1:N2
        Coordinates(num_nodes_outer+1+N2+i,1) = 0;
        Coordinates(num_nodes_outer+1+N2+i,2) = Coordinates(num_nodes_outer+1+i,2);
    end    
    for i = 1:1
        Coordinates(num_nodes_outer+1+N2+N2+i,1) = OIB(2,1);
        Coordinates(num_nodes_outer+1+N2+N2+i,2) = OIB(12,2)-1/2*(OIB(12,2)-OIB(5,2));
    end
    for i = 1:1
        Coordinates(num_nodes_outer+1+N2+N2+1+i,1) = OIB(6,1)+1/2*(OIB(12,2)-OIB(5,2));
        Coordinates(num_nodes_outer+1+N2+N2+1+i,2) = Coordinates(num_nodes_outer+1+N2+N2+1,2);
    end
    for i = 1:1
        Coordinates(num_nodes_outer+1+N2+N2+2+i,1) = OIB(6,1)+1/2*(OIB(12,2)-OIB(5,2));
        Coordinates(num_nodes_outer+1+N2+N2+2+i,2) = OHLb(2,2);
    end
else
    for i = 1:N2
        Coordinates(num_nodes_outer+i,1) = OIB(2,1);
        Coordinates(num_nodes_outer+i,2) = (((N2+1)/2)-i)*(Gcb+h(2));
    end
    for i = 1:N2
        Coordinates(num_nodes_outer+N2+i,1) = 0;
        Coordinates(num_nodes_outer+N2+i,2) = Coordinates(num_nodes_outer+i,2);
    end
    for i = 1:1
        Coordinates(num_nodes_outer+N2+N2+i,1) = OIB(2,1);
        Coordinates(num_nodes_outer+N2+N2+i,2) = OIB(12,2)-1/2*(OIB(12,2)-OIB(5,2));
    end
    for i = 1:1
        Coordinates(num_nodes_outer+N2+N2+1+i,1) = OIB(6,1)+1/2*(OIB(12,2)-OIB(5,2));
        Coordinates(num_nodes_outer+N2+N2+1+i,2) = Coordinates(num_nodes_outer+N2+N2+1,2);
    end
    for i = 1:1
        Coordinates(num_nodes_outer+N2+N2+2+i,1) = OIB(6,1)+1/2*(OIB(12,2)-OIB(5,2));
        Coordinates(num_nodes_outer+N2+N2+2+i,2) = OHLb(2,2);
    end
end

% Anchor
for i = 1:1
    Coordinates(num_nodes_outer+num_nodes_inner+i,1) = Coordinates(N1+1,1);
    Coordinates(num_nodes_outer+num_nodes_inner+i,2) = Coordinates(N1+1,2);
    Coordinates(num_nodes_outer+num_nodes_inner+i,3) = -1/2*TH;
end
for i = 1:1
    Coordinates(num_nodes_outer+num_nodes_inner+1+i,1) = Coordinates(N1+1,1)+Gmc;
    Coordinates(num_nodes_outer+num_nodes_inner+1+i,2) = Coordinates(N1+1,2);
    Coordinates(num_nodes_outer+num_nodes_inner+1+i,3) = -1/2*TH;
end
for i = 1:1
    Coordinates(num_nodes_outer+num_nodes_inner+2+i,1) = Coordinates(N1+1,1)+Gmc;
    Coordinates(num_nodes_outer+num_nodes_inner+2+i,2) = Coordinates(N1+1,2);
end
for i = 1:2
    Coordinates(num_nodes_outer+num_nodes_inner+3+i,1) = Coordinates(num_nodes_outer+num_nodes_inner+3,1) + i*abs((OIB(3,1)-OIB(2,1)))/3;
    Coordinates(num_nodes_outer+num_nodes_inner+3+i,2) = Coordinates(N1+1,2);
end

%% DOFs
DOFs = zeros(num_nodes,6);

% Outer
for i = 1:N1+1+N3
    DOFs(i,1:2) = 1;
    DOFs(i,6) = 1;
end
for i = 1:N1
    DOFs(N1+1+N3+i,1:2) = 1;
end
for i = 1:N3
    DOFs(N1+1+N3+N1+i,1:6) = 0;
end
for i = 1:1
    DOFs(N1+1+N3+N1+N3+i,2) = 1;
end
for i = 1:1
    DOFs(N1+1+N3+N1+N3+1+i,1:2) = 1;
end

% Inner
if mod(N2,2) == 0
    for i = 1:1+N2
        DOFs(num_nodes_outer+i,1:2) = 1;
        DOFs(num_nodes_outer+i,6) = 1;
    end
    for i = 1:N2
        DOFs(num_nodes_outer+1+N2+i,2) = 1;
    end
    for i = 1:2
        DOFs(num_nodes_outer+1+N2+N2+i,1:2) = 1;
        DOFs(num_nodes_outer+1+N2+N2+i,6) = 1;
    end
else
    for i = 1:N2
        DOFs(num_nodes_outer+i,1:2) = 1;
        DOFs(num_nodes_outer+i,6) = 1;
    end
    for i = 1:N2
        DOFs(num_nodes_outer+N2+i,2) = 1;
    end
    for i = 1:2
        DOFs(num_nodes_outer+N2+N2+i,1:2) = 1;
        DOFs(num_nodes_outer+N2+N2+i,6) = 1;
    end
end

% Anchor
for i = 1:5
    DOFs(num_nodes_outer+num_nodes_inner+i,1:2) = 1;
    DOFs(num_nodes_outer+num_nodes_inner+i,6) = 1;
end

%% Links
Links = zeros(num_links,2);

% Outer
for i = 1:(N1-1)+2+(N3-1)
    Links(i,1) = i;
    Links(i,2) = i+1;
end
for i = 1:N1
    Links((N1-1)+2+(N3-1)+i,1) = i;
    Links((N1-1)+2+(N3-1)+i,2) = N1+1+N3+i;
end
for i = 1:N3
    Links((N1-1)+2+(N3-1)+N1+i,1) = N1+1+i;
    Links((N1-1)+2+(N3-1)+N1+i,2) = N1+1+N3+N1+i;
end
for i = 1:N1-1
    Links((N1-1)+2+(N3-1)+N1+N3+i,1) = N1+1+N3+i;
    Links((N1-1)+2+(N3-1)+N1+N3+i,2) = N1+1+N3+i+1;
end
for i = 1:1
    Links((N1-1)+2+(N3-1)+N1+N3+(N1-1)+i,1) = N1+1+N3+N1;
    Links((N1-1)+2+(N3-1)+N1+N3+(N1-1)+i,2) = N1+1+N3+N1+N3+i;
end
for i = 1:1
    Links((N1-1)+2+(N3-1)+N1+N3+(N1-1)+1+i,1) = N1+1+N3+i;
    Links((N1-1)+2+(N3-1)+N1+N3+(N1-1)+1+i,2) = N1+1+N3+N1+N3+1+i;
end

% Inner
if mod(N2,2) == 0
    for i = 1:N2/2-1
        Links(num_links_outer+i,1) = num_nodes_outer+1+i;
        Links(num_links_outer+i,2) = num_nodes_outer+1+i+1;
    end
    for i = 1:1
        Links(num_links_outer+(N2/2-1)+i,1) = num_nodes_outer+1+N2/2;
        Links(num_links_outer+(N2/2-1)+i,2) = num_nodes_outer+1;
    end
    for i = 1:1
        Links(num_links_outer+(N2/2-1)+1+i,1) = num_nodes_outer+1;
        Links(num_links_outer+(N2/2-1)+1+i,2) = num_nodes_outer+1+N2/2+1;
    end
    for i = 1:N2/2-1
        Links(num_links_outer+(N2/2-1)+2+i,1) = num_nodes_outer+1+N2/2+i;
        Links(num_links_outer+(N2/2-1)+2+i,2) = num_nodes_outer+1+N2/2+i+1;
    end
    for i = 1:N2
        Links(num_links_outer+N2+i,1) = num_nodes_outer+1+i;
        Links(num_links_outer+N2+i,2) = num_nodes_outer+1+N2+i;
    end
    for i = 1:1
        Links(num_links_outer+N2+N2+i,1) = num_nodes_outer+1+N2;
        Links(num_links_outer+N2+N2+i,2) = num_nodes_outer+1+N2+N2+i;
    end
    for i = 1:2
        Links(num_links_outer+N2+N2+1+i,1) = num_nodes_outer+1+N2+N2+i;
        Links(num_links_outer+N2+N2+1+i,2) = num_nodes_outer+1+N2+N2+i+1;
    end
else
    for i = 1:N2-1
        Links(num_links_outer+i,1) = num_nodes_outer+i;
        Links(num_links_outer+i,2) = num_nodes_outer+i+1;
    end
    for i = 1:N2
        Links(num_links_outer+(N2-1)+i,1) = num_nodes_outer+i;
        Links(num_links_outer+(N2-1)+i,2) = num_nodes_outer+N2+i;
    end
    for i = 1:1
        Links(num_links_outer+(N2-1)+N2+i,1) = num_nodes_outer+N2;
        Links(num_links_outer+(N2-1)+N2+i,2) = num_nodes_outer+N2+N2+i;
    end
    for i = 1:2
        Links(num_links_outer+(N2-1)+N2+1+i,1) = num_nodes_outer+N2+N2+i;
        Links(num_links_outer+(N2-1)+N2+1+i,2) = num_nodes_outer+N2+N2+i+1;
    end
end

% Anchor
for i = 1:1
    Links(num_links_outer+num_links_inner+i,1) = N1+1;
    Links(num_links_outer+num_links_inner+i,2) = num_nodes_outer+num_nodes_inner+1;
end

for i = 1:4
    Links(num_links_outer+num_links_inner+1+i,1) = num_nodes_outer+num_nodes_inner+i;
    Links(num_links_outer+num_links_inner+1+i,2) = num_nodes_outer+num_nodes_inner+i+1;
end    

if mod(N2,2) == 0
    for i = 1:1
        Links(num_links_outer+num_links_inner+5+i,1) = num_nodes_outer+num_nodes_inner+5;
        Links(num_links_outer+num_links_inner+5+i,2) = num_nodes_outer+1;
    end
else
    for i = 1:1
        Links(num_links_outer+num_links_inner+5+i,1) = num_nodes_outer+num_nodes_inner+5;
        Links(num_links_outer+num_links_inner+5+i,2) = num_nodes_outer+(N2+1)/2;
    end
end

%% Geometry
Geometry = zeros(num_links,3); % L B H

for i = 1:num_links
    Geometry(i,1) = sqrt((Coordinates(Links(i,2),1)-Coordinates(Links(i,1),1))^2 + (Coordinates(Links(i,2),2)-Coordinates(Links(i,1),2))^2 + (Coordinates(Links(i,2),3)-Coordinates(Links(i,1),3))^2);
end
for i = 1:num_links
    Geometry(i,2) = TD;
end

% Outer
for i = 1:(N1-1)+2+(N3-1)
    Geometry(i,3) = abs(OOB(5,1)-OOB(15,1));
end
for i = 1:N1
    Geometry((N1-1)+2+(N3-1)+i,3) = h(1);
end
for i = 1:N3
    Geometry((N1-1)+2+(N3-1)+N1+i,3) = h(3);
end
for i = 1:N1-1
    Geometry((N1-1)+2+(N3-1)+N1+N3+i,3) = abs(OOB(17,1));
end
for i = 1:1
    Geometry((N1-1)+2+(N3-1)+N1+N3+(N1-1)+i,3) = abs(OOB(2,2)-OOB(17,2));
end
for i = 1:1
    Geometry((N1-1)+2+(N3-1)+N1+N3+(N1-1)+1+i,3) = abs(OOB(2,1));
end

% Inner
if mod(N2,2) == 0
    for i = 1:N2
    Geometry(num_links_outer+i,3) = (abs(OIB(2,1)-OIB(3,1)));
    end
    for i = 1:N2
    Geometry(num_links_outer+N2+i,3) = h(2);
    end    
    for i = 1:1
    Geometry(num_links_outer+N2+N2+i,3) = (abs(OIB(2,1)-OIB(3,1)));
    end      
    for i = 1:2
    Geometry(num_links_outer+N2+N2+1+i,3) = abs(OIB(12,2)-OIB(5,2));
    end
else
    for i = 1:N2-1
    Geometry(num_links_outer+i,3) = (abs(OIB(2,1)-OIB(3,1)));
    end
    for i = 1:N2
    Geometry(num_links_outer+(N2-1)+i,3) = h(2);
    end    
    for i = 1:1
    Geometry(num_links_outer+(N2-1)+N2+i,3) = (abs(OIB(2,1)-OIB(3,1)));
    end      
    for i = 1:2
    Geometry(num_links_outer+(N2-1)+N2+1+i,3) = abs(OIB(12,2)-OIB(5,2));
    end    
end

% Anchor
for i = 1:1
    Geometry(num_links_outer+num_links_inner+i,2) = abs(OHLl(1,1) - OHLl(10,1));
    Geometry(num_links_outer+num_links_inner+i,3) = abs(OHLl(1,2) - OHLl(2,2));
end
for i = 1:1
    Geometry(num_links_outer+num_links_inner+1+i,2) = TH;
    Geometry(num_links_outer+num_links_inner+1+i,3) = abs(OOB(3,2)-OOB(4,2));
end
for i = 1:1
    Geometry(num_links_outer+num_links_inner+2+i,2) = abs(OHLl(8,1) - OHLl(7,1));
    Geometry(num_links_outer+num_links_inner+2+i,3) = 1/2*(abs(OHLl(8,2) - OHLl(5,2))+abs(OHLl(7,2) - OHLl(6,2)));
end
for i = 1:3
    LL = OIB(3,2)-OIB(4,2)-(i-1)*(OIB(3,2)-OIB(4,2)-(OIB(2,2)-OIB(5,2)))/3;
    RR = OIB(3,2)-OIB(4,2)-(i)*(OIB(3,2)-OIB(4,2)-(OIB(2,2)-OIB(5,2)))/3;
    Geometry(num_links_outer+num_links_inner+3+i,3) = (LL+RR)/2;
end

% Units
Geometry = Geometry*10^-3;

%% Stiffness Matrix Parameters
% Effective Shear
ky = 1.2;

% Length
L = Geometry(:,1);

% Cross Sectional Area
A = Geometry(:,2).*Geometry(:,3);

% Moment of Inertia
I_z = 1/12*Geometry(:,2).*Geometry(:,3).^3;
I_y = 1/12*Geometry(:,3).*Geometry(:,2).^3;

% Torsion Constant
J = zeros(num_links,1);
Long = zeros(num_links,1); Short = zeros(num_links,1); 
for i = 1:num_links
    Long(i,1) = max(Geometry(i,2),Geometry(i,3)); Short(i,1) = min(Geometry(i,2),Geometry(i,3)); 
    J(i,1) = Long(i,1)*Short(i,1)^3*(1/3-0.21*Short(i,1)/Long(i,1)*(1-Short(i,1)^4/(12*Long(i,1)^4)));
end

% Direction Cosines
Coordinates = Coordinates*10^-3; % Change Coordinates Units

LXZ = zeros(num_links,1);
for i = 1:num_links
    LXZ(i,1) = sqrt((Coordinates(Links(i,2),1)-Coordinates(Links(i,1),1))^2 + (Coordinates(Links(i,2),3)-Coordinates(Links(i,1),3))^2);
end

CxX = zeros(num_links,1); CxY = zeros(num_links,1); CxZ = zeros(num_links,1); CzX = zeros(num_links,1); CzZ = zeros(num_links,1); CyX = zeros(num_links,1); CyZ = zeros(num_links,1); CyY = zeros(num_links,1); 

for i = 1:num_links
    CxX(i) = (Coordinates(Links(i,2),1) - Coordinates(Links(i,1),1))/L(i);
    CxY(i) = (Coordinates(Links(i,2),2) - Coordinates(Links(i,1),2))/L(i);
    CxZ(i) = (Coordinates(Links(i,2),3) - Coordinates(Links(i,1),3))/L(i);    
    CzX(i) = -(Coordinates(Links(i,2),3) - Coordinates(Links(i,1),3))/LXZ(i);
    CzZ(i) = (Coordinates(Links(i,2),1) - Coordinates(Links(i,1),1))/LXZ(i);
    CyX(i) = -CzZ(i)*CxY(i);
    CyZ(i) = CzX(i)*CxZ(i);
    CyY(i) = sqrt((Coordinates(Links(i,2),1) - Coordinates(Links(i,1),1))^2 + (Coordinates(Links(i,2),3) - Coordinates(Links(i,1),3))^2)/L(i);
end
    
Dir_Cos = zeros(3,3,num_links);

for i = 1:num_links
    if Coordinates(Links(i,2),1)-Coordinates(Links(i,1),1) == 0 && Coordinates(Links(i,2),3)-Coordinates(Links(i,1),3) == 0
        if Coordinates(Links(i,2),2) > Coordinates(Links(i,1),2)
            Dir_Cos(1,1,i) = 0; Dir_Cos(1,2,i) = 1; Dir_Cos(1,3,i) = 0;
            Dir_Cos(2,1,i) = -1; Dir_Cos(2,2,i) = 0; Dir_Cos(2,3,i) = 0;
            Dir_Cos(3,1,i) = 0; Dir_Cos(3,2,i) = 0; Dir_Cos(3,3,i) = 1;
        else
            Dir_Cos(1,1,i) = 0; Dir_Cos(1,2,i) = -1; Dir_Cos(1,3,i) = 0;
            Dir_Cos(2,1,i) = 1; Dir_Cos(2,2,i) = 0; Dir_Cos(2,3,i) = 0;
            Dir_Cos(3,1,i) = 0; Dir_Cos(3,2,i) = 0; Dir_Cos(3,3,i) = 1;
        end
    else
    Dir_Cos(1,1,i) = CxX(i); Dir_Cos(1,2,i) = CxY(i); Dir_Cos(1,3,i) = CxZ(i);
    Dir_Cos(2,1,i) = CyX(i); Dir_Cos(2,2,i) = CyY(i); Dir_Cos(2,3,i) = CyZ(i);
    Dir_Cos(3,1,i) = CzX(i); Dir_Cos(3,2,i) = 0; Dir_Cos(3,3,i) = CzZ(i);
    end
end

% Transformation Matrices
Transform_M = zeros(12,12,num_links);

for i = 1:num_links
    Transform_M(1:3,1:3,i) = Dir_Cos(:,:,i);
    Transform_M(4:6,4:6,i) = Dir_Cos(:,:,i);
    Transform_M(7:9,7:9,i) = Dir_Cos(:,:,i);
    Transform_M(10:12,10:12,i) = Dir_Cos(:,:,i);
end

% phi
phi_y = zeros(num_links,1);
phi_z = zeros(num_links,1);

for i = 1:num_links
    phi_y(i,1) = 12*E*I_z(i,1)*ky/(A(i)*G*(L(i))^2);
    phi_z(i,1) = 12*E*I_y(i,1)*ky/(A(i)*G*(L(i))^2);
end

%% Stiffness Matrices
% Timoshenko
X = A*E./L;

Y1 = 12*E*I_z./((1+phi_y).*L.^3);
Y2 = 6*E*I_z./((1+phi_y).*L.^2);
Y3 = (4+phi_y)*E.*I_z./((1+phi_y).*L);
Y4 = (2-phi_y)*E.*I_z./((1+phi_y).*L);

Z1 = 12*E*I_y./((1+phi_z).*L.^3);
Z2 = 6*E*I_y./((1+phi_z).*L.^2);
Z3 = (4+phi_z)*E.*I_y./((1+phi_z).*L);
Z4 = (2-phi_z)*E.*I_y./((1+phi_z).*L);

S = G*J./L;

% Element
K = zeros(12,12,num_links);
K_Global = zeros(12,12,num_links);

for i = 1:num_links
    K(1,:,i) = [X(i),0,0,0,0,0,-X(i),0,0,0,0,0];
    K(2,:,i) = [0,Y1(i),0,0,0,Y2(i),0,-Y1(i),0,0,0,Y2(i)];
    K(3,:,i) = [0,0,Z1(i),0,-Z2(i),0,0,0,-Z1(i),0,-Z2(i),0];
    K(4,:,i) = [0,0,0,S(i),0,0,0,0,0,-S(i),0,0];
    K(5,:,i) = [0,0,0,0,Z3(i),0,0,0,Z2(i),0,Z4(i),0];
    K(6,:,i) = [0,0,0,0,0,Y3(i),0,-Y2(i),0,0,0,Y4(i)];
    K(7,:,i) = [0,0,0,0,0,0,X(i),0,0,0,0,0];
    K(8,:,i) = [0,0,0,0,0,0,0,Y1(i),0,0,0,-Y2(i)];
    K(9,:,i) = [0,0,0,0,0,0,0,0,Z1(i),0,Z2(i),0];
    K(10,:,i) = [0,0,0,0,0,0,0,0,0,S(i),0,0];
    K(11,:,i) = [0,0,0,0,0,0,0,0,0,0,Z3(i),0];
    K(12,:,i) = [0,0,0,0,0,0,0,0,0,0,0,Y3(i)];
    
    K(:,:,i) = (K(:,:,i))'+(K(:,:,i));
    n = 1;
    for ii = 1:12
        K(n,n,i) = K(n,n,i)/2;
        n = n+1;
    end
 
    K_Global(:,:,i) = transpose(Transform_M(:,:,i))*K(:,:,i)*Transform_M(:,:,i);
end

% Global
KK = zeros(num_nodes*6,num_nodes*6,num_links);
K_ff = zeros(num_nodes*6,num_nodes*6);

for i = 1:num_links
    SN = Links(i,1); EN = Links(i,2);
    KK( ((SN-1)*6+1):SN*6 , ((SN-1)*6+1):SN*6 , i ) = K_Global(1:6,1:6,i);
    KK( ((SN-1)*6+1):SN*6 , ((EN-1)*6+1):EN*6 , i ) = K_Global(1:6,7:12,i);
    KK( ((EN-1)*6+1):EN*6 , ((SN-1)*6+1):SN*6 , i ) = K_Global(7:12,1:6,i);
    KK( ((EN-1)*6+1):EN*6 , ((EN-1)*6+1):EN*6 , i ) = K_Global(7:12,7:12,i);
    K_ff = K_ff + KK(:,:,i);
end    
    
% Rearrange
N = 1;
for i = 1:num_nodes
    for ii = 1:6
        if DOFs(i,ii) == 0
           K_ff(N,:) = [];
           K_ff(:,N) = [];
        else
           N = N+1;
        end
    end
end

%% Solve
run('Solve_Disp.m');
run('Solve_Stiff.m');

