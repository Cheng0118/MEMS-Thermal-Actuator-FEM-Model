function [Inequality, Equality] = GeoConstraint(GeoVariables)
% In this function we specify the geometric constraints and check their values. All Geometric dimensions has the unit of mm.
% Input "GeoVariables" is a 1x11 vector:
%   Ncb=GeoVariables(1);                    % Number of center beams
%   Lcb=GeoVariables(2);                    % Length of center beams
%   Wcb=GeoVariables(3);                    % Width of center beams
% 
%   Abb=deg2rad(GeoVariables(4));           % Angle of bottom beams
%   Nbb=GeoVariables(5);                    % Number of bottom beams
%   Lbb=GeoVariables(6);                    % Length of bottom beams
%   Wbb=GeoVariables(7);                    % Width of bottom beams
% 
%   Atb=deg2rad(GeoVariables(8));           % Angle of top beams
%   Ntb=GeoVariables(9);                    % Number of top beams
%   Ltb=GeoVariables(10);                   % Length of top beams
%   Wtb=GeoVariables(11);                   % Width of top beams

% Input Example:
% GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02] % Case 1
% GeoVariables=[2,0.7,0.04,     10,3,0.4,0.03,     20,3,0.9,0.02] % Case 2

% Output
% "Inequality" 
%       is a 5x1 vector corresponding to 5 constraints. It requires that all 5 values should be smaller than 0. Otherwise the input
%       "GeoVariables" cannot produce a a CAD design. See Page 11 of ducument.
% "Equality"
%       This is always an empty vector.




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Please DO NOT edit this part
GeoVariables=GeoVariables(:);
Resolution= [1,0.01,0.01,               1,1,0.01,0.01,       1,1,0.01,0.01];	% Feasible resolution
GeoVariables=Round2Res(GeoVariables,Resolution);
Equality = [];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% In this part you specify your variables and specify you inequality constraints.
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;     	% Anchor width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;            % Thickness of the device
Wtop=0.05;          % Width of the top plate
LMmin=0.625;       	% The minimum LM width
x=1;y=2;            % For a clearer notation we define x and y.

Ncb=GeoVariables(1);                    % Number of center beams
Lcb=GeoVariables(2);                    % Length of center beams
Wcb=GeoVariables(3);                    % Width of center beams

Abb=deg2rad(GeoVariables(4));           % Angle of bottom beams
Nbb=GeoVariables(5);                    % Number of bottom beams
Lbb=GeoVariables(6);                    % Length of bottom beams
Wbb=GeoVariables(7);                    % Width of bottom beams

Atb=deg2rad(GeoVariables(8));           % Angle of top beams
Ntb=GeoVariables(9);                    % Number of top beams
Ltb=GeoVariables(10);                   % Length of top beams
Wtb=GeoVariables(11);                   % Width of top beams

% Calcaulte the total width of three sets of beams
WB=Nbb*Wbb+(Nbb-1)*Gbb;
WC=Ncb*Wcb+(Ncb-1)*Gcb;
WT=Ntb*Wtb+(Ntb-1)*Gtb;

% Calculate the corner point of Inner beam
InnerBeamO1 =       [Lcb,           WC/2];
InnerBeamO22=       [Lcb,           -WC/2-Gmg-Wmc];        

% Calculate some basic geometric dimensions
Anchor=max([Anchormin,WT*sin(Atb)+Anchormin*cos(Atb),WB*sin(Abb)+Anchormin*cos(Abb)]);
CW=max([Ltb*cos(Atb)+Anchormin, Lbb*cos(Abb)+Wsc+Gmc,InnerBeamO1(x)+Anchormin+Gmc]);% Guarantee there is enough space for: 1. Top anchor width, 2. bottom Anchor gap 3. Center beams and anchor gap

Ld=CW+Anchor+Gmg+Wsb;
LM=Ld-Td;
WA=HH-Wtop-2*Gmg-Wmc-Ltb*sin(Atb)-WT*cos(Atb)-Ld-Td*tan(Abb)-Wsc*cos(Abb)+(Wsb-Wsc*sin(Abb))*tan(Abb)+Gmg*tan(Abb)-Gmg+(Anchor-WB*sin(Abb))*tan(Abb)-WB*cos(Abb);
% WA is calcualted by substracting all other vertical components from the total device height HH.
Ymax=WA/2+Wtop+2*Gmg+Wmc+Ltb*sin(Atb)+WT*cos(Atb);
Ymin=Ymax-HH;

OuterBeamO28x=CW-Lbb*cos(Abb)+WB*sin(Abb);
OuterBeamO27x=LM-Wsc-Gmc;

if OuterBeamO28x>=OuterBeamO27x
    Case=2;
else
    Case=1;
end


OuterBeamO19=      [CW,                  WA/2];
OuterBeamO20=      [CW,                 -WA/2];
OuterBeamO21=      [CW-Lbb*cos(Abb),    -WA/2-Lbb*sin(Abb)];
OuterBeamO31=      [CW+Anchor,          WA/2+Ltb*sin(Atb)+WT*cos(Atb)-(Anchor-WT*sin(Atb)+Ltb*cos(Atb))*tan(Atb)];
OuterBeamO30=      [CW+Anchor,          -WA/2-Lbb*sin(Abb)-WB*cos(Abb)+(Anchor-WB*sin(Abb)+Lbb*cos(Abb))*tan(Abb)];

InnerBeamO19x=(LM-Wsc-Gmc-Gmc/2)/2;
if Case==1
    Constraint6=-Anchormin;
else
    Constraint6=Anchormin-(OuterBeamO21(x)-InnerBeamO19x-Gmc)/cos(Abb);%<=0
end

Inequality =        [2*Gbb-(OuterBeamO19(y)+(CW-InnerBeamO1(x))*tan(Atb)-InnerBeamO1(y));%<=0                           %Constraint (1)                   
                    2*Gbb-(InnerBeamO22(y)-(OuterBeamO20(y)-(CW-InnerBeamO22(x))*tan(Abb)));%<=0                        %Constraint (2)
                    Anchormin-WA;%<=0                                                                                 	%Constraint (3)
                    LMmin-LM;%<=0                                                                                       %Constraint (4)
                    OuterBeamO30(y)-OuterBeamO31(y);%<=0                                                                %Constraint (5)
                    Constraint6
                    ];
                

end


function CombRound = Round2Res(Combination,Resolution)
% Thsi function discritize the input vector "Combination" according to the
% given "Resolution".
CombRound=nan(size(Combination));
for i =1: length(Combination)
    CombRound(i)=floor(Combination(i)/Resolution(i)+0.5)*Resolution(i);
end
end



















