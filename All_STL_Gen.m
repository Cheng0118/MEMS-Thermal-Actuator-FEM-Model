function [OIB,IIB, OOB,IOB, ORBr,ORBl,IRB, OHLt,OHLb,OHLl,OHLr,IHL, OMLtr,OMLtl,OMLbr,OMLbl,IML, Case]=All_STL_Gen( GeoVariables )
% This function creates the STL files of all parts. All Geometric dimensions has the unit of mm.
% The STL format requires all coordinates to be positive. Therefore we set the origin at a lower left point. 

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


% Outputs
% "O" outputs are outer points of a part. Numbered in counter-clockwise manner.
% "I" outputs are inner points of a part (holes between the beams). Numbered from top holes to bottom holes. The points in each hole are ordered in counter-clockwise manner.
%   OIB: Outer points of InnerBeam. See Page 1 of document.
%   IIB: Inner points of InnerBeam. See Page 1 of document.
%   OOB: Outer points of OuterBeam. See Page 2,3 of document.
%   IOB: Inner points of OuterBeam. See Page 2,3 of document. Order of the holes: [TopRight; TopLeft; BottomRight; BottomLeft]
%   ORBr,ORBl: Outer points of right and left RigidBeam. See Page 4,5 of document.
%   IRM: Empty matrix since there is no inner point in RigidBeam
%   OHLt,OHLb,OHLl,OHLr: Outer points of HandleLayer Top, Bottom, Left and Right. See Page 6,7 of document.
%   IHL: Empty matrix since there is no inner point in RigidBeam
%   OMLtr,OMLtl,OMLbr,OMLbl: Outer points of MetalLayer TopRight, TopLeft, BottomRight, and BottomLeft. See Page 8,9 of document.
%   IML: Empty matrix since there is no inner point in RigidBeam
%   Case: Different Design Cases. Case==1 or 2. See PAge 2,3; 6,7; 8,9 of document for details. 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Begin
[WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables);    % Calculate Basic Geometric Dimensions
rhoSi=0.0002; % Resistivyty of silicon, ohm.m
TD=0.4;     % Thickness of the device
TH=0.05;   % Thickness of the handle layer
TML=0.001;  % Thickness of the metal layer
dOrigin=[3 3 1];    % delta-origin vector

% Calculate InnerBeam Points and create STL
[OIB, IIB, NumberofIBLoop] = InnerBeam(GeoVariables,CW,LM,WA,Ymin,Case);
% CreateSTL(OIB, IIB, NumberofIBLoop, dOrigin, TD, 'InnerBeam.stl', 'ThermalDesign_InnerBeam_Unit[mm]')

% Calculate OuterBeam Points and create STL
[OOB, IOB, NumberofOBLoop] = OuterBeam(GeoVariables,WA, CW, LM, Anchor, Ymax, Ymin, Case);
% CreateSTL(OOB, IOB, NumberofOBLoop, dOrigin, TD, 'OuterBeam.stl', 'ThermalDesign_OuterBeam_Unit[mm]')

% Calculate RigidBeam Points and create STL
[ORBr,ORBl, IRB, NumberofRBLoop] = RigidBeam(GeoVariables, CW, Ld, LM, Ymax, Ymin);
% CreateSTL(ORBr, IRB, NumberofRBLoop, dOrigin, TD, 'RigidBeamR.stl', 'ThermalDesign_RigidBeamR_Unit[mm]')
% CreateSTL(ORBl, IRB, NumberofRBLoop, dOrigin, TD, 'RigidBeamL.stl', 'ThermalDesign_RigidBeamL_Unit[mm]')

% Calculate HandleLayer Points and create STL
[OHLt,OHLb, OHLl,OHLr,IHL, NumberofHLLoop] = HandleLayer(GeoVariables, Case, OIB, OOB, ORBr, ORBl);
% CreateSTL(OHLt, IHL, NumberofHLLoop, dOrigin-[0 0 TH], TH, 'HandleLayerTop.stl', 'ThermalDesign_HandleLayerTop_Unit[mm]')
% CreateSTL(OHLb, IHL, NumberofHLLoop, dOrigin-[0 0 TH], TH, 'HandleLayerBottom.stl', 'ThermalDesign_HandleLayerBottom_Unit[mm]')
% CreateSTL(OHLr, IHL, NumberofHLLoop, dOrigin-[0 0 TH], TH, 'HandleLayerRight.stl', 'ThermalDesign_HandleLayerRight_Unit[mm]')
% CreateSTL(OHLl, IHL, NumberofHLLoop, dOrigin-[0 0 TH], TH, 'HandleLayerLeft.stl', 'ThermalDesign_HandleLayerLeft_Unit[mm]')

% Calculate MetalLayer Points and create STL
[OMLtr,OMLtl, OMLbr,OMLbl,IML, NumberofMLLoop] = MetalLayer(GeoVariables, Case, OIB, OOB);
% CreateSTL(OMLtr, IML, NumberofMLLoop, dOrigin+[0 0 TD], TML, 'MetalLayerTopR.stl', 'ThermalDesign_MetalLayerTopR_Unit[mm]')
% CreateSTL(OMLtl, IML, NumberofMLLoop, dOrigin+[0 0 TD], TML, 'MetalLayerTopL.stl', 'ThermalDesign_MetalLayerTopL_Unit[mm]')
% CreateSTL(OMLbr, IML, NumberofMLLoop, dOrigin+[0 0 TD], TML, 'MetalLayerBottomR.stl', 'ThermalDesign_MetalLayerBottomR_Unit[mm]')
% CreateSTL(OMLbl, IML, NumberofMLLoop, dOrigin+[0 0 TD], TML, 'MetalLayerBottomL.stl', 'ThermalDesign_MetalLayerBottomL_Unit[mm]')

end

function [WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables)
% This function finds a possible combination of CW, WA, and LD

% Outputs:
% WA: Width of side Anchor
% CW: Center Width
% Ld: Width of the Device
% LM: Half width of the middle waist part
% Case: equal to 1 for Case 1, 2 for Case 2
% Ymax: The Y coordinate of highest point.
% Ymin: The Y coordinate of lowest point.

% Every Geo parameter has the unit of mm.
% Example GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02] % Case 1
% Example GeoVariables=[2,0.7,0.04,     10,3,0.4,0.03,     20,3,0.9,0.02] % Case 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;         % Anchor width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;           % Thickness of the device
Wtop=0.05;          % Width of the top plate


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
InnerBeamO1x = Lcb;

% Calculate some basic geometric dimensions
Anchor=max([Anchormin,WT*sin(Atb)+Anchormin*cos(Atb),WB*sin(Abb)+Anchormin*cos(Abb)]);
CW=max([Ltb*cos(Atb)+Anchormin, Lbb*cos(Abb)+Wsc+Gmc,InnerBeamO1x+Anchormin+Gmc]);% Guarantee there is enough space for: 1. Top anchor width, 2. bottom Anchor gap 3. Center beams and anchor gap
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
    
end

function CreateSTL(OuterLoop, InnerLoop, NumberofInnerLoop, dOrigin, Td, Filename, SolidTitle)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The origin is set towards the lower left of the device to make sure that
% all points are positive since STL file requires all points to be positive
% coordinates. 
% dOrigin=[3, 3, 1]; % Origin is (-3, -3)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OuterLoopnew=OuterLoop+ones(size(OuterLoop,1),2)*diag(dOrigin(1:2));
InnerLoopnew=InnerLoop+ones(size(InnerLoop,1),2)*diag(dOrigin(1:2));

% This function create a STL file from a set of OuterLoop and InnerLoop Points.
[AllP, TrueTriList] = TriangleGen(OuterLoopnew, InnerLoopnew, NumberofInnerLoop);
TriVertex=nan(3,2,size(TrueTriList,1));
for i=1:size(TrueTriList,1)
    TriVertex(:,:,i)=AllP(TrueTriList(i,:),:);
end


%Assumptions: *Inner loops are only made up of 4 edges
%*OuterLoop and InnerLoop list points in clockwise order

TriVertexBottom=nan(3,2,size(TrueTriList,1));
TrueTriListBottom=fliplr(TrueTriList); %flip the triangle list to reverse vertex selction (ccw-cw) for the bottom surface

for i=1:size(TrueTriList,1)
    TriVertexBottom(:,:,i)=AllP(TrueTriListBottom(i,:),:);
end
 
%TOP SURFACE FACET MATRIX
zT=dOrigin(3)*ones(3,1,size(TrueTriList,1))+Td;%z position (thickness) column matrix
nT=nan(1,3,size(TrueTriList,1)); %normal vector matrix (top surface)
for i=1:size(TrueTriList,1)
    nT(:,:,i)=[0 0 1];
end
facetsTop=cat(2,TriVertex,zT); %concatenates thickness column to each vertex matrix
facetsTop=cat(1,nT,facetsTop); %this then concatenates the normal row matrix 

%BOTTOM SURFACE FACET MATRIX
zB=dOrigin(3)*ones(3,1,size(TrueTriListBottom,1)); %z position column matrix
nB=nan(1,3,size(TrueTriListBottom,1)); %normal vector matrix (bottom surface)
for i=1:size(TrueTriListBottom,1)
    nB(:,:,i)=[0 0 -1];
end
facetsBottom=cat(2,TriVertexBottom,zB); %concatenates thickness column to each vertex matrix
facetsBottom=cat(1,nB,facetsBottom); %this then concatenates the normal row matrix 

%SIDE SURFACE MATRIX
%NOTE: THE OUTER AND INNER MATRICES MUST LIST POINTS IN THE Counter-CLOCKWISE DIRECTION FOR THIS SECTION OF CODE TO WORK

%create two matrices with the 3d coordiantes of the top and bottom vertices
zTop=dOrigin(3)+Td;
zBottom=dOrigin(3);
OuterPTop=cat(2,OuterLoopnew,zTop*ones(size(OuterLoopnew,1),1)); 
OuterPBottom=cat(2,OuterLoopnew,zBottom*ones(size(OuterLoopnew,1),1));
InnerPTop=cat(2,InnerLoopnew,zTop*ones(size(InnerLoopnew,1),1)); 
InnerPBottom=cat(2,InnerLoopnew,zBottom*ones(size(InnerLoopnew,1),1));

%create triangles bounding the outer sides
OuterSideTriVertex=nan(3,3,2*size(OuterLoopnew,1));
nOuter=nan(1,3,size(OuterLoopnew,1));
for i=1:size(OuterLoopnew,1);      %loops creates 2 triangles per iteration
    if  i==size(OuterLoopnew,1);
        OuterSideTriVertex(:,:,2*i-1)=[OuterPTop(1,:);OuterPTop(i,:);OuterPBottom(i,:)];
        OuterSideTriVertex(:,:,2*i)=[OuterPBottom(i,:);OuterPBottom(1,:);OuterPTop(1,:)];
        Crossvec=cross(OuterPTop(i,:)-OuterPTop(1,:),OuterPBottom(i,:)-OuterPTop(i,:));
    else
        OuterSideTriVertex(:,:,2*i-1)=[OuterPTop(i+1,:);OuterPTop(i,:);OuterPBottom(i,:)];
        OuterSideTriVertex(:,:,2*i)=[OuterPBottom(i,:);OuterPBottom(i+1,:);OuterPTop(i+1,:)];
        Crossvec=cross(OuterPTop(i,:)-OuterPTop(i+1,:),OuterPBottom(i,:)-OuterPTop(i,:));
    end
    nOuter(1,:,2*i-1)=Crossvec/norm(Crossvec);
    nOuter(1,:,2*i)=Crossvec/norm(Crossvec);
end

%create triangles bounding the inner sides
InnerSideTriVertex=nan(3,3,2*size(InnerLoopnew,1));
nInner=nan(1,3,size(InnerLoopnew,1));
for iNI=0:NumberofInnerLoop-1
    for i=iNI*4+1:iNI*4+4;      %loops creates 2 triangles per iteration
        if i==iNI*4+4;
            InnerSideTriVertex(:,:,2*i-1)=[InnerPTop(i,:);InnerPTop(iNI*4+1,:);InnerPBottom(iNI*4+1,:)];
            InnerSideTriVertex(:,:,2*i)=[InnerPBottom(iNI*4+1,:);InnerPBottom(i,:);InnerPTop(i,:)];
            Crossvec=cross(InnerPTop(iNI*4+1,:)-InnerPTop(i,:),InnerPBottom(iNI*4+1,:)-InnerPTop(iNI*4+1,:));
        else
            InnerSideTriVertex(:,:,2*i-1)=[InnerPTop(i,:);InnerPTop(i+1,:);InnerPBottom(i+1,:)];
            InnerSideTriVertex(:,:,2*i)=[InnerPBottom(i+1,:);InnerPBottom(i,:);InnerPTop(i,:)];
        	Crossvec=cross(InnerPTop(i+1,:)-InnerPTop(i,:),InnerPBottom(i+1,:)-InnerPTop(i+1,:)); 
        end
        nInner(1,:,2*i-1)=Crossvec/norm(Crossvec);
        nInner(1,:,2*i)=Crossvec/norm(Crossvec);
    end
end

%concatenate side normal vectors to the side trivertex matrices
facetsOuter=cat(1,nOuter,OuterSideTriVertex); %this then concatenates the normal row matrix 
facetsInner=cat(1,nInner,InnerSideTriVertex); %this then concatenates the normal row matrix 

%COMBINED FACET MATRIX
facets43=cat(3,facetsTop,facetsBottom,facetsOuter,facetsInner);
facets34=permute(facets43,[2 1 3]);
facets=round(facets34,6);






STLAcsii = fopen(Filename, 'w');
% Write HEADER
fprintf(STLAcsii,'solid %s\r\n',SolidTitle);
% Write DATA
fprintf(STLAcsii,[...
    '   facet normal %.6E %.6E %.6E\r\n' ...
    '      outer loop\r\n' ...
    '         vertex %.6E %.6E %.6E\r\n' ...
    '         vertex %.6E %.6E %.6E\r\n' ...
    '         vertex %.6E %.6E %.6E\r\n' ...
    '      endloop\r\n' ...
    '   endfacet\r\n'], facets);
% Write FOOTER
fprintf(STLAcsii,'endsolid %s\r\n',SolidTitle);
% Close the file
fclose(STLAcsii);


end

function [AllP, TrueTri] = TriangleGen(OuterLoop, InnerLoop, NumberofInnerLoop)
% This function triangulates the shape specified in OuterLoop and InnerLoop.

% The input OuterLoop and InnerLoop needs to be a 2 column matrix. No repeated points are allowed.
% (The end point should be the last point before connecting back to start point.)

% The output PointsList has 3 columns specifying the coordinates of vertex counter-clockwisely. 
AllPoints=[OuterLoop;InnerLoop];
Nout=size(OuterLoop,1);
if NumberofInnerLoop==0
    Nin=4;
else
    Nin=size(InnerLoop,1)/NumberofInnerLoop;
end
if Nin~=4
    disp('Warning: check inner loop')
end
Outeredge=[(1:Nout)',[2:Nout,1]'];
Inneredge=nan(NumberofInnerLoop*Nin,2);

for i=1:NumberofInnerLoop
    Inneredge(((i-1)*Nin+1):(i*Nin),1:2)=[((Nout+1:Nout+Nin)+(i-1)*Nin)',([Nout+2:Nout+Nin,Nout+1]+(i-1)*Nin)'];
end
C=[Outeredge;Inneredge];
DelauTri = delaunayTriangulation(AllPoints,C);
AllTriIndex=DelauTri.ConnectivityList;
AllP=DelauTri.Points;
TrueTri = AllTriIndex(isInterior(DelauTri),:);


%  triplot(TrueTri,AllP(:,1),AllP(:,2))
end

function [OIB, IIB, NumberofIBLoop] = InnerBeam(GeoVariables,CW,LM,WA,Ymin,Case)
% This file create the list of the points of the inner/center beam for the creation of an STL file.

% Every Geo parameter has the unit of mm.
% Example GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02], [WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;      % Anchor min width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;           % Thickness of the device
Wtop=0.05;          % Width of the top plate

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


WB=Nbb*Wbb+(Nbb-1)*Gbb;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The Outer "O" points are numbered starting from the highest right point
% counter-clockwisely. The initial origin is set at the center point of the
% center beams
OIB=nan(24,2);
x=1;y=2;            % For a clearer notation we define x and y.

i=1;  OIB(i,:)=[Lcb,(Ncb*Wcb+(Ncb-1)*Gcb)/2];                                           OIB(2,:)=[-OIB(i,x),OIB(i,y)];
% To avoid typos we define i. Each pair of O is symmetric about the center axis.
i=14; OIB(i,:)=[Lcb,-(Ncb*Wcb+(Ncb-1)*Gcb)/2];                                          OIB(13,:)=[-OIB(i,x),OIB(i,y)];
i=15; OIB(i,:)=[OIB(14,x),OIB(14,y)-Gmg];                                               OIB(12,:)=[-OIB(i,x),OIB(i,y)];
i=24; OIB(i,:)=[CW-Gmc,WA/2];                                                           OIB(3,:)=[-OIB(i,x),OIB(i,y)];
i=23; OIB(i,:)=[CW-Gmc,-WA/2];                                                          OIB(4,:)=[-OIB(i,x),OIB(i,y)];
i=16; OIB(i,:)=[Gmc/2,OIB(15,y)];                                                       OIB(11,:)=[-OIB(i,x),OIB(i,y)];
i=21; OIB(i,:)=[OIB(16,x)+Wmc,OIB(16,y)-Wmc];                                           OIB(6,:)=[-OIB(i,x),OIB(i,y)];
i=22; OIB(i,:)=[OIB(15,x),OIB(21,y)];                                                   OIB(5,:)=[-OIB(i,x),OIB(i,y)];
i=17; OIB(i,:)=[OIB(16,x),Ymin+Wsc+Gmc];                                                OIB(10,:)=[-OIB(i,x),OIB(i,y)];   
i=18; OIB(i,:)=[(LM-Wsc-Gmc-Gmc/2)/2,OIB(17,y)];                                        OIB(9,:)=[-OIB(i,x),OIB(i,y)];  

if Case==1
    OuterBeamO23y=-WA/2-Lbb*sin(Abb)-WB*cos(Abb)-Gbb-Anchormin;
    i=19; OIB(i,:)=[OIB(18,x),OuterBeamO23y-Gmc];                                       OIB(8,:)=[-OIB(i,x),OIB(i,y)];
    i=20; OIB(i,:)=[OIB(21,x),OIB(19,y)];                                               OIB(7,:)=[-OIB(i,x),OIB(i,y)];

else
    i=19; OIB(i,:)=[OIB(18,x),-WA/2-(CW-OIB(18,x))*tan(Abb)];                           OIB(8,:)=[-OIB(i,x),OIB(i,y)];
    i=20; OIB(i,:)=[OIB(21,x),OIB(19,y)];                                               OIB(7,:)=[-OIB(i,x),OIB(i,y)];
end



% The Inner "I" points are numbered starting from the highest right hole, and numbered counter-clockwisely for each hole.
IIB=nan((Ncb-1)*4,2);
for i=1:Ncb-1
    IIB(4*i-3,:)=[OIB(1,x),OIB(1,y)-i*(Wcb+Gcb)+Gcb];
    IIB(4*i-2,:)=[OIB(2,x),OIB(2,y)-i*(Wcb+Gcb)+Gcb];
    IIB(4*i-1,:)=[OIB(2,x),OIB(2,y)-i*(Wcb+Gcb)];
    IIB(4*i,:) = [OIB(1,x),OIB(1,y)-i*(Wcb+Gcb)];
end

NumberofIBLoop=Ncb-1;

end

function [OOB, IOB, NumberofOBLoop] = OuterBeam(GeoVariables,WA, CW, LM, Anchor, Ymax, Ymin, Case)
% This file create the list of the points of the inner/center beam for the creation of an STL file.

% Every Geo parameter has the unit of mm.
% Example GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02], [WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;      % Anchor width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;           % Thickness of the device
Wtop=0.05;          % Width of the top plate

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% The Outer "O" points are numbered starting from the highest right point
% counter-clockwisely. The initial origin is set at the center point 1of the
% center beams
OOB=nan(32,2);
x=1;y=2;            % For a clearer notation we define x and y.

i=19;  OOB(i,:)=[CW,WA/2];                                                          OOB(16,:)=[-OOB(i,x),OOB(i,y)];
% To avoid typos we define i. Each pair of O is symmetric about the center axis.
i=18;  OOB(i,:)=[CW-Ltb*cos(Atb),WA/2+Ltb*sin(Atb)];                                OOB(17,:)=[-OOB(i,x),OOB(i,y)];
i=32;  OOB(i,:)=[OOB(18,x)+WT*sin(Atb),OOB(18,y)+WT*cos(Atb)];                  	OOB(3,:)=[-OOB(i,x),OOB(i,y)];
i=1;   OOB(i,:)=[OOB(32,x),Ymax-Gmc];                                               OOB(2,:)=[-OOB(i,x),OOB(i,y)];
i=31;  OOB(i,:)=[CW+Anchor,OOB(32,y)-(CW+Anchor-OOB(32,x))*tan(Atb)];            	OOB(4,:)=[-OOB(i,x),OOB(i,y)];
i=20;  OOB(i,:)=[CW,-WA/2];                                                      	OOB(15,:)=[-OOB(i,x),OOB(i,y)];
i=21;  OOB(i,:)=[CW-Lbb*cos(Abb),-WA/2-Lbb*sin(Abb)];                               OOB(14,:)=[-OOB(i,x),OOB(i,y)];
i=29;  OOB(i,:)=[OOB(21,x)+WB*sin(Abb),OOB(21,y)-WB*cos(Abb)];                    	OOB(6,:)=[-OOB(i,x),OOB(i,y)];
i=28;  OOB(i,:)=[OOB(29,x),OOB(29,y)-Gbb];                                          OOB(7,:)=[-OOB(i,x),OOB(i,y)];
InnerBeamO19x=(LM-Wsc-Gmc-Gmc/2)/2;
i=25;  OOB(i,:)=[InnerBeamO19x+Gmc,Ymin+Wsc+Gmc];                                   OOB(10,:)=[-OOB(i,x),OOB(i,y)];
i=26;  OOB(i,:)=[LM-Gmc-Wsc,Ymin+Wsc+Gmc];                                          OOB(9,:)=[-OOB(i,x),OOB(i,y)];
i=30;  OOB(i,:)=[OOB(31,x),OOB(29,y)+(CW+Anchor-OOB(29,x))*tan(Abb)];            	OOB(5,:)=[-OOB(i,x),OOB(i,y)];

if Case==1
    i=22;  OOB(i,:)=[CW-(Lbb+Anchormin)*cos(Abb),-WA/2-(Lbb+Anchormin)*sin(Abb)];  	OOB(13,:)=[-OOB(i,x),OOB(i,y)];
    i=23;  OOB(i,:)=[OOB(22,x),OOB(28,y)-Anchormin];                                  	OOB(12,:)=[-OOB(i,x),OOB(i,y)];
    i=27;  OOB(i,:)=[OOB(26,x),OOB(28,y)];                                          OOB(8,:)=[-OOB(i,x),OOB(i,y)];
    i=24;  OOB(i,:)=[InnerBeamO19x+Gmc,OOB(23,y)];                                 	OOB(11,:)=[-OOB(i,x),OOB(i,y)];
else
    i=24;  OOB(i,:)=[InnerBeamO19x+Gmc,OOB(21,y)-(OOB(21,x)-InnerBeamO19x-Gmc)*tan(Abb)];     OOB(11,:)=[-OOB(i,x),OOB(i,y)];
    i=22;  OOB(i,:)=[2/3*OOB(21,x)+1/3*OOB(24,x),2/3*OOB(21,y)+1/3*OOB(24,y)];           	OOB(13,:)=[-OOB(i,x),OOB(i,y)];
    i=23;  OOB(i,:)=[1/3*OOB(21,x)+2/3*OOB(24,x),1/3*OOB(21,y)+2/3*OOB(24,y)];         	OOB(12,:)=[-OOB(i,x),OOB(i,y)];
    i=27;  OOB(i,:)=[OOB(26,x),OOB(28,y)-(OOB(28,x)-OOB(26,x))*tan(Abb)];               	OOB(8,:)=[-OOB(i,x),OOB(i,y)];
end



% The Inner "I" points are numbered starting from the highest right hole, and numbered counter-clockwisely for each hole.
Itr=nan((Ntb-1)*4,2);   % Top right holes
for i=1:Ntb-1
    Itr(4*i-3,:)=[OOB(32,x)-(i*(Wtb+Gtb)-Gtb)*sin(Atb),                 OOB(32,y)-(i*(Wtb+Gtb)-Gtb)*cos(Atb)];
    Itr(4*i-2,:)=[OOB(32,x)-i*(Wtb+Gtb)*sin(Atb),                       OOB(32,y)-i*(Wtb+Gtb)*cos(Atb)];
    Itr(4*i-1,:)=[OOB(32,x)-i*(Wtb+Gtb)*sin(Atb)+Ltb*cos(Atb),          OOB(32,y)-i*(Wtb+Gtb)*cos(Atb)-Ltb*sin(Atb)];
    Itr(4*i,:) = [OOB(32,x)-(i*(Wtb+Gtb)-Gtb)*sin(Atb)+Ltb*cos(Atb),    OOB(32,y)-(i*(Wtb+Gtb)-Gtb)*cos(Atb)-Ltb*sin(Atb)];
end

Itl=nan((Ntb-1)*4,2);   % Top left holes
for i=1:Ntb-1
    Itl(4*i-3,:)=[OOB(3,x)+i*(Wtb+Gtb)*sin(Atb),                        OOB(3,y)-i*(Wtb+Gtb)*cos(Atb)];
    Itl(4*i-2,:)=[OOB(3,x)+(i*(Wtb+Gtb)-Gtb)*sin(Atb),                  OOB(3,y)-(i*(Wtb+Gtb)-Gtb)*cos(Atb)];
    Itl(4*i-1,:)=[OOB(3,x)+(i*(Wtb+Gtb)-Gtb)*sin(Atb)-Ltb*cos(Atb),     OOB(3,y)-(i*(Wtb+Gtb)-Gtb)*cos(Atb)-Ltb*sin(Atb)];
    Itl(4*i,:) = [OOB(3,x)+i*(Wtb+Gtb)*sin(Atb)-Ltb*cos(Atb),           OOB(3,y)-i*(Wtb+Gtb)*cos(Atb)-Ltb*sin(Atb)];
end

Ilr=nan((Nbb-1)*4,2);   % Lower right holes
for i=1:Nbb-1
    Ilr(4*i-3,:)=[OOB(20,x)+i*(Wbb+Gbb)*sin(Abb),                         OOB(20,y)-i*(Wbb+Gbb)*cos(Abb)];
    Ilr(4*i-2,:)=[OOB(20,x)+(i*(Wbb+Gbb)-Gbb)*sin(Abb),                   OOB(20,y)-(i*(Wbb+Gbb)-Gbb)*cos(Abb)];
    Ilr(4*i-1,:)=[OOB(21,x)+(i*(Wbb+Gbb)-Gbb)*sin(Abb),                   OOB(21,y)-(i*(Wbb+Gbb)-Gbb)*cos(Abb)];
    Ilr(4*i,:) = [OOB(21,x)+i*(Wbb+Gbb)*sin(Abb),                         OOB(21,y)-i*(Wbb+Gbb)*cos(Abb)];

end

Ill=nan((Nbb-1)*4,2);   % Lower left holes
for i=1:Nbb-1
    Ill(4*i-3,:)=[OOB(15,x)-(i*(Wbb+Gbb)-Gbb)*sin(Abb),                   OOB(15,y)-(i*(Wbb+Gbb)-Gbb)*cos(Abb)];
    Ill(4*i-2,:)=[OOB(15,x)-i*(Wbb+Gbb)*sin(Abb),                         OOB(15,y)-i*(Wbb+Gbb)*cos(Abb)];
    Ill(4*i-1,:)=[OOB(14,x)-i*(Wbb+Gbb)*sin(Abb),                         OOB(14,y)-i*(Wbb+Gbb)*cos(Abb)];
    Ill(4*i,:) = [OOB(14,x)-(i*(Wbb+Gbb)-Gbb)*sin(Abb),                   OOB(14,y)-(i*(Wbb+Gbb)-Gbb)*cos(Abb)];
end

IOB=[Itr;Itl;Ilr;Ill];

NumberofOBLoop=(Ntb-1)*2+(Nbb-1)*2;

end

function [ORBr,ORBl, IRB, NumberofRBLoop] = RigidBeam(GeoVariables, CW, Ld, LM, Ymax, Ymin)
% This file create the list of the points of the inner/center beam for the creation of an STL file.

% Every Geo parameter has the unit of mm.
% Example GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02], [WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;    	% Anchor min width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;           % Thickness of the device
Wtop=0.05;          % Width of the top plate

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% The Outer "O" points are numbered starting from the highest right point
% counter-clockwisely. The initial origin is set at the center point of the
% center beams
ORBr=nan(16,2);
x=1;y=2;            % For a clearer notation we define x and y.

ORBr(1,:)=[Ld,Ymax];                                                        
% To avoid typos we define i. Each pair of O is symmetric about the center axis.
ORBr(16,:)=[Ld,Ymax-Wtop];
ORBr(2,:)=[Gmc+(CW-Ltb*cos(Atb)+WT*sin(Atb)),Ymax];
ORBr(15,:)=[ORBr(2,x)+Wmc,ORBr(16,y)];
ORBr(14,:)=[ORBr(15,x),ORBr(15,y)-2*Gmg];
ORBr(13,:)=[Ld,ORBr(14,y)];
ORBr(3,:)=[ORBr(2,x),ORBr(14,y)-Wmc];
ORBr(4,:)=[Ld-Wsb,ORBr(3,y)];
ORBr(12,:)=[Ld,Ymin+Ld+Td*tan(Abb)];
ORBr(5,:)=[ORBr(12,x)-Wsb,ORBr(12,y)+Wsc*cos(Abb)-(Wsb-Wsc*sin(Abb))*tan(Abb)];
ORBr(11,:)=[Ld-Td,Ymin+Ld];
ORBr(6,:)=[Ld-Td-Wsc,ORBr(5,y)-(Td+Wsc-Wsb)*tan(Abb)];
ORBr(10,:)=[LM,Ymin];
ORBr(9,:)=[Gmc/2,Ymin];
ORBr(8,:)=[Gmc/2,Ymin+Wsc];
ORBr(7,:)=[LM-Wsc,Ymin+Wsc];

ORBl=flipud([-ORBr(:,x),ORBr(:,y)]);

IRB=ones(size([],1),2);
NumberofRBLoop=0;

end

function [OHLt,OHLb, OHLl,OHLr,IHL, NumberofHLLoop] = HandleLayer(GeoVariables, Case, OIB, OOB, ORBr, ORBl)
% This file create the list of the points of handle layer structure for the creation of an STL file.

% Every Geo parameter has the unit of mm.
% Example
% GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02],
% [WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables)
% [ORBr,ORBl, IRB, NumberofRBLoop] = RigidBeam(GeoVariables, CW, Ld, LM, Ymax, Ymin);
% [OOB, IOB, NumberofOBLoop] = OuterBeam(GeoVariables,WA, CW, LM, Ymax, Ymin, Case);
% [OIB, IIB, NumberofIBLoop] = InnerBeam(GeoVariables,CW,LM,WA,Ymin,Case);
% [OHLt,OHLb, OHLl,OHLr,IHL, NumberofHLLoop] = HandleLayer(GeoVariables, OIB, OOB, ORBr, ORBl)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;     	% Anchor min width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;           % Thickness of the device
Wtop=0.05;          % Width of the top plate

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% The Outer "O" points are numbered starting from the highest right point
% counter-clockwisely. The initial origin is set at the center point of the
% center beams


x=1;y=2;            % For a clearer notation we define x and y.
% Set OHLt
i=1; OHLt(i,:)=ORBr(1,:);               	OHLt(2,:)=[-OHLt(i,x),OHLt(i,y)];
i=3; OHLt(i,:)=ORBl(1,:);                   OHLt(12,:)=[-OHLt(i,x),OHLt(i,y)];
i=4; OHLt(i,:)=ORBl(2,:);                   OHLt(11,:)=[-OHLt(i,x),OHLt(i,y)];
i=5; OHLt(i,:)=[ORBl(3,x),ORBl(14,y)];      OHLt(10,:)=[-OHLt(i,x),OHLt(i,y)];
i=6; OHLt(i,:)=OOB(3,:);                    OHLt(9,:)=[-OHLt(i,x),OHLt(i,y)];
i=7; OHLt(i,:)=OOB(17,:);                   OHLt(8,:)=[-OHLt(i,x),OHLt(i,y)];

% Set OHLr
OHLr(2,:)=OOB(19,:);
OHLr(1,:)=[OHLr(2,x)+WT*sin(Atb),OHLr(2,y)+WT*cos(Atb)];
OHLr(3,:)=OIB(24,:);
OHLr(4,:)=OIB(1,:);
OHLr(5,:)=OIB(22,:);
OHLr(6,:)=OIB(23,:);
OHLr(7,:)=OOB(20,:);
OHLr(8,:)=[OHLr(7,x)+WB*sin(Abb),OHLr(7,y)-WB*cos(Abb)];
OHLr(9,:)=OOB(30,:);
OHLr(10,:)=OOB(31,:);

% Set OHLl
OHLl=flipud([-OHLr(:,x),OHLr(:,y)]);

% Set OHLb
i=19; OHLb(i,:)=OOB(27,:);              	OHLb(4,:)=[-OHLb(i,x),OHLb(i,y)];
i=18; OHLb(i,:)=ORBr(6,:);                  OHLb(5,:)=[-OHLb(i,x),OHLb(i,y)];
i=17; OHLb(i,:)=ORBr(5,:);                  OHLb(6,:)=[-OHLb(i,x),OHLb(i,y)];
i=16; OHLb(i,:)=[ORBr(4,x),ORBr(13,y)];    	OHLb(7,:)=[-OHLb(i,x),OHLb(i,y)];
i=15; OHLb(i,:)=ORBr(13,:);                 OHLb(8,:)=[-OHLb(i,x),OHLb(i,y)];
i=14; OHLb(i,:)=ORBr(12,:);                 OHLb(9,:)=[-OHLb(i,x),OHLb(i,y)];
i=13; OHLb(i,:)=ORBr(11,:);                 OHLb(10,:)=[-OHLb(i,x),OHLb(i,y)];
i=12; OHLb(i,:)=ORBr(10,:);                 OHLb(11,:)=[-OHLb(i,x),OHLb(i,y)];

if Case==1
    i=20; OHLb(i,:)=OOB(28,:);             	OHLb(3,:)=[-OHLb(i,x),OHLb(i,y)];
    i=1;  OHLb(i,:)=[OOB(23,x),OOB(28,y)]; 	OHLb(2,:)=[-OHLb(i,x),OHLb(i,y)];
else
    i=20; OHLb(i,:)=OOB(22,:);             	OHLb(3,:)=[-OHLb(i,x),OHLb(i,y)];
    i=1;  OHLb(i,:)=OIB(19,:);              OHLb(2,:)=[-OHLb(i,x),OHLb(i,y)]; 
end




IHL=ones(size([],1),2);
NumberofHLLoop=0;

end

function [OMLtr,OMLtl, OMLbr,OMLbl,IML, NumberofMLLoop] = MetalLayer(GeoVariables, Case, OIB, OOB)
% This file create the list of the points of handle layer structure for the creation of an STL file.

% Every Geo parameter has the unit of mm.
% Example
% GeoVariables=[2,0.7,0.04,     10,3,0.7,0.03,     20,3,0.9,0.02],
% [WA, CW, Ld, LM, Anchor, Ymax, Ymin, Case] = BasicGeoCalc(GeoVariables)
% [ORBr,ORBl, IRB, NumberofRBLoop] = RigidBeam(GeoVariables, CW, Ld, LM, Ymax, Ymin);
% [OOB, IOB, NumberofOBLoop] = OuterBeam(GeoVariables,WA, CW, LM, Ymax, Ymin, Case);
% [OIB, IIB, NumberofIBLoop] = InnerBeam(GeoVariables,CW,LM,WA,Ymin,Case);
% [OMLtr,OMLtl, OMLbr,OMLbl,IML, NumberofMLLoop] = MetalLayer(GeoVariables, Case, OIB, OOB)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize all geometric variables.
Gcb=0.01;           % Center Beam Gap
Gbb=0.01;           % Bottom Beam Gap
Gtb=0.01;           % Top Beam Gap
Gmg=0.04;           % "Moving Gap"
Anchormin=0.1;      % Anchor min width
Wmc=0.03;           % Width of middle connection
Gmc=0.05;           % Gap of middle connection
HH=2.5;             % Half height of the cube
Wsb=0.15;           % Width of the side beam(rigid upper part).
Wsc=0.2;            % Width of the side beam(lower part).
Td=0.45;           % Thickness of the device
Wtop=0.05;          % Width of the top plate

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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




% The Outer "O" points are numbered starting from the highest right point
% counter-clockwisely. The initial origin is set at the center point of the
% center beams


x=1;y=2;            % For a clearer notation we define x and y.

% Set OMLtr
OMLtr(1,:)=[OIB(24,x)-0.005,OIB(24,y)-0.01];
OMLtr(2,:)=[OIB(1,x)+0.005,OIB(1,y)-0.005];
OMLtr(3,:)=[OIB(14,x)+0.005,OIB(14,y)+0.005];
OMLtr(4,:)=[OIB(23,x)-0.005,OIB(23,y)+0.01];

% Set OMLtl
OMLtl=flipud([-OMLtr(:,x),OMLtr(:,y)]);

% Set OMLbr
OMLbr(1,:)=[OOB(21,x)+0.005*(sin(Abb)-cos(Abb)),OOB(21,y)-0.005*(sin(Abb)+cos(Abb))];
OMLbr(4,:)=[OOB(29,x)-0.005*(sin(Abb)+cos(Abb)),OOB(29,y)-0.005*(sin(Abb)-cos(Abb))];

if Case==1
    OMLbr(2,:)=[OOB(22,x)+0.005*(sin(Abb)+cos(Abb)),OOB(22,y)+0.005*(sin(Abb)-cos(Abb))];
else
    OMLbr(2,:)=[OOB(24,x)+0.005*(sin(Abb)+cos(Abb)),OOB(24,y)+0.005*(sin(Abb)-cos(Abb))];
end
OMLbr(3,:)=OMLbr(2,:)+OMLbr(4,:)-OMLbr(1,:);

% Set OMLbl
OMLbl=flipud([-OMLbr(:,x),OMLbr(:,y)]);


IML=ones(size([],1),2);
NumberofMLLoop=0;

end


