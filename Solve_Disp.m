%% Thermal
% Temperature Profile

T(:,1) = Results_E;
T(:,2) = Results_C;

Temperature = zeros(num_links,1);

for ii = 1:2
    
% Outer
for i = 1:N1-1
    Temperature(i,ii) = T(2,ii);
end
for i = 1:2
    Temperature((N1-1)+i,ii) = T(12,ii);
end
for i = 1:N3-1
    Temperature((N1-1)+2+i,ii) = T(4,ii);
end
for i = 1:N1
    Temperature((N1-1)+2+(N3-1)+i,ii) = T(3,ii);
end
for i = 1:N3
    Temperature((N1-1)+2+(N3-1)+N1+i,ii) = T(6,ii);
end
for i = 1:N1-1
    Temperature((N1-1)+2+(N3-1)+N1+N3+i,ii) = T(1,ii);
end
for i = 1:1
    Temperature((N1-1)+2+(N3-1)+N1+N3+(N1-1)+i,ii) = T(11,ii);
end
for i = 1:1
    Temperature((N1-1)+2+(N3-1)+N1+N3+(N1-1)+1+i,ii) = T(10,ii);
end

% Inner
if mod(N2,2) == 0
   for i = 1:N2
       Temperature(num_links_outer+i,ii) = T(13,ii);
   end
   for i = 1:N2
       Temperature(num_links_outer+N2+i,ii) = T(7,ii);
   end
   for i = 1:1
       Temperature(num_links_outer+N2+N2+i,ii) = T(13,ii);
   end
   for i = 1:1
       Temperature(num_links_outer+N2+N2+1+i,ii) = T(8,ii);
   end
   for i = 1:1
       Temperature(num_links_outer+N2+N2+1+1+i,ii) = T(9,ii);
   end   
else
   for i = 1:(N2-1)
       Temperature(num_links_outer+i,ii) = T(13,ii);
   end
   for i = 1:N2
       Temperature(num_links_outer+(N2-1)+i,ii) = T(7,ii);
   end
   for i = 1:1
       Temperature(num_links_outer+(N2-1)+N2+i,ii) = T(13,ii);
   end
   for i = 1:1
       Temperature(num_links_outer+(N2-1)+N2+1+i,ii) = T(8,ii);
   end
   for i = 1:1
       Temperature(num_links_outer+(N2-1)+N2+1+1+i,ii) = T(9,ii);
   end
end

% Anchor
for i = 1:1
    Temperature(num_links_outer+num_links_inner+i,ii) = 23;
end
for i = 1:1
    Temperature(num_links_outer+num_links_inner+1+i,ii) = T(14,ii);
end
for i = 1:1
    Temperature(num_links_outer+num_links_inner+2+i,ii) = 23;
end
for i = 1:3
    Temperature(num_links_outer+num_links_inner+3+i,ii) = T(13,ii)/2;
end

end

% Temperature
dT = zeros(num_links,1);
T_K = zeros(num_links,1);
alpha = zeros(num_links,1);
for ii = 1:2
for i = 1:num_links
    dT(i,ii) = Temperature(i,ii)-23;
    T_K(i,ii) = Temperature(i,ii)+273;
    alpha(i,ii) = (3.725*(1-exp(-5.88*10^(-3)*(T_K(i,ii)-124)))+5.548*10^(-4)*T_K(i,ii))*10^(-6);
end
end

%% Loads
ff_t = zeros(6*num_nodes,2);
for ii = 1:2
    for i = 1:num_links
        SN = Links(i,1); EN = Links(i,2);
        if SN == N1+1 || SN == num_nodes_outer+num_nodes_inner+2
        else
            ff_t((SN-1)*6+1,ii) = ff_t((SN-1)*6+1,ii) - (Coordinates(EN,1)-Coordinates(SN,1))/sqrt((Coordinates(EN,1)-Coordinates(SN,1))^2+(Coordinates(EN,2)-Coordinates(SN,2))^2)*A(i,1)*E*alpha(i,ii)*dT(i,ii);
            ff_t((SN-1)*6+2,ii) = ff_t((SN-1)*6+2,ii) - (Coordinates(EN,2)-Coordinates(SN,2))/sqrt((Coordinates(EN,1)-Coordinates(SN,1))^2+(Coordinates(EN,2)-Coordinates(SN,2))^2)*A(i,1)*E*alpha(i,ii)*dT(i,ii);
            ff_t((EN-1)*6+1,ii) = ff_t((EN-1)*6+1,ii) + (Coordinates(EN,1)-Coordinates(SN,1))/sqrt((Coordinates(EN,1)-Coordinates(SN,1))^2+(Coordinates(EN,2)-Coordinates(SN,2))^2)*A(i,1)*E*alpha(i,ii)*dT(i,ii);
            ff_t((EN-1)*6+2,ii) = ff_t((EN-1)*6+2,ii) + (Coordinates(EN,2)-Coordinates(SN,2))/sqrt((Coordinates(EN,1)-Coordinates(SN,1))^2+(Coordinates(EN,2)-Coordinates(SN,2))^2)*A(i,1)*E*alpha(i,ii)*dT(i,ii);
        end
    end
end

N_3 = 1;
for i = 1:num_nodes
    for ii = 1:6
        if DOFs(i,ii) == 0
           ff_t(N_3,:) = [];
        else
           N_3 = N_3+1;
        end
    end
end


%% Solve
Displacement = zeros(numel(ff_t(:,1)),2);
Disp_Out = zeros(numel(ff_t(:,1)),2);

for ii = 1:2
Displacement(:,ii) = K_ff\ff_t(:,ii);
Disp_Out(:,ii) = Displacement(:,ii);
end

Outputs = zeros(6*num_nodes,1);
Outputs(6*(num_nodes_outer)-4,1) = 1;

N_4 = 1;
for i = 1:num_nodes
    for ii = 1:6
        if DOFs(i,ii) == 0
            Outputs(N_4,:) = [];
        else
            N_4 = N_4+1;
        end
    end
end

N_5 = 1;
for i = 1:size(Outputs)
    if Outputs(i,1) == 0
        Disp_Out(N_5,:) = [];
    else
        N_5 = N_5+1;
    end
end

Stroke = abs(Disp_Out(1)-Disp_Out(2));