%% Compression
% Links
Links_C = zeros(num_links_inner,2);
Links_C(:,:) = Links(num_links_outer+1:num_links_outer+num_links_inner,:);

% Lengths
L_C = L(num_links_outer+1:num_links_outer+num_links_inner,1);

% Cross-Section Areas
A_C = A(num_links_outer+1:num_links_outer+num_links_inner,1);

% Circuit
V_node = zeros(2,2);
V_node(1,1) = 17;
V_node(1,2) = Vin_C;
V_node(2,1) = 21;
V_node(2,2) = 0;

%% Solving Circuit
% Resistance
R = zeros(num_links_inner,1);
for i = 1:num_links_inner
    R(i,1) = row*L_C(i,1)/A_C(i,1);
end

Links_R = zeros(num_links_inner,3);
Links_R(:,1:2) = Links_C(:,1:2);
Links_R(:,3) = R;

AA = zeros(num_nodes_inner,num_nodes_inner);
BB = zeros(num_nodes_inner,1);

% AA Matrix Creation
for i = 1:num_links_inner
    n1 = Links_R(i,1)-num_nodes_outer;
    n2 = Links_R(i,2)-num_nodes_outer;
    if n1 == n2
    else
        AA(n1,n2) = AA(n1,n2) - 1/Links_R(i,3);
        AA(n2,n1) = AA(n2,n1) - 1/Links_R(i,3);
        AA(n1,n1) = AA(n1,n1) + 1/Links_R(i,3);
        AA(n2,n2) = AA(n2,n2) + 1/Links_R(i,3);
    end
end

% BB Matrix Creation 
for i = 1:size(V_node,1)
    AA(V_node(i,1)-num_nodes_outer,:) = zeros(1,num_nodes_inner);
    AA(V_node(i,1)-num_nodes_outer,V_node(i,1)-num_nodes_outer) = 1;
    BB(V_node(i,1)-num_nodes_outer,1) = V_node(i,2);
end

Vo = AA\BB;

Links_V = zeros(num_links_inner,1);
for i = 1:num_links_inner
    n1 = Links_R(i,1)-num_nodes_outer;
    n2 = Links_R(i,2)-num_nodes_outer;
    Links_V(i,1) = Vo(n1,1)-Vo(n2,1);
end

Links_P = Links_V.*Links_V./R;

%% New Full Structure
num_links_N = (n_n+1)*num_links;
num_nodes_N = num_nodes+n_n*num_links;

% Links
New_nodes = zeros(n_n*num_links,1);
New_nodes(:,1) = num_nodes+1:num_nodes_N;

New_links = zeros(num_links_N,3);
New_links(:,1) = 1:num_links_N;

for i = 1:num_links
    for ii = 1:n_n+1
        if ii == 1
            New_links((i-1)*(n_n+1)+1,2) = Links(i,1);
            New_links((i-1)*(n_n+1)+ii,3) = New_nodes((i-1)*n_n+ii,1);
        elseif ii == n_n+1
            New_links((i-1)*(n_n+1)+ii,2) = New_nodes((i-1)*n_n+ii-1,1);            
            New_links(i*(n_n+1),3) = Links(i,2);
        else
            New_links((i-1)*(n_n+1)+ii,2) = New_nodes((i-1)*n_n+ii-1,1);
            New_links((i-1)*(n_n+1)+ii,3) = New_nodes((i-1)*n_n+ii,1);
        end
    end
end
    
% Lengths
New_L = zeros(num_links_N,1);
for i = 1:num_links
    New_L((i-1)*(n_n+1)+1:i*(n_n+1)) = L(i)/(n_n+1);
end

% Cross-Section Areas
New_A = zeros(num_links_N,1);
for i = 1:num_links
    New_A((i-1)*(n_n+1)+1:i*(n_n+1)) = A(i);
end

% Conductivity
New_Cond = zeros(num_links_N,1);
for i = 1:num_links
    New_Cond((i-1)*(n_n+1)+1:i*(n_n+1)) = Cond(i);
end

%% FEM
% Temperature
Temperature = zeros(num_nodes_N,2);
Temperature(:,1) = 1:num_nodes_N;

Temperature(8,2) = 22;
Temperature(15,2) = 22;
Temperature(num_nodes_outer+4,2) = 22;
Temperature(num_nodes_outer+num_nodes_inner,2) = 22;

% q
q = zeros(num_nodes_N,2);
q(:,1) = 1:num_nodes_N;
for i = 1:num_links_inner
    q(num_nodes+(num_links_outer+i-1)*n_n+1:num_nodes+(num_links_outer+i)*n_n,2) = Links_P(i)/n_n;
end

%% Stiffness Matrix
% Beams
K = zeros(2,2,num_links_N);
for i = 1:num_links_N
    K(1,:,i) = [New_A(i)*New_Cond(i)/New_L(i),-New_A(i)*New_Cond(i)/New_L(i)];
    K(2,:,i) = [-New_A(i)*New_Cond(i)/New_L(i),New_A(i)*New_Cond(i)/New_L(i)];
end

% Global
KK = zeros(num_nodes_N,num_nodes_N,num_links_N);
K_G = zeros(num_nodes_N,num_nodes_N);

for i = 1:num_links_N
    SN = New_links(i,2); EN = New_links(i,3);
    KK(SN,SN,i) = K(1,1,i);
    KK(SN,EN,i) = K(1,2,i);
    KK(EN,SN,i) = K(2,1,i);
    KK(EN,EN,i) = K(2,2,i);
    K_G = K_G + KK(:,:,i);
end

% Rearrange
K_ff = K_G;
K_fs = K_G;

N = 1;
for i = 1:num_nodes_N
    if Temperature(i,2) ~= 0
        K_ff(N,:) = [];
        K_ff(:,N) = [];
    else
        N = N+1;
    end
end

N = 1;
NN = 1;
for i = 1:num_nodes_N
    if Temperature(i,2) ~= 0
        K_fs(N,:) = [];
    else
        N = N+1;
    end
end
for i = 1:num_nodes_N
    if Temperature(i,2) == 0
        K_fs(:,NN) = [];
    else
        NN = NN+1;
    end
end

%% Temperatures

D_s = Temperature(:,2);

N = 1;
for i = 1:num_nodes_N
    if Temperature(i,2) == 0
        D_s(N,:) = [];
    else
        N = N+1;
    end
end

%% q

F_f = q(:,2);

N = 1;
for i = 1:num_nodes_N
    if Temperature(i,2) ~= 0
        F_f(N,:) = [];
    else
        N = N+1;
    end
end

%% Solve

D_f = K_ff\(F_f - K_fs*D_s);

Temp_Out_C = zeros(num_nodes_N,2);
Temp_Out_C(:,1) = Temperature(:,1);

N = 1;
M = 1;
for i = 1:num_nodes_N
    if Temperature(M,2) ~= 0
        Temp_Out_C(M,2) = Temperature(M,2);
        M = M+1;
    else
        Temp_Out_C(M,2) = D_f(N,1);
        M = M+1;
        N = N+1;
    end
end

%% Output
Results_C(1) = Temp_Out_C(2,2);
Results_C(2) = Temp_Out_C(3,2);
Results_C(3) = (Results_C(1)+Results_C(2)+sum(Temp_Out_C(num_nodes+n_n*1+1:num_nodes+n_n*2,2)))/(2+n_n);
Results_C(4) = Temp_Out_C(5,2);
Results_C(5) = Temp_Out_C(6,2);
Results_C(6) = (Results_C(4)+Results_C(5)+sum(Temp_Out_C(num_nodes+n_n*(1+N1+2)+1:num_nodes+n_n*(1+N1+3),2)))/(2+n_n);
Results_C(7) = (Temp_Out_C(17,2)+Temp_Out_C(21,2)+sum(Temp_Out_C(num_nodes+n_n*(num_links_outer+num_links_inner-1)+1:num_nodes+n_n*(num_links_outer+num_links_inner),2)))/(2+n_n);
Results_C(8) = (Temp_Out_C(17,2)+Temp_Out_C(18,2)+sum(Temp_Out_C(num_nodes+n_n*num_links_outer+1:num_nodes+n_n*(num_links_outer+1),2)))/(2+n_n);
Results_C(9) = (Temp_Out_C(18,2)+Temp_Out_C(19,2)+sum(Temp_Out_C(num_nodes+n_n*(num_links_outer+1)+1:num_nodes+n_n*(num_links_outer+2),2)))/(2+n_n);
Results_C(10) = Temp_Out_C(16,2);
Results_C(11) = (Temp_Out_C(1,2)+Temp_Out_C(2,2)+sum(Temp_Out_C(num_nodes+1:num_nodes+n_n,2)))/(2+n_n);
Results_C(12) = Temp_Out_C(4,2);
Results_C(13) = Temp_Out_C(num_nodes_outer+1,2);
Results_C(14) = (sum(Temp_Out_C(num_nodes+n_n*(num_links_outer+num_links_inner+2)+1:num_nodes+n_n*(num_links_outer+num_links_inner+3),2)))/(n_n);