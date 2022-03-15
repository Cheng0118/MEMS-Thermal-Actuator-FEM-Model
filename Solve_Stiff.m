%% Mechanical
% Loads
Mechanical_Load = zeros(num_nodes,6);
Mechanical_Load(num_nodes_outer,2) = -0.02;

% Mechanical
N_2 = 1;
ff_m = zeros(6*num_nodes,1);
for i = 1:num_nodes
    for ii = 1:6
        if DOFs(i,ii) == 0
           ff_m(N_2,:) = [];
        else
           ff_m(N_2,1) = Mechanical_Load(i,ii);
           N_2 = N_2+1;
        end
    end
end

%% Mechanical
Displacement = K_ff\ff_m;
Disp_Out_stiff = Displacement;

Outputs = zeros(6*num_nodes,1);
Outputs(6*(num_nodes_outer-1)+2,1) = 1;

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
        Disp_Out_stiff(N_5,:) = [];
    else
        N_5 = N_5+1;
    end
end

%% Results

Mech_stiffness = (-0.04/Disp_Out_stiff);

F_Output = Stroke*Mech_stiffness;