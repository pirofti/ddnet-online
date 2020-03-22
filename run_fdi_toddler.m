% Copyright (c) 2020 Paul Irofti <paul@irofti.net>
% Copyright (c) 2020 Florin Stoican <florin.stoican@acse.pub.ro>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

%% Test FDI results with TODDLeR algorithm
clear all; close all; fclose all; format compact;

%% ATTENTION:
% You will need a copy of the dictionary learning toolbox before being able
% to run the scripts. Execute setupDL.m from there before running this.

%% Load dataset
% water network given as an INP file, necessary when using the Graph-GS method for sensor placement

% network = 'networks/generic_2.inp'; 
network = 'networks/hanoi.inp'; 

% load residuals and auxiliary variables, this is the data used for
% pretrain, train and test for the DL algorithms
% relevant variables are P_all (matrix gathering all the residuals); 
% the numbers of pretrain, train and test realisations (pretrain_no_realisations, train_no_realisations, test_no_realisations)
% the realisations variable which characterizes each of the residuals by
% node under fault, active profile and emitter value

% load('./data/residues_generic_2.mat');  % mat file obtained from running 'get_residuals_random_profile.m'
load('./data/residues_hanoi.mat');  % mat file obtained from running 'get_residuals_random_profile.m'

% load communities within which the network is divided
if exist('./data/generic_2-communities.mat')
    load('./data/generic_2-communities.mat');
end

%% parameters
% % save_name='generic_2';  % root name for the mat file in which the obtained data will be saved
% save_name='hanoi';  % root name for the mat file in which the obtained data will be saved
fault_at_node_level=1;  % 1 (default) means that each fault is associated to a node; any other value means that the fault is associated to a community

if fault_at_node_level~=1
    save_name=[save_name '_comm_'];
end
%% associate each node to a fault class
nodes=[];
nodes_index=[];
if fault_at_node_level==1
    % this is the 'each node is a fault class' case
    nodes=1:junctions;
    nodes_index=1:junctions;
else
    % this the community case
    for i=1:length(comm_nodes)
        nodes=[nodes comm_nodes{i}+1];
        nodes_index=[nodes_index i*ones(1,length(comm_nodes{i}))];
    end
end
%% construct matrices and their labels for pre-train, train and test
interval=1:pretrain_no_realisations;
R_pre_train=P_all(interval,:)';%-repmat(P_nominal',1,length(interval));
labels_pre=[];
for i=interval
    labels_pre=[labels_pre nodes_index(find(nodes==realisations(i,2)))];
end
% labels_pre=realisations(interval,2); % simpler assignation which works in the case node=fault case

interval=pretrain_no_realisations+1:pretrain_no_realisations+train_no_realisations;
R_post_train=P_all(interval,:)';%-repmat(P_nominal',1,length(interval));
labels_post=[];
for i=interval
    labels_post=[labels_post nodes_index(find(nodes==realisations(i,2)))];
end
% labels_post=realisations(interval,2); % simpler assignation which works in the case node=fault case

interval=pretrain_no_realisations+train_no_realisations+1:pretrain_no_realisations+train_no_realisations+test_no_realisations;
R_test=P_all(interval,:)';%-repmat(P_nominal',1,length(interval));
labels_test=[];
for i=interval
    labels_test=[labels_test nodes_index(find(nodes==realisations(i,2)))];
end
% labels_test=realisations(interval,2); % simpler assignation which works in the case node=fault case

%% construct dictionary, run classification for varying number of sensors
ss=5:50;

for sensor=ss
    
    % apply the graph-GS methods for sensor placement; MSC and MCT methods
    % may be used as well but they are much slower or eve
    s_mat_proj_lambda = 10000; % node distance penalty    
    [D_laplacian,A_laplacian,G] = build_distances(network);
    [selected_sensors_GGS, err]=s_alloc_mat_proj(R_pre_train,sensor,s_mat_proj_lambda,D_laplacian);
    selected_sensors=sort(selected_sensors_GGS);
    
    %%
    Y_pre=normc(double(R_pre_train(selected_sensors,:)));
    Y=normc(double([R_post_train(selected_sensors,:) R_test(selected_sensors,:)]));
    
    %% Pretraining
    % Run Label Consistent DL Classification (with shared dictionary) on
    % labeled signals
    [D, W, A] = toddler_pretrain(Y_pre, labels_pre);
    
    %% Run TODDLeR
    index=length(labels_test); % index = 3287 for the figures used in the journal paper
    [estimate, D, W, A, X] = toddler(Y(:,1:index+length(labels_post)), D, W, A, 'Classes', numel(unique(nodes_index)), 'Method','G2', ...
        'Constraint',{4,16,8,8}, 'Forget', 0.9, ...
        'Labels', labels_post);
    estimate(1:length(labels_post))=[];
    X(:,1:length(labels_post))=[];

    % Succesful classification percentage
    
    success0=0;     % classification is successful only if the estimation is equal with the node under fault or is inside the cummunity containing the node
    success1=0;     % classification is successful if the estimation is equal to the node or any of its neighbors
    success2=0;     % classification is successful if the estimation is equal to the node or any of its neighbors or any of the neighbirs' neighbors
    for i=1:length(estimate)
        if labels_test(i)==estimate(i)
            success0=success0+1;
        end
        if any(unique([labels_test(i); neighbors(G,labels_test(i))])==estimate(i))
            success1=success1+1;
        end
        if any(unique(cell2mat(arrayfun(@(x) neighbors(G,x), [labels_test(i); neighbors(G, labels_test(i))],'UniformOutput',false)))==estimate(i))
            success2=success2+1;
        end
    end
    % transform into percentages
    success0=success0*100/length(estimate)
    success1=success1*100/length(estimate)
    success2=success2*100/length(estimate)
    
    if ~exist(['./data/absolute residual ' save_name])
        mkdir(['./data/absolute residual ' save_name])
    end
    
    % save the data for further analysis/plotting
    save(['./data/absolute residual ' save_name '/data_dl_' save_name '_s=' num2str(sensor)],...
        'R_pre_train', 'R_post_train', 'R_test', ...
        'labels_pre', 'labels_post', 'labels_test',...
        'selected_sensors', 'success0', 'success1', 'success2',...
        'estimate', 'D', 'W', 'A', 'X')
end
