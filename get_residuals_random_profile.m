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

clear all; close all; fclose all; format compact; clc

%% define network parameters
% network = 'networks/generic_2.inp'; % INP file definintg the network
network = 'networks/generic_1.inp';   % INP file definintg the network
report = 'logs/generic.txt';          % EPANET log file
K=32:42;                              % interval over which to average the values of the node heads
train_no_realisations=2400;           % total number of residuals to be computed for the train set
test_no_realisations=3200;            % total number of residuals to be computed for the test set      
save_transitional_data=0;             % in addition to the steady-state data, also save the transitional data; use sparringly, it may lead to large files
save_name='generic_1';

network = 'networks/hanoi.inp';       % INP file definintg the network
report = 'logs/hanoi.txt';            % EPANET log file
K=12:18;                              % interval over which to average the values of the node heads
train_no_realisations=372;            % total number of residuals to be computed for the train set
test_no_realisations=745;             % total number of residuals to be computed for the test set      
save_transitional_data=1;             % in addition to the steady-state data, also save the transitional data; use sparringly, it may lead to large files
save_name='hanoi';


% The number of junctions in a network equals the number of nodes minus the number of tanks and reservoirs.
d=epanet(network);
junctions =  double(d.getNodeCount() - d.getNodeTankReservoirCount());
d.closeNetwork();

no_profiles=10;                       % the first is nominal (stored in the INP file), the others randomized wrt the nominal one
index=1;                              % index for the active pattern selected from the INP file

max_emitter_value=32;                 % number of emitter values
 
%% for large scale networks we cannot really take all posible combinations of profiles, faults and estimators, so we need to select a subset (randomly or not)
% each row represents a particular combination of profile, node under fault and emitter value

% first we consider pre-train data which we take to mean profile=1 (the
% nominal one), all possible faults and emitter values taken at 8:2:30
pretrain_profiles=1;
pretrain_faults=1:junctions;
pretrain_emitters=8:2:30;
realisations=...
    [kron(pretrain_profiles',ones(length(pretrain_emitters)*length(pretrain_faults),1))...
    repmat(kron(pretrain_faults',ones(length(pretrain_emitters),1)),length(pretrain_profiles),1)...
    repmat(repmat(pretrain_emitters',length(pretrain_faults),1),length(pretrain_profiles),1)];
pretrain_no_realisations=size(realisations,1);

% next we consider the train datasest, we cannot take all possible
% combinations so we make a random selection but still we assume that we
% consider the same number of cases for each node
train_profiles=1:no_profiles;
train_emitters=9:2:31;
train_faults=1:junctions;

realisations=[realisations;...
    [train_profiles(randi(length(train_profiles),train_no_realisations,1))' ...
    kron(train_faults',ones(floor(train_no_realisations/length(train_faults)),1)) ...
    train_emitters(randi(length(train_emitters),train_no_realisations,1))']];

% for the test data the same procedure repeats but this time we let
% everything to vary, including the nodes that may or may not be under fault
test_profiles=1:no_profiles;
test_emitters=8:31;
test_faults=1:junctions;
realisations=[realisations;...
    [test_profiles(randi(length(test_profiles),test_no_realisations,1))' ...
    test_faults(randi(length(test_faults),test_no_realisations,1))' ...
    test_emitters(randi(length(test_emitters),test_no_realisations,1))']];

% each time this cell is run a different combination of realisations is
% obtained; for further testing we may consider a realisations matrix a
% priori computed and saved in the realisations_data.mat file
% load ./data/realisations_data

%% get the nominal profile and create randomized perturbations around its nominal value

% read the nominal pattern
d = epanet(network);
% prepare the pattern profiles
d.openHydraulicAnalysis();
% get the values of the active pattern (the one of index "index")
pattern=[];
periods=d.getPatternLengths(index);
for p=1:periods
    pattern=[pattern d.getPatternValue(index,p)];
end

% randomize the active pattern by adding uniform noise in the [-2.5%, 2.5%] range of the nominal value
patterns=pattern;
for p=2:no_profiles
    patterns=[patterns; pattern.*(1-0.05*(1-2*rand(1,periods)))];
end
d.closeHydraulicAnalysis();
d.closeNetwork();

%% Get residues via EPANET emulation

% run for each combination of profile, fault and emitter value
P_all_trans={};
index=1;
P_all=zeros(size(realisations,1),junctions);
old_profile=-1;

%{
% for plotting reasons we need the node values obtained under nominal
% profile and nominal functioning (profile index = 1, emitter value =0); to
% obtain these values we hard code the realisations and save_name values
realisations=[1 1 0];
save_name='hanoi_nominal'
%}
for r=1:size(realisations,1)
    p=realisations(r,1);
    f=realisations(r,2);
    e=realisations(r,3);
    
    d = epanet(network);
    d.openHydraulicAnalysis();
    
    if p~=old_profile
        % set the active pattern for each of the junction nodes
        for m=1:junctions
            d.setNodeDemandPatternIndex(m, index);
        end
        d.setPattern(index,patterns(p,:))
        old_profile=p;
    end
    % node under fault
    d.setNodeEmitterCoeff(f, e);
    
    % Get element of interest in each junction;
    [P,Ptrans] = all_junctions_get_element(d,junctions,K);
    
    d.setNodeEmitterCoeff(f, 0);  % reset to healthy functioning
    
    disp(['r=' num2str(r) '   profile=' num2str(p)  '  node under fault=' num2str(f) '  emitter value=' num2str(e)])
    
    P_all(r,:)=P;
    if save_transitional_data==1
        P_trans_all{r}=Ptrans;
    end
    
    d.closeHydraulicAnalysis();
    d.closeNetwork();
end

if ~exist('./data')
    mkdir('./data')
end
% save the data into a mat file for further use
save(['./data/residues_' save_name]);
