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

close all; clear all; clc

%% scripts to be run beforehand

% the mat files used here are obtained through repeated runs of the
% get_residuals_random_profile.m and run_fdi_toddler.m, as follows:

% 1. run get_residuals_random_profile.m with the variables containing the
% 'hanoi' word uncommented; the result is the residues_hanoi.mat file,
% saved in the ./data subfolder

% 2. run get_residuals_random_profile.m with the variables containing the
% 'hanoi_nominal' word uncommented (lines 121-122); the result is the residues_hanoi_nominal.mat file,
% saved in the ./data subfolder

% 3. run get_residuals_random_profile.m with the variables containing the
% 'generic_2' word uncommented; the result is the residues_generic_2.mat file,
% saved in the ./data subfolder

% 4. run run_fdi_toddler.m with the variables containing the 'generic 2' word
% uncommented; the result are mat files in the format data_dl_generic_2_s.mat where 's' is the replaced by the sensor number (by default, a number from 5 to 50), saved in the
% '/data/absolute residual generic 2/' subfolder

% 5. change fault_at_node_level=0 and run the previous step; the result are mat files in the format data_dl_generic_2_s_comm.mat where 's' is the replaced by the sensor number (by default, a number from 5 to 50), saved in the
% '/data/absolute residual generic 2 comm/' subfolder; the file
% ./data/generic_2-communities.mat has to exist (obtained through the
% python communities script)


%% plot the nominal profile for the hanoi network (figure 2a)

load('./data/residues_hanoi.mat')
figure; hold on; grid on
stairs(patterns(1,:))


%% plot the nodes' heads while under the nominal profile and with nominal functioning (figure 2b)
load('./data/residues_hanoi_nominal.mat')
figure; hold on; grid on
plot(P_trans_all{1}')

%% plot for a given node the residuals for a subset of fault magnitudes (figure 3)
figure; hold on; grid on

node_index=17;
emitter_values=[4 8 12 20];

load('./data/residues_hanoi.mat')
ii=find(realisations(:,2)==node_index); % get all realisations corresponding to node_index
% plot only these emitter values that we are interested in 
for e=emitter_values
    jj=find(realisations(ii,3)==e);
    if ~isempty(jj)
        for ind=1:length(ii(jj))            
            plot(P_trans_all{ii(jj(ind))}')
        end
    end
end

% add also the nominal case
load('./data/residues_hanoi_nominal.mat')
ii=find(realisations(:,2)==node_index); % get all realisations corresponding to node_index
% plot only these emitter values that we are interested in 
e=0;
jj=find(realisations(ii,3)==e);
if ~isempty(jj)
    for ind=1:length(ii(jj))
        plot(P_trans_all{ii(jj(ind))}')
    end
end

%% plot the residuals for a certain fault and emitter selection obtained from the steady-state information (figure 4)

figure; hold on; grid on

% first get the value of the nominal head values
load('./data/residues_hanoi_nominal.mat')
P_nominal=P_all(1,:);

load('./data/residues_hanoi.mat')
ii=find(realisations(:,2)==node_index); % get all realisations corresponding to node_index
% plot only these emitter values that we are interested in 
for e=emitter_values
    jj=find(realisations(ii,3)==e);
    if ~isempty(jj)
        for ind=1:length(ii(jj))            
            plot((P_all(ii(jj(ind)),:)-P_nominal)')
        end
    end
end

%% plot fault signature matrix (figure 5)

load('./data/residues_hanoi.mat')
% construct the absolute residuals, based on the nominal profile data
R=P_all(1:pretrain_no_realisations,:)-repmat(P_nominal,pretrain_no_realisations,1);
R=R';

% take the threshold tau and construct the binarized signature matrix
tau=3;
M=zeros(size(R));
M=abs(R)>=tau;

% M has a column for each of the fault magnitude, we wish to simplify and
% have a column for each fault (hence all magnitudes are grouped together)

% since the data is saved for the nominal profile in continuous
% sub-matrices (all emitter values are enumerated for a given fault we can
% easily find these groupings)

group=pretrain_no_realisations/junctions;
Mg=zeros(size(M,1),size(M,2)/group);

for j = 1:size(M,2)/group % numberof_rows    
    Mg(:,j)=double(sum(M(:,(j-1)*group+1:j*group)')'>=1);
end

figure
imshow(1-Mg)

%% sensor selections with various methods (figure 6)

ss=2:10;
for si=1:length(ss)
    s=ss(si);
    s_mat_proj_lambda = 10000; % node distance penalty
    
    [D_laplacian,A_laplacian,G] = build_distances(network);
    [selected_sensors_GGS, err]=s_alloc_mat_proj(R,s,s_mat_proj_lambda,D_laplacian);
    selected_sensors_GGS_all{si}=sort(selected_sensors_GGS);
    
    tau=3;
    M=zeros(size(R));
    M=abs(R)>=tau;
    group=pretrain_no_realisations/junctions;
    Mg=zeros(size(M,1),size(M,2)/group);
    
    for j = 1:size(M,2)/group % number of_rows
        Mg(:,j)=double(sum(M(:,(j-1)*group+1:j*group)')'>=1);
    end
    [selected_sensors_MSC,failures_MSC]=sensor_selection_MI(Mg,s);
    selected_sensors_MSC_all{si}=selected_sensors_MSC;
    
    combinations=nchoosek(1:size(Mg,2),2);
    Mgtil=zeros(size(Mg,1),size(combinations,1));
    for i=1:size(combinations,1)
        Mgtil(:,i)=Mg(:,combinations(i,1))+Mg(:,combinations(i,2))- 2*Mg(:,combinations(i,1)).*Mg(:,combinations(i,2)); % XOR
    end
    
    [selected_sensors_MTC,failures_MTC]=sensor_selection_MI(Mgtil,s);
    selected_sensors_MTC_all{si}=selected_sensors_MTC;
end

figure; hold on; grid on
for si=1:length(ss)
    scatter(selected_sensors_MSC_all{si},ss(si)*ones(1,length(selected_sensors_MSC_all{si})),'bo')
    scatter(selected_sensors_MTC_all{si},ss(si)*ones(1,length(selected_sensors_MTC_all{si})),'bx')
    scatter(selected_sensors_GGS_all{si},ss(si)*ones(1,length(selected_sensors_GGS_all{si})),'bs')
end

%% illustration of fault detection and isolation for the hanoi network (figure 7)
figure; hold on; grid on

node_index=22;
plot(R(:,find(realisations(:,2)==node_index,1))','b*-')
% and show the selectd nodes when we take 5 (corresponds to the 4th element
% inthe cell array)
scatter(selected_sensors_GGS_all{4},R(selected_sensors_GGS_all{4},find(realisations(:,2)==node_index,1))','o','filled','r')
a=axis;
axis([1 junctions a(3) a(4)])

%% plot success rates for the generic_2 network (figure 9)

% pass through each of the mat files contained in the given folder and
% extract the success rates
files=dir('./data/absolute residual generic 2');
mat_files=numel(find([files.isdir]==0));

success0_all=zeros(1,mat_files);
success1_all=zeros(1,mat_files);
success2_all=zeros(1,mat_files);

for i=1:length(files)
    file=files(i);
    if file.isdir
        continue
    end
    load([file.folder '/' file.name])
    success0_all(length(selected_sensors))=success0;
    if success1>100
        success1_all(length(selected_sensors))=success1*100/length(estimate); %<== because I initially forgot to norm it; now is corrected; this should not appear in newly obtained mat files
    else
        success1_all(length(selected_sensors))=success1;
    end
    if success2>100
        success2_all(length(selected_sensors))=success2*100/length(estimate); %<== because I initially forgot to norm it; now is corrected; this should not appear in newly obtained mat files
    else
        success2_all(length(selected_sensors))=success2;
    end    
end

%
hold on; grid on
plot(5:50,success0_all(5:end),'b*-')
plot(5:50,success1_all(5:end),'r*-')
plot(5:50,success2_all(5:end),'g*-')

files=dir('./data/absolute residual generic 2 comm');
mat_files=numel(find([files.isdir]==0));
success0_all=zeros(1,mat_files);

for i=1:length(files)
    file=files(i);
    if file.isdir
        continue
    end
    load([file.folder '/' file.name])
    success0_all(length(selected_sensors))=success0;   
end
plot(5:50,success0_all(5:end),'m*-')

%% plot sparsity and details for the dictionary learning for the generic_2 network (figure 10)

%load a particular run (eg, for s=10) and anayze the classification results
%(by checking the sparse classification matrix X)

load('./data/residues_generic_2.mat')
load('./data/absolute residual generic 2/data_dl_generic_2_s=10.mat')

figure; 

subplot(4,1,[1 2]); grid on; hold on
% plot the active atoms for each realisation
[~,ind]=sort(labels_test);
[ii,jj]=find(X(:,ind)>1e-3);
scatter(jj,ii,'filled','b')
axis([1 size(X,2) 1 size(X,1)])

subplot(4,1,3); grid on; hold on
% plot active atoms for class 140
ll=labels_test(ind);
class=140;%129;
ind2=find(ll==class);
XX=X(:,ind);
[ii,jj]=find(XX(:,ind2)>1e-3);
scatter(ii,ind2(1)+jj-1,'filled','b')
a=axis;
axis([1 size(D,2) a(3:4)])

% plot classification values for class 140
subplot(4,1,4); grid on; hold on
tmp=W*XX(:,ind2(1)+jj-1);
tmp2=1:junctions;
for i=1:size(tmp,2)
    ii=find(abs(tmp(:,i))>=2e-4);    
    stem(tmp2(ii),tmp(ii,i),'b','filled')
end
