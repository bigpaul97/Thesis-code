% clear all;
% close all;
% clc;
filename = 'GSM3656899_HiC_SSY14_ndt80D_2A2_8h_2000_iced.matrix'; % name of the .matrix file of interest
A = load(filename);
A = A'; % delete this if use iced data. (this is just format stupidity...)
B = A(1,:);
C = A(2,:);
At = A.';

%% Scan
close
nflat = size(At(:,1));
% ndim = max(B);
cmap = zeros(max(max(B),max(C)),max(max(B),max(C)));
for i=1:nflat
    cmap(At(i,1),At(i,2)) = At(i,3);
    cmap(At(i,2),At(i,1)) = At(i,3);
end

%% Plot

% hmap_log = imagesc_log(cmap,0.2);