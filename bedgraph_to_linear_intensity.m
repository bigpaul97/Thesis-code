function [output,full_list] = bedgraph_to_linear_intensity(filename,start_p,end_p,rsln)
% filename: .bedgraph file output by "bigWigToBedGraph" utility (e.g. GSM1277163_Rec8.wt.4h.saccer2.bedgraph)
% start_p: start basepair no.
% end_p: end basepair no.
% rsln: the bin size (needs to divide end_p-start_p)
% output: a linear array of intensity in the resolution specified by rsln


% bin boundary interpretation: e.g. 2-5 represents data point of 3, 4, and 5 (right binned)


%% convert bedgraph to array (change chromosome name here)
file_path = filename;
fileID = fopen(file_path, 'r');
formatSpec = '%s %d %d %f';
rawdata = textscan(fileID, formatSpec, 'Delimiter', '\t');
% chromosome = rawdata{1};
startPos = rawdata{2};
endPos = rawdata{3};
value = rawdata{4};


% sum(strcmp(chromosome, 'chrII'))
data = double([startPos,endPos]);
data(:,3) = value; % this is important, so that the value is not rounded!!!

% get rid of data outside of RoI
% data(data(:,2) < start_p,:) = [];
% data(data(:,1) > end_p,:) = [];

data(1:find(data(:,2) < start_p,1,'last'),:) = [];
data(find(data(:,1) > end_p,1,'first'):end,:) = [];

%% Initialization
lattice_num = (end_p-start_p)/rsln;
output = zeros(lattice_num,2);
output(:,1) = 1:lattice_num;

full_list = zeros(end_p - start_p,1);

%% Main loop

disp(size(data,1)) % display total labor size

for i = 1:size(data,1)
    
    % output progress
    if mod(i,10000) == 0
        disp(i)
    end

    idx_range = (data(i,1)+1:data(i,2))-start_p;
    
    if (sum(idx_range <= 0) == 0) && (sum(idx_range > end_p-start_p) == 0)

        full_list(idx_range) = data(i,3);

    elseif sum(idx_range <= 0) ~= 0

        idx_range = idx_range(idx_range>0);
        full_list(idx_range) = data(i,3);

    elseif sum(idx_range > end_p-start_p) ~= 0

        idx_range = idx_range(idx_range <= end_p-start_p);
        full_list(idx_range) = data(i,3);
    end

end

for i = 1:lattice_num
    output(i,2) = sum(full_list((i-1)*rsln+1:i*rsln));
end

% output(:,2) = output(:,2)./max(output(:,2)).*100; % scale it to max=100

end

% for i = start_p+1:end_p
%     
%     % output progress
%     if mod(i,10000) == 0
%         disp(i)
%     end
% 
%     upp_bound = find(data(:,2) >= i,1,'first');
%     low_bound = find(data(:,1) <= i,1,'last');
% %     disp(upp_bound)
% %     disp(low_bound)
%     if isempty(low_bound) || isempty(upp_bound)
%         continue
%     elseif (upp_bound ~= low_bound) && (upp_bound + 1 ~= low_bound)
% %         disp(upp_bound)
% %         disp(low_bound)
%         continue
%     end
%     
% %     disp(i)
%     indx = i-start_p;
% 
%     full_list(indx) = data(upp_bound,3);
% end


% upp_bound = find(data(:,2) <= start_p + i * rsln);
%     if data(upp_bound,2) ~= start_p + rsln * i
%         if data(upp_bound,2) >= start_p + rsln * i
%             error('dd')
%         end
%         tail_v = (start_p + rsln * i - data(upp_bound+1,1)) * data(upp_bound+1,3);
%     end