function [data] = add_gap_point(data,binsize)
% Put chip-seq signal = 0 for those missing points
% Also, convert log2 signal to relative coverage
% data: the two-column signal (e.g. chr2_mis4_array.mat)
% binsize: binsize of the first column of original <data>



base_number = data(:,1)/binsize;

data(:,1) = base_number;


%% remake the filled data with 0's at the missing spots


base_diff = data(2:end,1)-data(1:end-1,1);
pos = find(base_diff ~= 1,1,'first');

while ~isempty(pos)
    gap = base_diff(pos);
    
    inserted_rows = [pos+1:pos+gap-1;zeros(gap-1,1)']';
    data = [data(1:pos,:); inserted_rows; data(pos+1:end,:)];
    
    base_diff = data(2:end,1)-data(1:end-1,1);
    pos = find(base_diff ~= 1,1,'first');

end

data(:,1) = data(:,1).*binsize;
data = [[0,0]; data];


% data(data(:,2)==0,2) = -Inf;
% data(:,2) = 2.^data(:,2);

% filled_data = [];
% filled_data(1,:) = data(1,:);
% 
% for i = 1:length(base_diff)
%     base_diff(i)
% end



end