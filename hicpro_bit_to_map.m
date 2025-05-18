function [hic_count] = hicpro_bit_to_map(data,bin_range)
% convert bit-by-bit hic data in to matrix with given resolution
% data: .matrix file ouput by HiC-Pro (e.g. 4000671_500.matrix)
% bit_range: the range of bins of for the chromosome of interest (e.g. mouse chr12: 44828-56840)

%% load .matrix file
data = load(data);

%% find the corresponding row_range for the given bit_range
% row_low_bound = find(data(:,1)>=bin_range(1),1,'first');
% row_upp_bound = find(data(:,1)<=bin_range(end),1,'last');

data_chr = data(logical(ismember(data(:,1),bin_range) .* ismember(data(:,2),bin_range)),:);

% Only select pairs within the same chromosome (intra-contact)
% data_chr = data(row_low_bound:row_upp_bound,:);
% data_chr = data_chr(ismember(data_chr(:,2), bin_range),:);

% re-index the bin pairs
data_chr(:,1) = data_chr(:,1) - bin_range(1) + 1;
data_chr(:,2) = data_chr(:,2) - bin_range(1) + 1;


%% Matrix initialization
hic_matrix = zeros(length(bin_range),length(bin_range));

%% Assign contact number to hic_matrix

for i = 1:size(data_chr,1)
    hic_matrix(data_chr(i,1),data_chr(i,2)) = data_chr(i,3);
end

hic_count = hic_matrix + hic_matrix'; % make it symmetric


%% Normalization by the largest contact
% cp_map = hic_count ./ max(reshape(hic_count,[],1));


end