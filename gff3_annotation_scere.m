function [gene_info_sorted,conv_index,conv_interval,conv_signal_norm] = ...
    gff3_annotation_scere(filename,genome_size,startPt,endPt,rsln)
% Calculate the convergent gene variable of targeted interval

% filename: .gff3 file containing annotated gene sequence features of an entire chromosome (e.g. All Annotated Sequence Features-chrII-1..813184.gff3)
% genome_size: # of basepair of the processed chromosome
% startPt: start location of the targeted interval
% endPt: end location of the targeted interval
% rsln: bin size of the convergent gene variable


%% Pick out the gene position and direction for analysis
gene_info = readtable(filename,'filetype','text','delimiter','\t');

gene_row = strcmp(gene_info.Var3, 'gene');

gene_info_sorted = sortrows(gene_info(gene_row,:),4); 

gene_start_sorted = str2double(string(gene_info_sorted.Var4));
gene_end_sorted = str2double(string(gene_info_sorted.Var5));
gene_dir_sorted = cell(gene_info_sorted.Var7);


%% Initialization

conv_interval = [];
conv_index = []; % give the index of the gene of the + strand gene of each convergent region

%% Main loop

for i = 1:length(gene_start_sorted)-1

    if strcmp(gene_dir_sorted(i),'-')
        continue
    end

    if strcmp(gene_dir_sorted(i+1),'+')
        continue
    end
    
    if gene_end_sorted(i) <= gene_start_sorted(i+1)
        conv_index = [conv_index; i];
        conv_interval = [conv_interval; gene_end_sorted(i), gene_start_sorted(i+1)];
    end
    
end

%% Give 1's for convergent regions

conv_signal = zeros(genome_size,1);
for i = 1:length(conv_index)
    conv_signal(conv_interval(i,1):conv_interval(i,2)) = 1;
end

%% Bin the convergent signal to desired resolution (rsln), limited by startPt and endPt

conv_signal_binned = zeros((endPt-startPt)/rsln,1);

for i = 1:length(conv_signal_binned)
    conv_signal_binned(i) = sum(conv_signal(rsln*(i-1)+1+startPt:rsln*i+startPt));
end


%% Make signal binary (comment out if don't want)
for i = 1:length(conv_signal_binned)
    if conv_signal_binned(i) > 0
        conv_signal_binned(i) = 1;
    end
end

%% Make signal normalized by convergent event count

conv_interval = ceil(conv_interval / rsln); % give the lattice index of each conv start and end

conv_signal_norm = zeros(ceil(genome_size/rsln),1); % initialize

% loop to assign weight to lattice
for i = 1:size(conv_interval,1)

    if conv_interval(i,1) == conv_interval(i,2)
        conv_signal_norm(conv_interval(i,1)) = conv_signal_norm(conv_interval(i,1)) + 1;
    else
        dilution = conv_interval(i,2) - conv_interval(i,1) + 1;
        for j = conv_interval(i,1):conv_interval(i,2)
            conv_signal_norm(j) = 1/dilution;
        end
    end

end

% pick the target region (startPt to endPt)
start_idx = startPt/rsln + 1;
end_idx = endPt/rsln;
conv_signal_norm = conv_signal_norm(start_idx:end_idx);


end


