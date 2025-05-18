function [binned_data] = bin_chipseq_narrowpeak(filename,chr,target_bin,start_p,end_p)
% Bin the rawdata to target bin by summing the signal in that whole bin

% filename: .narrowpeak file (e.g. GSM6052391_scc2_noIAA_1_peaks.narrowPeak)
% chr: chromosome number in roman numeral (e.g. 'chrIII')
% target_bin: bin size of the final ChIP-seq array
% start_p: start position of the target interval
% end_p: end position of the target interval

%% narrowpeak file processing

% Open the file
fid = fopen(filename, 'r');
% Define the format of the .narrowPeak file
formatSpec = '%s%f%f%s%f%s%f%f%f%f';
% Read the file using textscan
data = textscan(fid, formatSpec, 'Delimiter', '\t');
% Close the file
fclose(fid);
% Convert the cell array to a table for easier manipulation
dataTable = table(data{1}, data{2}, data{3}, data{4}, data{5}, data{6}, data{7}, data{8}, data{9}, data{10}, ...
'VariableNames', {'chrom', 'chromStart', 'chromEnd', 'name', 'score', 'strand', 'signalValue', 'pValue', 'qValue', 'peak'});

rawdata_chr = dataTable(strcmp(dataTable.chrom, chr),:);

% Make sure the data is sorted from lower end to higher end
sorted_data = table2array(rawdata_chr(:,[2,3,7]));
bp_data = zeros(end_p-start_p,1);


%% Main loop to bin the data

for i = 1:size(sorted_data,1)
    
    l_idx = sorted_data(i,1);
    r_idx = sorted_data(i,2);
    bp_data(l_idx:r_idx) = sorted_data(i,3);

end

%% Bin data to target_bin size

binned_data = zeros((end_p-start_p)/target_bin,2);
binned_data(:,1) = 1:size(binned_data,1);

for i = 1:size(binned_data,1)
    binned_data(i,2) = sum(bp_data(1+target_bin*(i-1):target_bin*i));
end

% Normalization
binned_data(:,2) = binned_data(:,2)/max(binned_data(:,2)).*100;

end