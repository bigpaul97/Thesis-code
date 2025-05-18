function [sorted_data] = chipseq_data_process_pombe(file,chr)
% Input 
% file: raw pombe psc3 chip-seq data in some messy form (e.g. GSM1370192_SG12054172_251601010199_S001_ChIP_1100_Jul11_1_2.txt
% chr: chromosome number, e.g. 'chr2'
% Note: file and chr are strings.

% Process: get the position (raw_table.text) and the log2 signal value (raw_table.float_2)
% Output: chromosome's position vs. chip signal (in resolution of 300 bp)

%% Read raw files 
raw_table = readtable(file);

pos = raw_table.text;
signal = raw_table.float_2;


%% Get chromosome and genomic position, as well as the corresponding signal
pos_str = string(pos);
pos_split = string;

for i = 1:size(pos_str,1)
    entry = strsplit(pos_str(i),'_');
    entry_size = size(entry,2);
    pos_split(i,1:entry_size) = entry;
end

% Extract out the position for the chr
index = strcmp(pos_split(:,1),chr); % the row index with the specific chromosome number
pos_data_str = pos_split(index,3); 
signal_data = signal(index);


pos_data = str2double(pos_data_str);
% pos_data = sortrows(pos_data);

%% Get the position with chip-seq signal, correspondingly, and sort
combined_data = [pos_data signal_data];
sorted_data = sortrows(combined_data,1);

end
