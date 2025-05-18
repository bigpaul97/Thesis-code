function [output] = shift_rescale(data,min_value,shift_unit)
% data is in the form of two-column vector
% first column = integer lattice index, second column = chip-seq
% shift may be downward or upward: down = shift_unit(-); up = (+)
% after the shift, the max is rescaled to 100.
% If there are values below min_value, make it a min_value


%% Initialization

output = zeros(size(data,1),2);
data(:,2) = data(:,2)/max(data(:,2))*100; % make max chipseq 100

%% rescale
shifted = data(:,2)+shift_unit;
shifted = shifted/(100+shift_unit)*100;
shifted(shifted < min_value) = min_value;

output(:,2) = shifted;
output(:,1) = 1:size(data,1);


end