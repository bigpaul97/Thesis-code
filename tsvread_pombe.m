function varargout = tsvread_pombe( varargin )
% Input (varargin): a .tsv file with pairwise contact frequency data (e.g. GSM1379427_wt_999a-corrected-matrix_hic.tsv)
% Output(varargout): a matrx of Hi-C data.


if nargin == 0
    disp( 'importing most recent tsv file in local directory' );
    d = dir( '*.tsv' );
    [~,i] = sort( [d.datenum], 'descend' );
    fprintf( 'tsvread( ''%s'', 1:10 )\n', d(i(1)).name );
    eval( sprintf( 'tsvread( ''%s'', 1:10 )', d(i(1)).name ) );
    return;
end

fid = fopen( varargin{1}, 'r' );
if nargout == 0
    fprintf( 'loading %s, and displaying sideways\n', varargin{1} );
end
varargin(1) = [];
stuff = textscan( fid, '%s', 'delimiter', '\n');
stuff = stuff{1};
fclose(fid);

numrows = 1:size(stuff,1);
ind = cellfun( 'isclass', varargin, 'double' );
if any( ind )
    numrows = varargin{ind};
    varargin(ind) = [];
    if numel(numrows) == 1
        numrows = 1:min(numrows, size( stuff, 1) );
    end
end
numrows = intersect( numrows, 1:size( stuff, 1 ) );
header = regexprep( regexp( stuff{1}, '[^\t]*\t', 'match' ), '\t', '' );
raw = repmat( {}, numel(numrows), numel(header) );
for i=numrows
    stuff{i}(end+1) = 9;
    tmp = regexprep( regexp( stuff{i}, '[^\t]*\t', 'match' ), '\t', '');
    raw(i,1:numel(tmp)) = tmp;
end
header = raw(1,:);
data = nan(size(raw));
for i=numrows
    for j=1:size( raw, 2 )
        if ~isempty( raw{i,j} )
            [a, count, errmsg] = sscanf( raw{i,j}, '%f' );
            if ~isempty( a )
                data(i,j) = a;
            end
        end
    end
end
%%
if numel( varargin ) == 0
    j=1:size(data, 2);
else
    j = [];
    for i=1:numel(varargin)
        j = [j, find( strncmp( header, varargin{i}, numel(varargin{i})) )];
    end
end
if nargout == 0
    disp( [ strvcat( header(j)), num2str( data(numrows,j)' ) ] );
    return;
end
if numel(varargin)==0
    varargout = {data, header, raw};
    varargout = varargout(1:nargout);
    return;
end
varargout = {};
for i=j
    varargout{end+1} = data( :, i );
end
varargout{end+1} = header(:,j);
varargout{end+1} = raw(:,j);
