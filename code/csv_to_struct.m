%% CSV_TO_STRUCT
%
% Converts download data index CSV file to a MATLAB structure.
%
% [CSV_STRUCT] = CSV_TO_STRUCT(FILENAME)
%
% See also: DOWNLOAD_DATA, DOWNLOAD_DATA_ALL, DOWNLOAD_DATA_VERIFY, DOWNLOAD_DATA_VERIFY_ALL, MD5
%
function [csv_struct] = csv_to_struct(filename)
    
    fid = fopen(filename);
    all_entries = textscan(fid,'%s','delimiter',',');
    fclose(fid);

    fid = fopen(filename);
    non_header_entries = textscan(fid,'%s','delimiter',',','HeaderLines',1);
    fclose(fid);

    num_cols = numel(all_entries{1}) - numel(non_header_entries{1});
    col_names = all_entries{1}(1:num_cols);
    num_rows = numel(non_header_entries{1})/num_cols;
    csv_struct = cell2struct(reshape(non_header_entries{1},num_cols,num_rows),col_names,1);
    