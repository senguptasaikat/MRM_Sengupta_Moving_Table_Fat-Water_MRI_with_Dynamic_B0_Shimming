%% DOWNLOAD_DATA_ALL
%
% Downloads all data files from on-line storage.
%
% [SUCCESS_COUNT, FAIL_COUNT] = DOWNLOAD_DATA_ALL(FILENAME_INDEX)
%
%   Downloads files using the URL and destination folder
%   specified in the comma-separated-varaiable FILENAME_INDEX, which should
%   have the columns:
%
%   filename,URL,folder,MD5
%
% [SUCCESS_COUNT, FAIL_COUNT] = DOWNLOAD_DATA_ALL()
%
%   Uses the default './download_data_index.csv' for FILENAME_INDEX
%
%   Specified paths are assumed to be relative to the detected location of
%   this m-file script, i.e., the location of download_data_all.m
%
%   If the MD5 hash of the downloaded file matches, the function increments
%   SUCCESS_COUNT otherwise it increments FAIL_COUNT
%
%   Warning messages are generated when:
%       * file already exists locally
%       * file did not download
%       * file MD5 hash value does not match expected value
%
%   If the file storage service changes, the URLs in the index file
%   will be updated to reflect the new location.
%
% See also: DOWNLOAD_DATA, DOWNLOAD_DATA_VERIFY, DOWNLOAD_DATA_VERIFY_ALL, MD5, CSV_TO_STRUCT
%
function [success_count, fail_count] = download_data_all(filename_index)

    %% initialize counts
    success_count = 0;
    fail_count = 0;
    
	%% detect location of this script
    [pathstr_mfile, name_mfile, ext_mfile] = fileparts(mfilename('fullpath'));

    %% default filename_index
    default_filename_index = sprintf('%s/download_data_index.csv', pathstr_mfile);
    
    %% set default_filename_index if necessary
    if nargin==0,
        filename_index = default_filename_index;
    end
    
    %% load file index
    if exist(filename_index, 'file'),
        try
            download_data_index = csv_to_struct(filename_index);
        catch
            error(  sprintf('Unable to load download data index %s', filename_index) );
        end
    else    
        try
            download_data_index = csv_to_struct(default_filename_index);
        catch
            error(  sprintf('Unable to load download data index %s', default_filename_index) );
        end
    end 
    
    %% loop over all filenames in the download_data_index
    for k=1:length(download_data_index),
        filename = download_data_index(k).filename;
        success = download_data(filename,filename_index);
        if success,
            success_count = success_count + 1;
        else
            fail_count = fail_count + 1;
        end
    end