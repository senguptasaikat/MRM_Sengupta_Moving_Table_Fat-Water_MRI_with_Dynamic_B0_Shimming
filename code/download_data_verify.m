%% DOWNLOAD_DATA_VERIFY
%
% Verifies download data file.
%
% [VERIFIED] = DOWNLOAD_DATA_VERIFY(FILENAME_VERIFY,FILENAME_INDEX)
%
%   Verifies files using the MD5 checksum and destination folder
%   specified in the comma-separated-varaiable FILENAME_INDEX, which should
%   have the columns:
%
%   filename,URL,folder,MD5
%
% [VERIFIED] = DOWNLOAD_DATA_VERIFY(FILENAME)
%
%   Uses the default './download_data_index.csv' for FILENAME_INDEX
%
%   Specified paths are assumed to be relative to the detected location of
%   this m-file script, i.e., the location of download_data_verify.m
%
%   Warning messages are generated when:
%       * file requested is not in the index
%       * file does not exist
%       * file MD5 hash value does not match expected value
%
%   If the MD5 hash of the downloaded file matches, the function returns
%   VERIFIED=true otherwise VERIFIED=false
%
% See also: DOWNLOAD_DATA, DOWNLOAD_DATA_ALL, DOWNLOAD_DATA_VERIFY_ALL, MD5, CSV_TO_STRUCT
%
function [verified] = download_data_verify(filename_verify,filename_index)

    %% initialize verified to false
    verified = false;
    
    %% detect location of this script
    [pathstr_mfile, name_mfile, ext_mfile] = fileparts(mfilename('fullpath'));

    %% default filename_index
    default_filename_index = sprintf('%s/download_data_index.csv', pathstr_mfile);
    
    %% set default_filename_index if necessary
    if nargin==1,
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
        
    %% verify if the file is a valid option
    [valid, idx] = ismember(filename_verify, {download_data_index(:).filename} );
    
    if valid,

        %% file details
        folder_relative_path = download_data_index(idx).folder;
        md5_expected =  download_data_index(idx).MD5;
        
        destination_filename_fullpath = sprintf('%s/%s/%s', pathstr_mfile, folder_relative_path, filename_verify);

        %% warn if the file still does not exist
        if exist(destination_filename_fullpath, 'file'),
            
            %% success only if MD5 checksum matches
            md5_download = MD5(destination_filename_fullpath);
            if all(md5_expected==md5_download),
                verified = true;
                disp( sprintf('Success! MD5 checksum verified for %s : %s (expected) == %s (downloaded)', filename_verify, md5_expected, md5_download) );
            else
                warning( sprintf('Invalid MD5 checksum for %s: %s (expected) != %s (downloaded)', filename_verify, md5_expected, md5_download) );
            end
            
        else
            warning( sprintf('%s does not exist', filename_verify) );
        end

    else
        warning( sprintf('%s is not a valid file to verify', filename_verify) );
    end