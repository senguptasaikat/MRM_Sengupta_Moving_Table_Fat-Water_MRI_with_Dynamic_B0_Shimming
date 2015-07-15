%% DOWNLOAD_DATA
%
% Downloads data files from on-line storage.
%
% [SUCCESS] = DOWNLOAD_DATA(FILENAME_DOWNLOAD,FILENAME_INDEX)
%
%   Downloads FILENAME_DOWNLOAD using the URL and destination folder
%   specified in the comma-separated-varaiable FILENAME_INDEX, which should
%   have the columns:
%
%   filename,URL,folder,MD5
%
% [SUCCESS] = DOWNLOAD_DATA(FILENAME)
%
%   Uses the default './download_data_index.csv' for FILENAME_INDEX
%
%   Specified paths are assumed to be relative to the detected location of
%   this m-file script, i.e., the location of download_data.m
%
%   If the MD5 hash of the downloaded file matches, the function returns
%   SUCCESS=true otherwise SUCCESS=false
%
%   Warning messages are generated when:
%       * file requested is not in the index
%       * file already exists locally
%       * file did not download
%       * file MD5 hash value does not match expected value
%
%   If the file storage service changes, the URLs in the index file
%   will be updated to reflect the new location.
%
% See also: DOWNLOAD_DATA_ALL, DOWNLOAD_DATA_VERIFY, DOWNLOAD_DATA_VERIFY_ALL, MD5, CSV_TO_STRUCT
%
function [success] = download_data(filename_download,filename_index)

    %% initialize success to false
    success = false;

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
        
    %% download if the file is a valid option
    [valid, idx] = ismember(filename_download, {download_data_index(:).filename} );
    
    if valid,

        %% download details
        URL = download_data_index(idx).URL;
        folder_relative_path = download_data_index(idx).folder;
        md5_expected =  download_data_index(idx).MD5;
        
        destination_filename_fullpath = sprintf('%s/%s/%s', pathstr_mfile, folder_relative_path, filename_download);
        
        %% warn if the file already exists
        if exist(destination_filename_fullpath, 'file'),
            warning(  sprintf('%s already exists', destination_filename_fullpath) );            
        end

        %% download file
        disp( sprintf('Downloading file %s from %s', filename_download, URL) );
        urlwrite(URL, destination_filename_fullpath);

        %% warn if the file still does not exist
        if exist(destination_filename_fullpath, 'file'),
            
            %% success only if MD5 checksum matches
            md5_download = MD5(destination_filename_fullpath);
            if all(md5_expected==md5_download),
                success = true;
                disp( sprintf('Success! MD5 checksum verified for %s : %s (expected) == %s (downloaded)', filename_download, md5_expected, md5_download) );
            else
                warning( sprintf('Invalid MD5 checksum for %s: %s (expected) != %s (downloaded)', filename_download, md5_expected, md5_download) );
            end
            
        else
            warning( sprintf('%s did not download', filename_download) );
        end

    else
        warning( sprintf('%s is not a valid file to download', filename_download) );
    end