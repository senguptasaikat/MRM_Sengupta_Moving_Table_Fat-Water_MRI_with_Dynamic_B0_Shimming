#!/usr/bin/env python
"""
Downloads all data files from on-line storage using the URL and destination
folder specified in the comma-separated-varaiable './download_data_index.csv'
file, which should have the columns:

    filename,URL,folder,MD5

Specified paths are assumed to be relative to the detected location of
this python script, i.e., the location of download_data_all.py

If the MD5 hash of the downloaded file matches, the function increments the
reported SUCCESS_COUNT otherwise it increments the reported FAIL_COUNT

Warning messages are generated when:
       * file already exists locally
       * file did not download
       * file MD5 hash value does not match expected value

   If the file storage service changes, the URLs in the index file
   will be updated to reflect the new location.

 See also: download_data_verify_all.py
"""
import pandas as pd
import os
import urllib
import hashlib

# detect location of this script
pathstr_python_script = os.path.dirname(os.path.realpath(__file__))

# full path to download_data_index.csv
filename_index = os.path.join(pathstr_python_script,'download_data_index.csv')

# load download_data_index.csv to a pandas dataframe
file_index = pd.read_csv(filename_index)

# initialize counts
success_count = 0;
fail_count = 0;

# loop over files
for idx, filename_download in enumerate(file_index['filename']):
    
    print("")    
    
    # download details
    URL = file_index['URL'][idx]
    folder_relative_path = file_index['folder'][idx]
    md5_expected = file_index['MD5'][idx]
    
    destination_filename_fullpath = os.path.join(pathstr_python_script, folder_relative_path, filename_download)
    
    # warn if the file already exists
    if os.path.isfile(destination_filename_fullpath):
        print("%s already exists" % (destination_filename_fullpath) )
    
    # download file
    print("Downloading file %s from %s" % (filename_download, URL) )
    filehandle = urllib.URLopener()
    filehandle.retrieve(URL, destination_filename_fullpath)
    
    if os.path.isfile(destination_filename_fullpath):
        print("%s already exists" % (destination_filename_fullpath) )
        
        # success only if MD5 checksum matches
        md5_download = hashlib.md5(open(destination_filename_fullpath, 'rb').read()).hexdigest()
        if md5_expected==md5_download:
            print("Success! MD5 checksum verified for %35s : %s (expected) == %s (downloaded)" % (filename_download, md5_expected, md5_download) )
            success_count += 1
        else:
            print("INVALID MD5 checksum for            %35s: %s (expected) != %s (downloaded)" % (filename_download, md5_expected, md5_download) )   
            fail_count += 1
    else:
        print("%s did not download" % (filename_download) )
        fail_count += 1
        
print("\nSuccessfully downloaded %.2f%% : success_count = %d, fail_count = %d\n" % (100.0 * success_count/float(success_count+fail_count), success_count, fail_count) )