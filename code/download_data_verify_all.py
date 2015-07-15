#!/usr/bin/env python
"""
Verfies all data files from on-line storage using the MD5 and destination
folder specified in the comma-separated-varaiable './download_data_index.csv'
file, which should have the columns:

    filename,URL,folder,MD5

Specified paths are assumed to be relative to the detected location of
this python script, i.e., the location of download_verify_all.py

If the MD5 hash of the file matches, the function increments the
reported VERIFIED_COUNT otherwise it increments the reported FAIL_COUNT

Warning messages are generated when:
       * file does not exist
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
verified_count = 0;
fail_count = 0;

# loop over files
for idx, filename_verify in enumerate(file_index['filename']):
    
    #print("\n")    
    
    # download details
    folder_relative_path = file_index['folder'][idx]
    md5_expected = file_index['MD5'][idx]
    
    destination_filename_fullpath = os.path.join(pathstr_python_script, folder_relative_path, filename_verify)
      
    if os.path.isfile(destination_filename_fullpath):
  
        # success only if MD5 checksum matches
        md5_download = hashlib.md5(open(destination_filename_fullpath, 'rb').read()).hexdigest()
        if md5_expected==md5_download:
            print("Success! MD5 checksum verified for %35s : %s (expected) == %s (downloaded)" % (filename_verify, md5_expected, md5_download) )
            verified_count += 1
        else:
            print("INVALID MD5 checksum for           %35s: %s (expected) != %s (downloaded)" % (filename_verify, md5_expected, md5_download) )   
            fail_count += 1
    else:
        print("%s does not exist" % (filename_verify) )
        fail_count += 1
        
print("\nSuccessfully verfied %.2f%% : verified_count = %d, fail_count = %d\n" % (100.0 * verified_count/float(verified_count+fail_count), verified_count, fail_count) )