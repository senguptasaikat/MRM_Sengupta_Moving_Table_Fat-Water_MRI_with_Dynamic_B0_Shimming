%% MD5
%
% Calculate MD5 checksum on a file.
%
% [MD5_HASH] = MD5(FILENAME)
%
% See also: DOWNLOAD_DATA, DOWNLOAD_DATA_ALL, DOWNLOAD_DATA_VERIFY, DOWNLOAD_DATA_VERIFY_ALL, CSV_TO_STRUCT
%
function [md5_hash] = MD5(filename)

    import java.security.*;
    import java.io.*;
    import java.math.*;

    mddigest   = MessageDigest.getInstance('MD5'); 
    filestream = FileInputStream(File(filename)); 
    digestream = DigestInputStream(filestream, mddigest);

    target_buffer = java.util.Arrays.copyOf(int8(zeros(8192,1)),8192);

    while(digestream.read(target_buffer) ~= -1) end
    
    hash = mddigest.digest();
    bi = BigInteger(1, hash);
    md5_hash = char(bi.toString(16));
    
    % prepend zeroes if necessary
    while length(md5_hash)<32,
        md5_hash = ['0', md5_hash];
    end
