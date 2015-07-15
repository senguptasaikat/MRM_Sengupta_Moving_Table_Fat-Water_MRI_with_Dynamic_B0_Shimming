function [ pmap id] = Threshold(vol,vol_pmap)
%%%%%%%% Thresolding the Radon Recon Field Data %%%%%%%%%%%%%

id = vol < 0.18*max(vol(:)) ;
vol_pmap(id)  = 0;

pmap = vol_pmap;