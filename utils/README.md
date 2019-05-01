# A short description of the scripts in this directory.

- `write_galcat_ddf.py` : this is effectively `write_galcat.py` with a couple of changes. First, it has extra columns so that we can do validation tests. Second, the healpixel list (which was extracted from the properties of the catalog are now set by hand to be the healpixels covering the DDF region (and neighbouring ones)
- `join_gals_SN.py` : A script that quickly joins these galaxies with the SN. The result can be written by supplying a non-None argument to the script specifying the absolute path.  
