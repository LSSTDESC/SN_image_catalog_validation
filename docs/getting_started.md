# SN host surface brightness:

## Objective: To find a distribution of the surface brightness of the host galaxy at the location of the SNe. 

- The SN we care about are on nersc in the form of a sqlite database. This might get updated, but we want to be able to run the validation tests on a similar sqlite database. The database is located [here](/global/projecta/projectdirs/lsst/groups/SSim/DC2/cosmoDC2_v1.1.4/sne_cosmoDC2_v1.1.4_MS_DDF.db). 

- To read the database on the terminal, we can go to the directory,  and use sqlite to read the top database:
```
$.open sne_cosmoDC2_v1.1.4_MS_DDF.db
$.headers on
$SELECT * FROM sne_params LIMIT 1;
>htmid_level_6|galaxy_id|c_in|mB|t0_in|x0_in|x1_in|z_in|snid_in|snra_in|sndec_in
9021|10562502111|0.0352006814105152|18.5164887395395|60772.3685148144|0.000617374435812333|1.47737398494902|0.0772781372070312|MS_10199_0|66.1155874638724|-40.8660552286036
```
So the parameters are
`htmid_level_6`, `galaxy_id`, `c_in`, `mB`, `t0_in`, `x0_in`, `x1_in`, `z_in`, `snid_in`, `snra_in`, `sndec_in`

Of these we need the `galaxy_id` which is the `id` of the galaxy in the extra-galactic catalog, but we do have to keep in mind that only values greater than `1250000000` are real galaxies (lower values should be cases where the SN are really hostless), so we restrict to the SN with `galaxy_id` greater than this value.

Next, we need to find the surface brightness at the location of the SN, given by the `snra_in` and `sndec_in` values. Ideally, this should be the surface brightness from all galaxies convolved with the PSF, etc. but that is too hard for now. We are planning to get the surface brightness contribution from a single galaxy (the host) at the location of interest.




