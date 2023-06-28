# climidx

A python package for deriving climate indices based on climate data.

## content:

```
 ./api:
BiC                    cfg_hwmi_other-2.yml    eur-11_import_.yml
BiH                    cfg_hwmi_other.yml      eur-11_import__.yml
cfg_hwmi_cmip-2.yml    cfg_hwmi_stkh.yml       eur-11_import___.yml
cfg_hwmi_cmip.yml      cmip5_import_cp.yml     eur-11_smhi-rca4.yml
cfg_hwmi_cordex-2.yml  cmip5_import.yml        gcm_gwls_old.yml
cfg_hwmi_cordex-e.yml  cmip5_import_.yml       gcm_gwls.yml
cfg_hwmi_cordex.yml    cmip6_smhi_len.yml      gcm_gwls_.yml
cfg_hwmi_eobs.yml      eur-11_import_eval.yml  norcp_.yml
cfg_hwmi_erai.yml      eur-11_import.yml

./climidx:
diction_.py  exe_.py  func_.py  group_.py  hwmid_.py  hwmi_.py  io_.py zzzz_.py

```

## help info for

- BiH

```
usage: BiH [-h] {expl,calc} ...

Calculate hwmi(d) based on data set on BI.nsc.liu.se

options:
  -h, --help   show this help message and exit

subcommands:
  valid subcommands

  {expl,calc}  additional help
```

```
usage: BiH expl [-h] opt

show controlfile examples

positional arguments:
  opt         options for controlfile examples: obs | cmip5 | cordex | other

options:
  -h, --help  show this help message and exit
```

```
usage: BiH calc [-h] [-l LOG] controlfile

calculate HWMId/HWMI [& even WSDI]

positional arguments:
  controlfile        yaml file with metadata

options:
  -h, --help         show this help message and exit
  -l LOG, --log LOG  logfile
```

- BiH

```
usage: BiC [-h] {calc,prnt} ...

Calcualte climate indices based on data set on BI.nsc.liu.se

options:
  -h, --help   show this help message and exit

subcommands:
  valid subcommands

  {calc,prnt}  additional help
```

```
usage: BiC prnt [-h] [opt]

Information printer

positional arguments:
  opt         options for print content: None->list all available Indices |
              INDEX->definition of INDEX | grp_X[Y...]->list Indices in group
              X[Y...]

options:
  -h, --help  show this help message and exit
```

```
usage: BiC calc [-h] [-x INDICES] [-X INDICES_EXCL] [-w GWL] [-s START]
                [-e END] [-i IDXP] [-u USER_CFG] [-p PATH_CFG] [-t TINT]
                [--y0y1 Y0Y1] [--lll LLL] [-d DOMAIN] [-o ODIR] [-g G_MTHD]
                [-n VNX] [-l LOG]
                opt

Calculator

positional arguments:
  opt                   options for dataset on BI: eobs | erai | cdx_eval;
                        cdx_eval_smhi | cdx; cdx_smhi | norcp | cmp5 | cmp6

options:
  -h, --help            show this help message and exit
  -x INDICES, --indices INDICES
                        indices to be calculated. Formats: file name (yaml):
                        read list from yaml file | indexA,indexB,indexC (no
                        space after comma) | grp_a[bc]: indices belong to a
                        [and b and c] currently available group labels: w | t
                        | p | r | c
  -X INDICES_EXCL, --indices_excl INDICES_EXCL
                        indices to not be calculated. Format:
                        indexA,indexB,indexC (no space after comma)
  -w GWL, --gwl GWL     warming levels: current | gwl15 | gwl2 | gwl25 | gwl3
                        | gwl35 | gwl4 | xx-xx | xxxx-xxxx
  -s START, --start START
                        simulation-loop start.
  -e END, --end END     simulation-loop end
  -i IDXP, --idxp IDXP  simulation-loop index. expl: 0,1,3 meaning simulation
                        #1,2,4 in the lists (specified in a yaml file) to be
                        calculated
  -u USER_CFG, --user_cfg USER_CFG
                        yaml file that stores user configuration
  -p PATH_CFG, --path_cfg PATH_CFG
                        yaml file that stores paths (simulations)
  -t TINT, --tint TINT  temporal resolution(s) of input data: mon,day
  --y0y1 Y0Y1           only for period: y0-y1
  --lll LLL             longitude/latitude limits: lo0,lo1,la0,la1
  -d DOMAIN, --domain DOMAIN
                        name of domain
  -o ODIR, --odir ODIR  directory where the results to be stored!
  -g G_MTHD, --g_mthd G_MTHD
                        method for grouping indices for a single call: None
                        (default) | 'v' | 'i'
  -n VNX, --vnx VNX     maximum number of variables: maximum input variables
                        for a single call
  -l LOG, --log LOG     exclusive log identifier

```

## what you will need:

```python
import iris
import dask.array as da
import atmos
```

