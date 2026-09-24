# Archive

Writeups and the scripts that produced the published nights. A new night uses
`pipeline/` instead of these drivers.

`<GRATING> pypeit run/` was the old layout. The published curves, `.pypeit`
files, `.coadd3d` files, and `.sens` recipes now live in `results/<GRATING>/`.
The drivers, shell scripts, and logs are in `archive/<GRATING>/`. Raw frames
and the PypeIt `Calibrations/` and `Science/` trees were removed. They can be
downloaded again and reduced again.

`archive/pypeit_test/` holds the one-off checkers, template builders, and plots
used while the method was being settled. The templates those builders wrote
are in `templates/`. The scripts a new night still runs are in `pipeline/`.
