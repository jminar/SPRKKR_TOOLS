Author: Ondrej Sipr

Purpose:
--------

Generates a KKR potential by picking various potentials from various
files.  Usefull to get a pre-converged potential as input for SCF run.



Description:
------------

You start with an old KKR potential file <Old-potfile>, with the
proper geometry and so on, which will serve as a template.   

Newly assembled potential file is the second argument.


How to run:

compose-kkr-potential.x   <Old-potfile>   <New-potfile>   [E_Fermi]

When asked for, supply a KKR potential filename  <pick-potfile>  where
the potential you want to pick is stored,  type index IT and RWS of
that potential, followed by IT to which this new potential should be
written and RWS of the radial grid for this new positon:

<pick-potfile>  <IT-from-where>  <RWS-from-where>  <IT-to-which>  <RWS-to-which>

Enter endo-fofile condition to finish.

Potentials for which no replacement is given will stay intact.

Optionally, a new common Fermi energy  <E_Fermi>  can be specified.
In that case all potentials are adjusted to it.


