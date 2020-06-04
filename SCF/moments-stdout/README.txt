Author: Ondrej Sipr

General usage:
--------------
  moments-stdout.sc   <KKR-stdout-file>



Description:
------------

Reads stdout of a KKRSCF run and prints valence electron charge, spin
magnetic moment and orbital magnetic moment for each atomic type.

If the stdout file corresponds to a running job, it prints results for
the last "complete" cycle.

The name of the standard output file is the only argument.  It can be
specified using wildcards as long as it is unique.  So you may type
just
    moments-stdout.sc  *.out
and that's it (as long as there is only one *.out file in the current
directory). 


If there are less than 30 atomic types, it prints also the spin
magnetization and the total magnetization per unit cell.


The script is written in CSH.  It makes use of the calculator bc.  If
bc is not available, magnetization per cell is not produced (but
otherwise the script runs).

For full-potential calculation, unit cell volume cannot be calculated
so a benign message "Problems with finding average Wigner-Seitz radius
RWS ..." is produced.  It only means that the magnetization per unit
cell will not be transformed into SI units.  No real problem.




Example:
--------

$ moments-stdout.sc  ~/temp/afm-along-YZ45_SCF.out
Problems with finding average Wigner-Seitz radius RWS ...
 ------------------------------------------------
                        NOS     m_spin    m_orb  
 ------------------------------------------------
  IT=   5  Vc_5       0.0028   -0.0000  -0.00000 
  IT=   6  Vc_6       0.0028    0.0000   0.00000 
  IT=   7  Vc_7       0.5011   -0.0001  -0.00008 
  IT=   8  Vc_8       0.5011    0.0001   0.00008 
  IT=   9  W_1        5.5642   -0.0226   0.00242 
  IT=  10  W_2        5.5642    0.0226  -0.00242 
  IT=  11  W_3        5.9863    0.0010   0.00027 
  IT=  12  W_4        5.9863   -0.0010  -0.00027 
  IT=  13  W_5        6.0004    0.0151  -0.00342 
  IT=  14  W_6        6.0004   -0.0151   0.00342 
  IT=  15  W_7        5.9988    0.0120  -0.00290 
  IT=  16  W_8        5.9988   -0.0120   0.00290 
  IT=  17  W_9        5.9974    0.0601  -0.00769 
  IT=  18  W_10       5.9974   -0.0601   0.00769 
  IT=  19  W_11       5.7702    0.2626   0.01427 
  IT=  20  W_12       5.7702   -0.2626  -0.01427 
  IT=  21  Mn_1       6.8826    3.6677   0.04245 
  IT=  22  Mn_2       6.8826   -3.6677  -0.04245 
  IT=  23  Mn_3       7.9494   -2.7945  -0.20205 
  IT=  24  Vc_9       0.2938   -0.0757   0.00056 
  IT=  25  Vc_10      0.2938    0.0757  -0.00056 
  IT=  26  Vc_11      0.0024    0.0005   0.00000 
  IT=  27  Vc_12      0.0024   -0.0005  -0.00000 
  IT=  28  Vc_13      0.0000   -0.0000  -0.00000 
  IT=  29  Vc_14      0.0000    0.0000   0.00000 
Spin magnetization:   -2.794000 mu_B/cell
Total magnetization:  -2.996050 mu_B/cell
