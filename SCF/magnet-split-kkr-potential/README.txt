Author: Ondrej Sipr 

Purpose:
--------

(i) Add magnetic component to an otherwise non-magnetic potential:
    use it to generate an initial spin-polarized potential for KKRSCF
    run.  Usefull if you deal with clusters or surfaces of elements
    which are non-magnetic in the bulk.

(ii) Manipulate the exchange B-field - for example, to get
     an antiferromagnetic-potential from a ferromagnetic potential.


If the potential is in a FULLPOT mode, it manipulates only the
spherical part of it.


It is a fortran program, has to be compiled (any compiler will do).



Usage:
------

magnet-split-kkr-potentials.x   <old-pot-file>  <new-pot-file>
  <IT>      [ +<addmag> | *<facmag> ]
<EOF>




Description:
------------

Two command-line arguments needed.  First must be an existing SPRKKR
potential file.  The potential is assumed to be in the
      FMT07 = '(1P,5E22.14)'
format.  If it is in different format, change the appropriate
statement at the beginning of the code.

The second argument is the name of the new (manipulated) potential.
    
Then, you will be asked to enter the IT-subscript for the potential to
be manipulated and either
   +<number>
if the B-field for atomic type IT is to be changed additively by a
value <number> or
   *<number>
if the B-field  for atomic type IT is to be multiplied by <number>.

If end-of-file condition is met (type <ctrl>+D in a terminal), reading
stops.  Potentials for atomic types not specified in this way are left intact.


The number <number> can be negative.  So, e.g.,
   2   *-1.0
is a valid instruction: it means that for atomic type 2, the B-field
will be multiplied by (-1).  Which is something you want to do to
obtained anti-ferromagnetic configuration.



See template script  template-script.sc  for an example.
