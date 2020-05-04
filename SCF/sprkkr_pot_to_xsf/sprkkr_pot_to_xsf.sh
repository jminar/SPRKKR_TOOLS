#!/bin/bash
# Convert SPRKKR *.pot file (to be passed as first argument) into a XCrysDen -compatible *.xsf file.
#
# Notice that only actual potentials, resulting from at least one iteration of kkrscf, are supported.
#
# "SUPPRESS SYMMETRY" (from i.e. XBand button) is still assumed.
#
# No support yet for substitutional disorder i.e. CPA.
#
# A.Marmodoro, Frankfurt, TU-Darmstadt, 2016;
#              Muenchen, LMU, 2017;
#              Praha, FZU, 2020.

mag_mom_scalefactor=200  # 400

alat=$(grep ALAT $1 | tr -s " " | cut -d " " -f 2 | awk '{ print $1 * 0.52917}')

echo CRYSTAL
echo PRIMVEC

grep -A3 ALAT $1 | tail -n3 | awk -v alat=$alat '{printf "%13.9f %13.9f %13.9f\n", $2*alat, $3*alat, $4*alat}'
echo CONVVEC
grep -A3 ALAT $1 | tail -n3 | awk -v alat=$alat '{printf "%13.9f %13.9f %13.9f\n", $2*alat, $3*alat, $4*alat}'
echo PRIMCOORD

# Number of sublattices:
nsubs=$(awk '
/QBAS/ {flag=1;next}                 # Initial pattern found --> turn on the flag and read the next line
/OCCUPATION/    {flag=0}             # Final pattern found   --> turn off rhe flag
flag     {print}                     # Flag on --> print the current line
' $1 | head -n -1 | wc -l)
echo $nsubs 1

# Sublattices' origin, in cartesian coordinates.
awk '
/QBAS/ {flag=1;next}                 # Initial pattern found --> turn on the flag and read the next line
/OCCUPATION/    {flag=0}             # Final pattern found   --> turn off rhe flag
flag     {print}                     # Flag on --> print the current line
' $1 | head -n -1 | awk -v alat=$alat '{printf "%13.9f %13.9f %13.9f  \n", $2*alat, $3*alat, $4*alat}' \
> .tmp_sublattices 

awk -v nsubs=$nsubs '
/ZT/ {flag=1;next}                 # Initial pattern found --> turn on the flag and read the next line
/POTENTIAL/    {flag=0}             # Final pattern found   --> turn off rhe flag
flag  {IT_index=IT_index+1 ; if (IT_index <= nsubs) {print " ",$3," "}}               # Flag on --> print the current line, as long as all sublattices are parsed.
' $1 \
> .tmp_atomicnumbers

# Spin mag.moment value taken from the bottom of potential file (at convergence!):
#
# Magnitude rescaled by $mag_mom_scalefactor only for convenience of graphical depiction in XCrysDen (press 'f' key for 'forces'; move the figure to refresh it)
#
rm -f .tmp_mag_moments
grep "MOMENTS        QEL  NOS  SMT  OMT  HFF" -A9999  $1 | grep -v MOMENTS | grep -v === | grep -v TYPE | \
 awk -v mag_mom_scalefactor=$mag_mom_scalefactor '{ printf " %8.5f\n", $3 / mag_mom_scalefactor }' >> .tmp_mag_moments

# Rotated magnetization direction:
awk '
/IQ       MTET_Q                MPHI_Q/ {flag=1;next}                    # Initial pattern found --> turn on the flag and read the next line
/IT       MTET_T                MPHI_T                MGAM_T/ {flag=0}   # Final pattern found --> turn off the flat
flag    {print}                                                          # Flag on --> print the current line
' $1 | awk '{ print $2, $3 }' \
> .tmp_mag_angles

paste .tmp_mag_moments .tmp_mag_angles | awk '{ PI = 3.14159265 ;
 printf " %13.9f %13.9f %13.9f \n", $1 * sin( PI * $2 / 180 ) * cos ( PI * $3 / 180 ),
                                    $1 * sin( PI * $2 / 180 ) * sin ( PI * $3 / 180 ),
                                    $1 * cos ( PI * $2 / 180 ) }' \
> .tmp_mag_vectors

paste .tmp_atomicnumbers .tmp_sublattices .tmp_mag_vectors

# Clean up temporary files:
rm -f .tmp_mag_angles .tmp_mag_moments
rm -f .tmp_sublattices .tmp_atomicnumbers .tmp_mag_vectors

