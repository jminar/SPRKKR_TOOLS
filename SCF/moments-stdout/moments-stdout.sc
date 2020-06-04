#!/bin/csh -f

if ( $#argv != 1 ) then
   echo "One argument:  SPRKKR standart output file \!"
   exit(2)
endif



# ------------------------------------------------------
# Checks whether the bc calculator is available

which bc >& /dev/null
if ( $status == 0 ) then
   set availbc = 1
else
   set availbc = 0
   echo "bc calculator not available, moments per cell will NOT be calculated"
endif





# ------------------------------------------------------
# Deals only with "complete" part of SCF standart output 

set stdfile = $1

grep -n ' E_band ' $stdfile >& /dev/null
if ( $status != 0 ) then
   echo "Problems with finding closing  E_band  target \!"
   exit(2)
endif
set complete = `grep -n ' E_band ' $stdfile | tail -1 | cut -d':' -f1`

sed ${complete},\$d  $stdfile > com_$$.tmp



# ------------------------------------------------------
# Checks consistency of file with SPRKKR results

set resfile = com_$$.tmp

grep -n 'CHARGE MISFIT' $resfile >& /dev/null
if ( $status != 0 ) then
   echo "Problems with finding initiating  CHARGE MISFIT  target \!"
   exit(2)
endif
set prvni = `grep -n 'CHARGE MISFIT' $resfile | tail -1 | cut -d':' -f1`

set posledni = `wc -l $resfile | awk '{ print $1 }'`
if ( $status != 0 ) then
   echo "Problems with counting the lines in $resfile \!"
   exit(2)
endif

if ( $posledni < 20 ) then
   echo "Suspiciously low number of lines in $resfile:   $posledni"
   exit(2)
endif

if ( ${prvni} < 2 || ${prvni} >= $posledni ) then
   echo "Strange CHARGE MISFIT target line number in $resfile:   ${prvni}"
   exit(2)
endif


# ------------------------------------------------------
# Generates auxiliary files with moments, charges etc.

sed  -n  ${prvni},${posledni}p  ${resfile}  >  konec_$$.tmp

grep  ' E=.*        IT='  konec_$$.tmp > it_$$.tmp
grep  '  P_spin '         konec_$$.tmp > hd_$$.tmp
grep  ' sum.* v+c'        konec_$$.tmp > sm_$$.tmp

set pocet = `wc -l it_$$.tmp | awk '{ print $1 }'`


# ------------------------------------------------------
# Initiates arrays for type-specific magnetic moments

set maxnt = 30

set spn = ( 0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 )
set orb = ( 0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 )
set wgt = ( 0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0  0 0 0 0 0 0 0 0 0 0 )


if ( $pocet > $maxnt ) then
   set cellmom = 0
   echo "Magnetization for the whole cell will NOT be evaluated (too many atomic types)"
else
   set cellmom = 1
endif


if ( $cellmom == 1 ) then
  grep -n '^ *type *TXTT *NL *mesh ' $resfile >& /dev/null
  if ( $status == 0 ) then
     set lintyp = `grep -n '^ *type *TXTT *NL *mesh ' $resfile | tail -1 | cut -d':' -f1`
     set formconc = 1
     goto konconc 
  endif

  grep -n '^ *type *TXTT *NL *VAL *COR *mesh ' $resfile >& /dev/null
  if ( $status == 0 ) then
     set lintyp = `grep -n '^ *type *TXTT *NL *VAL *COR *mesh ' $resfile | tail -1 | cut -d':' -f1`
     set formconc = 2
     goto konconc 
  endif

  echo "Problems with finding type TXTT NL mesh target \!"
  exit(2)
  konconc:
  
  @ zac = ${lintyp} + 1
  @ kon = ${lintyp} + ${pocet}
  
  sed  -n  ${zac},${kon}p  ${resfile}  >  vaha_$$.tmp


# Reads concentrations for each type
  
  set it = 1  
  while ( ${it} <= ${pocet} )

    if ( $formconc == 1 ) then
	set itx = "`sed  -n  ${it}p  vaha_$$.tmp | cut -c3-7`"
	if ( $itx != $it ) then
	echo "Reading of concentrations failed:  ITX, IT = ${itx}, ${it}"
	exit(2)
	endif
	
	set nat = `sed  -n  ${it}p  vaha_$$.tmp | cut -c42-45`
	set con = `sed  -n  ${it}p  vaha_$$.tmp | cut -c46-53`
    else if ( $formconc == 2 ) then
	set itx = "`sed  -n  ${it}p  vaha_$$.tmp | cut -c3-7`"
	if ( $itx != $it ) then
	echo "Reading of concentrations failed:  ITX, IT = ${itx}, ${it}"
	exit(2)
	endif
	
	set nat = `sed  -n  ${it}p  vaha_$$.tmp | cut -c47-51`
	set con = `sed  -n  ${it}p  vaha_$$.tmp | cut -c51-58`
    else
       echo "Illegal CONC format  formconc =  $formconc \!"
       exit(2)
    endif
  
    if ( $availbc == 1 ) then
      set wgt[${it}] = `echo "${nat}*${con}" | bc -l | awk '{printf "%f", $0}'`
      if ( $wgt[${it}] =~  *[A-Za-z]* ) then
        echo "wgt[${it}] must be a number while it is ${wgt[${it}]} \!"
        exit(2)
      endif
    endif
  
    @ it = ${it} + 1
  end
  
# Unit cell volume
  
  set volume = 1.0

  grep '^ *number of sites *NQ '   ${resfile} >& /dev/null
  if ( $status != 0 ) then
    echo "Problems with finding the number of sites NQ ..."
    set avermom = 0
    goto koncvol
  endif
  set nq = `grep '^ *number of sites *NQ '   ${resfile} | cut -c48-`
  if ( $nq =~  *[A-Za-z]* ) then
    echo "nq must be a number while it is $nq \!"
    exit(2)
  endif

  grep '^ *average Wigner-Seitz radius '  ${resfile} >& /dev/null
  if ( $status != 0 ) then
    echo "Problems with finding average Wigner-Seitz radius RWS ..."
    set avermom = 0
    goto koncvol
  endif
  set rws = `grep '^ *average Wigner-Seitz radius '   ${resfile} | cut -c43-`
  if ( $rws =~  *[A-Za-z]* ) then
    echo "rws must be a number while it is $rws \!"
    exit(2)
  endif

  if ( $availbc == 1 ) then
     set volume = `echo "${nq}*4/3*3.14159*${rws}^3" | bc -l | awk '{printf "%f", $0}'`
  else
     set volume = 1.0
  endif     
  if ( $volume =~  *[A-Za-z]* ) then
    echo "volume must be a number while it is $volume \!"
    exit(2)
  endif
  set avermom = 1

endif

koncvol:



# ------------------------------------------------------
# Writes the moments to stdout


echo ' ------------------------------------------------'

set it = 1
while ( ${it} <= ${pocet} )
  set aaa = "`sed  -n  ${it}p  it_$$.tmp | cut -c30-`"
  set bbb = "`sed  -n  ${it}p  hd_$$.tmp | cut -c14-23,33-42,52-61`"
  set ccc = "`sed  -n  ${it}p  sm_$$.tmp | cut -c14-23,33-42,52-61`"

  if ( $it == 1 ) then
    echo "                   ${bbb}"
     echo ' ------------------------------------------------'
  endif

  echo -n  "${aaa}"
  echo "${ccc}"

  if ( $cellmom == 1 ) then
    set spn[${it}] = `sed  -n  ${it}p  sm_$$.tmp | cut -c33-41`
    set orb[${it}] = `sed  -n  ${it}p  sm_$$.tmp | cut -c52-60`
  endif

  @ it = ${it} + 1
end



# ------------------------------------------------------
# Calculates the magnetization

if ( $cellmom == 1 && $availbc == 1 ) then
  set sumspn = 0.0
  set sumorb = 0.0
  
  set it = 1
  while ( ${it} <= ${pocet} )
    set sumspn  = `echo "${sumspn} + ${wgt[${it}]}*$spn[${it}]" | bc -l | awk '{printf "%f", $0}'`
    set sumorb  = `echo "${sumorb} + ${wgt[${it}]}*$orb[${it}]" | bc -l | awk '{printf "%f", $0}'`
    @ it = ${it} + 1
  end
    
  set sumtot  = `echo "${sumspn} + ${sumorb}" | bc -l | awk '{printf "%f", $0}'`
  
  set magspn  = `echo "625.841 * ${sumspn} / ${volume}" | bc -l | awk '{printf "%f", $0}'`
  set magtot  = `echo "625.841 * ${sumtot} / ${volume}" | bc -l | awk '{printf "%f", $0}'`

  if ( avermom == 1 ) then
    echo "Spin magnetization:   ${sumspn} mu_B/cell,   ${magspn} 10^{5} A/m"
    echo "Total magnetization:  ${sumtot} mu_B/cell,   ${magtot} 10^{5} A/m"
  else
    echo "Spin magnetization:   ${sumspn} mu_B/cell"
    echo "Total magnetization:  ${sumtot} mu_B/cell"
  endif
endif  
  
  
  
foreach f (  konec_$$.tmp  it_$$.tmp  hd_$$.tmp  sm_$$.tmp vaha_$$.tmp com_$$.tmp )
  if ( -f $f ) then
    rm $f
  endif
end

exit
