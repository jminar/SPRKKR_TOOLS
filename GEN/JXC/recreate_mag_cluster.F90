!
! Program to depict the magnetic ground state and show the Heisenberg 
! exchange interactions around a chosen atomic type
!
! AND/OR
!
! prepare this output in suitable form as input for Atomistic Spin 
! Dynamics pacakges like UppASD.
!
! A.Marmodoro, FZU-Praha, April 2020: initial version.
!                           May 2020: initial support for UppASD.
!
program recreate_mag_cluster
  implicit none
  integer,parameter :: iounit_potfile=654,iounit_xyzfile=656, &
    & iounit_xyzfile2=657,iounit_jmolfile=658,iounit_uppasd=659
  real(kind(0.d0)),parameter :: AU_TO_ANGSTROM = 0.52917721067D0, &
    & SMT_rescale = 200, &
    & SMT_threshold = 0.1D0, &
    & SMT_minimal = 0.1D0, &
    & min_trans = 0.20D0 ! Translucent 1.0 is invisible, this prevents it.
!
  integer IT,IQ, JT,JQ, iline, NCORT,NVALT,NSEMCORSHLT, NT,NQ, &
    & IT_CENTER, NT_CLUSTER, N_J_ij_large,i_J_ij_large,length_potfilename, &
    & IREFQ,IMQ,NOQ,nargs,max_N_J_ij_large,IT_RELEVANT
  integer,allocatable :: Z_ATOMIC(:),J_ij_target_large_from_IT(:),J_ij_endpoints(:), &
    & ITOQ(:,:),IT_ONLY_RELEVANT(:)
  real(kind(0.d0)) ALAT,ABAS(3,3),BASQ_IT_CENTER(3), &
    & QEL,NOS,OMT,HFF,MAX_SMT,MAX_J_IJ,J_ij_threshold,POS_JT(3),CONC, MAX_DR, &
    & MAGVECTOR(3),deg_to_rad
  real(kind(0.d0)),allocatable :: QBASQ(:,:),SMT(:),J_ij_value_large_from_IT(:), &
    & MTET_Q(:),MPHI_Q(:)
  character overwrite
  character(300) textline,scratchtext,headerline,potfilename,jfilename, &
    & xyzfilename,jmolfilename
  character(300) string_to_lowercase
  character(7),allocatable :: TXT_T(:)
  logical file_exists
  logical prepare_UppASD_input,SPRKKR_supported_format
! Obsolete, non-standard:
!   external iargc
!
  interface initparse_JXC_SP_SREL
    subroutine initparse_JXC_SP_SREL(jfilename,iounit_xyzfile,NT,NQ,ITOQ,IT_CENTER,Z_ATOMIC, &
    & IT_ONLY_RELEVANT, &
    & prepare_UppASD_input, &
    & QBASQ,SMT,SMT_minimal,MTET_Q,MPHI_Q,ALAT,TXT_T, &
    & J_ij_threshold,MAX_DR, &
    & NT_CLUSTER,N_J_ij_large,MAX_J_IJ, &
    & max_N_J_ij_large,J_ij_endpoints,J_ij_target_large_from_IT,J_ij_value_large_from_IT)
    implicit none
    character(300), intent(IN) :: jfilename
! No support yet for substitutional CPA
    integer,intent(IN) :: iounit_xyzfile,IT_CENTER,NT,NQ, &
      & ITOQ(1,NQ),Z_ATOMIC(NT),IT_ONLY_RELEVANT(NT)
    character(7),intent(IN) :: TXT_T(NT)
    logical,intent(IN) :: prepare_UppASD_input
    real(kind(0.d0)),intent(IN) :: QBASQ(3,NQ),SMT(NT),MTET_Q(NQ),MPHI_Q(NQ),ALAT, &
      & J_ij_threshold,MAX_DR,SMT_minimal
!
    integer,intent(OUT) :: max_N_J_ij_large,NT_CLUSTER,N_J_ij_large
    integer,allocatable,intent(OUT) :: J_ij_target_large_from_IT(:), J_ij_endpoints(:)
    real(kind(0.d0)),allocatable,intent(OUT) :: J_ij_value_large_from_IT(:)
    real(kind(0.d0)),intent(OUT) :: MAX_J_IJ
    end subroutine initparse_JXC_SP_SREL
! 
  end interface initparse_JXC_SP_SREL
!
! Obsolete, non-standard:
!   nargs = iargc()
  nargs = command_argument_count()
!
  if( (nargs.ne.4).and.(nargs.ne.5).and.(nargs.ne.6) ) then
!
    write(6,'("Expected arguments: CONVERGED potential (without SYMMETRY nor CPA) filename;")')
    write(6,'("                    JXC *.dat (for IREL=2) filename;")')
    write(6,'("                    IT_CENTER i.e. magnetic IT for which to depict Jij''s")')
    write(6,'("                    minimal J_ij_threshold (in meV) of Jij''s to be depicted")')
    write(6,'(" (optional) maximal cluster radius DR shown DR (in multiples of ALAT)")')
    write(6,'(" (optional) UppASD [case-insensitive]: input for Uppsala Atomistic Spin Dynamics package")')
    stop "Start again with 4 or 5 or 6 arguments"
!
  else
!
! Obsolete, non-standard:
!     call getarg(1,potfilename)
    call get_command_argument(1,potfilename)
!
!     call getarg(2,jfilename)
    call get_command_argument(2,jfilename)
!
!     call getarg(3,textline)
    call get_command_argument(3,textline)
    read(textline,*,ERR=5) IT_CENTER
    goto 6
5   continue
    stop "IT_CENTER must be an integer \in [1,NT]. Aborting."
6   continue
!
!     call getarg(4,textline)
    call get_command_argument(4,textline)
    read(textline,*) J_ij_threshold
!
! Safe defaults:
!
    MAX_DR = 9999D0
    prepare_UppASD_input = .FALSE.
!
    if( nargs.eq.5 ) then
!       call getarg(5,textline)
      call get_command_argument(5,textline)
      read(textline,*,ERR=9) MAX_DR
      goto 10
! Check if optional argument MAX_DR was omitted, but an atomic spin dynamics package name was specified instead:
9     continue
      textline = trim(string_to_lowercase(textline))
!
      if( index(textline,'uppasd').gt.0 ) then
        prepare_UppASD_input = .TRUE.
! For atomic spin dynamics, export the output from every atomic type as possible center:
        IT_CENTER = 1
!
        write(6,'("Ignoring input IT_CENTER, and scanning through all...")')
!
      else
        prepare_UppASD_input = .FALSE.
        write(6,'("requested Atomistic Spin Dynamics package :",a," not supported. Continuing...")') trim(textline)
      endif
    endif
10  continue
    if( nargs.eq.6 ) then
!       call getarg(6,textline)
      call get_command_argument(6,textline)
      textline = trim(string_to_lowercase(textline))
!
      if( index(textline,'uppasd').gt.0 ) then
        prepare_UppASD_input = .TRUE.
        IT_CENTER = 1
      else
        prepare_UppASD_input = .FALSE.
        write(6,'("requested Atomistic Spin Dynamics pacakge :",a," not supported. Continuing...")') textline
      endif
    endif
!
  endif
!   
  write(6,'("Reading unit cell geometry, occupancy and spin mag.moment from potfilename=",a)') trim(potfilename)
  write(6,'("Reading Jij geometry and coupling energy from jfilename=",a)') trim(jfilename)
  if( .not.prepare_UppASD_input ) then
    write(6,'("Depicting Jij for the Jij cluster around IT_CENTER=",i0)') IT_CENTER
  else
    write(6,'("Depicting Jij for the Jij cluster around EVERY 1<=IT<=NT in the unit cell")')
  endif
  write(6,'(" only if larger (in absolute value) than J_ij_threshold=",f15.8)') J_ij_threshold
  write(6,'(" and only if within MAX_DR=",f15.8," (ALAT) radius from IT_CENTER")') MAX_DR
!
! Safeguard:
  inquire(file=potfilename,EXIST=file_exists)
!
  if( .not.file_exists ) then
    write(6,'("potfilename=",a," not found, aborting.")') trim(potfilename)
    stop "check converged potential filename."
  endif
!
  open(iounit_potfile,file=potfilename,action='READ')
!
  SPRKKR_supported_format = .false.
!
! Scan for the payload:
  iline = 0
50 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=75,END=80) textline
  if( textline(1:17).eq.'PACKAGE   SPR-KKR' ) then
    read(iounit_potfile,'(a)',ERR=75,END=80) textline
    write(6,'("textline=",a)') textline
    if( textline(1:12).eq.'FORMAT     9' ) SPRKKR_supported_format = .true.
  elseif( textline(1:4).eq.'ALAT') then
    read(textline,*,ERR=75) scratchtext,alat
  elseif( textline(1:4).eq.'A(1)') then
    read(textline,*,ERR=75) scratchtext,ABAS(1:3,1)
  elseif( textline(1:4).eq.'A(2)') then 
    read(textline,*,ERR=75) scratchtext,ABAS(1:3,2)
  elseif( textline(1:4).eq.'A(3)') then
    read(textline,*,ERR=75) scratchtext,ABAS(1:3,3)
  elseif( textline(1:3).eq.'NQ ') then
    read(textline,*,ERR=75) scratchtext,NQ
  elseif( textline(1:3).eq.'NT ') then
    read(textline,*,ERR=75) scratchtext,NT
  endif
  goto 50

75 continue
  write(6,'("ERR/EOF at iline=",i0," after parsing textline=",a, &
    & " of potfilename=",a)') iline,textline,trim(potfilename)
  stop "potential file problems? Aborting"
!
80 continue
!
  if( .not.SPRKKR_supported_format ) stop "PACKAGE   SPR-KKR, FORMAT 9 not found in POTFILE"
!
  write(6,'("Parsed from: ",a,": ALAT=",f15.8)')  trim(potfilename),ALAT
  write(6,'("        and: A(1)=",3f15.8)') ABAS(1:3,1)
  write(6,'("             A(2)=",3f15.8)') ABAS(1:3,2)
  write(6,'("             A(3)=",3f15.8)') ABAS(1:3,3)
!
  if( MAX_DR.le.0.D0 ) stop "MAX_DR must be >= 0. Aborting."
!
! Safeguards:
  if( (.not.prepare_UppASD_input).and.(IT_CENTER.gt.NT) ) stop "IT_CENTER must be <= NT. Aborting."
!
!   if( NT.NE.NQ ) THEN
!     write(6,'("NQ=",i0," NT=",i0," not supported yet")') NQ,NT
!     stop "use SUPPRESS SYMMETRY / no CPA"
!   endif
!
  allocate( QBASQ(3,NQ) )
!
! Scan for the payload:
  rewind(iounit_potfile)
  headerline="        IQ        QBAS(X)               QBAS(Y)               QBAS(Z)"
  iline = 0
81 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=82,END=82) textline
  if( index(textline,headerline).gt.0 ) goto 83
  goto 81
82 continue
  write(6,'("ERR/EOF in potfilename=",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "QBAS not parsed. Aborting"
83 continue
  do JQ=1,NQ
    read(iounit_potfile,*) IQ,QBASQ(1:3,IQ)
  enddo
  write(6,'("Parsed unit cell sublattices origins.")')
!   write(6,'("QBASQ(IQ=",i0,")=",3f15.8)') (IQ,QBASQ(1:3,IQ),IQ=1,NQ)
! 
! No CPA for now: only 1 alternative on any sublattice:
  allocate( ITOQ(1,NQ) )
! Scan for the payload:
  rewind(iounit_potfile)
  headerline = "        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC"
  iline = 0
84 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=85,END=85) textline
  if( index(textline,headerline).gt.0 ) goto 86
  goto 84
85 continue
  write(6,'("ERR/EOF in potfilename=",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "ITOQ not parsed. Aborting"
86 continue
!   JT = 0
  do JQ=1,NQ
!     IO = 0
    read(iounit_potfile,*) IQ,IREFQ,IMQ,NOQ,IT,CONC
! Safeguard:
    if( NOQ.ne.1 ) then
      stop "substitutional disorder not supported yet: only 1 alternative NOQ per site IQ. Aborting"
    elseif( CONC.lt.1D0 ) then
      stop "substitutional disorder not supported yet: CONC must be 100%. Aborting."
    elseif( (IT.lt.1).or.(IT.gt.NT) ) then
      stop "inconsistent IT < 1 or > NT. Aborting."
    endif
! Scan for same IT on more than one IQ:
!     if( IT.NE.JT ) THEN
!       IO = IO + 1
!     endif
! No CPA for now: only 1 alternative on any sublattice.
    ITOQ(1,JQ) = IT
  enddo
  write(6,'("Parsed unit cell occupation by equivalent IT''s:")')
!   write(6,'("ITOQ(IO=",i0,",IQ=",i0,")=",i0)') (1,IQ,ITOQ(1,IQ),IQ=1,NQ)
!  
  allocate( Z_ATOMIC(NT) )
  allocate( TXT_T(NT) )
  allocate( IT_ONLY_RELEVANT(NT) )
! Scan for the payload:
  rewind(iounit_potfile)
  headerline="   IT     TXT_T            ZT     NCORT     NVALT    NSEMCORSHLT"
  iline = 0
90 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=91,END=91) textline
  if( index(textline,headerline).gt.0 ) goto 95
  goto 90
91 continue
  write(6,'("ERR/EOF in potfilename=",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "Z(IT) not parsed. Aborting"
!
95 continue
  do JT=1,NT
    read(iounit_potfile,*) IT,TXT_T(IT),Z_ATOMIC(IT), NCORT,NVALT,NSEMCORSHLT
    TXT_T(IT) = TRIM(TXT_T(IT))
  enddo
  write(6,'("Parsed unit cell occupancy.")')
!   write(6,'("Z(IT=",i0,")=",i0)') (IT,Z_ATOMIC(IT),IT=1,NT)

! Extract local quantization direction for the magnetization:
  headerline='        IQ       MTET_Q                MPHI_Q                MGAM_Q'
!
  allocate( MTET_Q(NQ),MPHI_Q(NQ) )
  rewind(iounit_potfile)
!
116 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=117,END=117) textline
  if( index(textline,headerline).gt.0 ) goto 118
!   write(6,*) textline
  goto 116
117 continue
  write(6,'("ERR/EOF in potfilename=",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "MPHI_Q, MTET_Q not parsed (older POTFILE version?). Aborting."
118 continue
  do JQ=1,NQ
    read(iounit_potfile,*) IQ, MTET_Q(IQ),MPHI_Q(IQ)
  enddo
! Safeguard:
  if( IQ.NE.NQ ) stop "Inconsistency in parsing angles. Aborting."
!
  write(6,'("Parsed MTET_Q, MPHI_Q.")')
  write(6,'("TXT_T(IT=",i0,")=",a, &
    & ": MTET_Q(IQ=",i0,")=",f15.8," MPHI_Q(IQ=",i0,")=",f15.8)') (ITOQ(1,IQ), &
    & TXT_T(ITOQ(1,IQ)), IQ,MTET_Q(IQ), IQ,MPHI_Q(IQ), IQ=1,NQ)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Extract the spin mag.moment (only assumed along z, not parsing angles yet):
  allocate( SMT(NT) )
!
  headerline="MOMENTS        QEL  NOS  SMT  OMT  HFF"
  iline = 0
146 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=147,END=147) textline
  if( index(textline,headerline).gt.0 ) goto 148
  goto 146
147 continue
  write(6,'("ERR/EOF in potfilename=",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "SMT not parsed. Aborting."
148 continue
  do IT=1,NT
    read(iounit_potfile,'(a)') textline 
    read(iounit_potfile,*) QEL,NOS,SMT(IT),OMT,HFF
    read(iounit_potfile,'(a)') textline
  enddo
  write(6,'("Parsed SMT.")')
  write(6,'("TXT_T(IT=",i0,")=",a," SMT(IT=",i0,")=",f15.8)') (IT, &
    & TXT_T(IT),IT,SMT(IT),IT=1,NT)
!
  close(iounit_potfile)
!
! Depict more translucent, depending on mag.moment:
!
  MAX_SMT = maxval( abs(SMT(1:NT)) )
  write(6,'("Largest spin mag.moment (absolute value in \mu_B): MAX_SMT=",f15.8)') MAX_SMT
!
! Account for the possibility that the first IT_CENTER not to be neglected
! is not the first in the listing: create a remapping, allowing to skip
! some IT's:
!
  IT_RELEVANT = 0
!
  do IT=1,NT
    if( abs(SMT(IT)).ge.SMT_minimal ) then
      IT_RELEVANT = IT_RELEVANT + 1
      IT_ONLY_RELEVANT(IT) = IT_RELEVANT
    else
!
! Safeguard: trigger an error if attempts are made to use skipped IT's:
!
      IT_ONLY_RELEVANT(IT) = 0
!
    endif
  enddo
!
  write(6,'("Due to SMT_minimal=",f15.8)') SMT_minimal
  write(6,'("IT_ONLY_RELEVANT(IT=",i0,")=",i0)') (IT,IT_ONLY_RELEVANT(IT),IT=1,NT)
!
1050 continue
! 
! Safeguard:
!
  if( abs(SMT(IT_CENTER)).lt.SMT_minimal ) then
!
    write(6,'("For TXT_T(IT_CENTER=",i0,")=",a," SMT(IT_CENTER=",i0,")=", &
      & f15.8," is too tiny, there would be no Jij to visualize.")') &
      &  IT_CENTER,TXT_T(IT_CENTER),IT_CENTER,SMT(IT_CENTER)
    write(6,'("Lower SMT_minimal if this IT_CENTER was really desired.")')
!
    if( prepare_UppASD_input ) then
!
      write(6,'("Skipping IT_CENTER=",i0," from UppASD jfile: SMT.lt.SMT_minimal=",f15.8)') IT_CENTER,SMT_minimal
!
      IT_CENTER = IT_CENTER + 1
      if( IT_CENTER.eq.NT ) goto 2000
!
      goto 1050
    else
      stop "No Jij to show if IT_CENTER has too tiny mag.moment, choose another one. Aborting."
    endif
!
  endif
!
! Validation: unit cell dump:
!
  write(xyzfilename,'("Jij-from-IT_",i0,".xyz")') IT_CENTER
  length_potfilename = len(trim(potfilename))
  xyzfilename = trim(potfilename)
!
  xyzfilename(length_potfilename-3:length_potfilename) = ".xyz"
!
  open(iounit_xyzfile,file=xyzfilename,action='WRITE')
!
  write(iounit_xyzfile,'(i0)') NQ
!
! XCrysDen defaults for forces is too long, further reduce by factor 2:
!
  write(iounit_xyzfile,'("Dump from potfilename=",a," with SMT_rescale=",f15.8)') trim(potfilename),SMT_rescale/2D0
  do IQ=1,NQ
    IT = ITOQ(1,IQ)
!
    MAGVECTOR(1) = SMT(IT) * sin( deg_to_rad(MTET_Q(IQ)) ) * cos( deg_to_rad(MPHI_Q(IQ)) )
    MAGVECTOR(2) = SMT(IT) * sin( deg_to_rad(MTET_Q(IQ)) ) * sin( deg_to_rad(MPHI_Q(IQ)) )
    MAGVECTOR(3) = SMT(IT) * cos( deg_to_rad(MTET_Q(IQ)) )
!
! XCrysDen defaults for forces is too long, further reduce by factor 2:
!
    MAGVECTOR(1:3) = MAGVECTOR(1:3) / SMT_rescale / 2D0
!    
    write(iounit_xyzfile,'(" ",i0," ",3f15.8," ",3f15.8)') Z_ATOMIC(IT),ALAT*AU_TO_ANGSTROM*QBASQ(1:3,IQ), &
      & MAGVECTOR(1:3)
  enddo
!
  close(iounit_xyzfile)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!
  write(6,'("Calling initparse_JXC_SP_SREL for IT_CENTER=",i0)') IT_CENTER
!
  call initparse_JXC_SP_SREL(jfilename,iounit_xyzfile,NT,NQ,ITOQ,IT_CENTER,Z_ATOMIC, &
    & IT_ONLY_RELEVANT, &
    & prepare_UppASD_input, &
    & QBASQ,SMT,SMT_minimal,MTET_Q,MPHI_Q,ALAT,TXT_T, &
    & J_ij_threshold,MAX_DR, &
    & NT_CLUSTER,N_J_ij_large,MAX_J_IJ, &
    & max_N_J_ij_large,J_ij_endpoints,J_ij_target_large_from_IT,J_ij_value_large_from_IT)
!
! Copy with correct header:
  write(xyzfilename,'("Jij-from-IT_",i0,".xyz")') IT_CENTER
!
  open(iounit_xyzfile2,file=xyzfilename,action='WRITE')
!
! The atom at the center has been added by hand:
  write(iounit_xyzfile2,'(i0)') NT_CLUSTER + 1
!
  write(iounit_xyzfile2,'("J_ij dump from iounit_jfile=",a, &
    & ", around TXT_T(IT_CENTER=",i0,")=",a," above SMT_threshold=", &
    & f15.8," MAX_DR=",f15.8)') trim(jfilename), &
    & IT_CENTER,TXT_T(IT_CENTER),MAX_SMT,MAX_DR
!
! Copy line-by-line from the scratch file, skipping 2 header lines:
!
  rewind(iounit_xyzfile)
!
  read(iounit_xyzfile,*) textline
  read(iounit_xyzfile,*) textline
1100 continue
  read(iounit_xyzfile,'(a)',END=1200) textline
  write(iounit_xyzfile2,'(a)') textline
  goto 1100
1200 continue
  close(iounit_xyzfile2)
!
  close(iounit_xyzfile)
!
! JMol file - open via: jmol -s file
  write(jmolfilename,'("Jij-from-IT_",i0,".jmol")') IT_CENTER
! 
  open(iounit_jmolfile,file=jmolfilename,action='WRITE')
!
  write(iounit_jmolfile,'("# Magnetic exchange couplings around TXT_T(IT_CENTER=", &
    & i0,")=",a)') IT_CENTER,TXT_T(IT_CENTER)
  write(iounit_jmolfile,'("# located within the unit cell at POS(1:3)=",f15.8)') BASQ_IT_CENTER(1:3)
  write(iounit_jmolfile,'("# above J_ij_threshold=",f15.8)') J_ij_threshold
!
  write(iounit_jmolfile,'("# ")')
  write(iounit_jmolfile,'("# Interactively customize: open ''File'' menu, then ''Console''")')
  write(iounit_jmolfile,'("# Documentation from: https://chemapps.stolaf.edu/jmol/docs/")')
! 
  write(iounit_jmolfile,'("load ''",a,"'' ")') trim(xyzfilename)
!
! The central starting point atom is always relevant, hence the +1 below too:
  J_ij_endpoints(1) = 1
!
  do i_J_ij_large=1,N_J_ij_large
!
    J_ij_endpoints(i_J_ij_large+1) = J_ij_target_large_from_IT(i_J_ij_large)
!
    write(iounit_jmolfile, &
      & '("draw arrow",i0," arrow ""> Jij = ",f6.2, &
      & " (meV)"" (atomno=",i0,") (atomno=",i0,") ")',advance='NO') &
      & i_J_ij_large, J_ij_value_large_from_IT(i_J_ij_large), &
      & 1, J_ij_target_large_from_IT(i_J_ij_large)
    write(iounit_jmolfile,'(" width .20 scale .80 ;")') 
!
! Use different colors for FM and AF couplings:
! Color names from: http://jmol.sourceforge.net/jscolors/
    if( J_ij_value_large_from_IT(i_J_ij_large).lt.0D0 ) then
      write(iounit_jmolfile,'(" color $arrow",i0," [130,130,210] ")',advance='NO') i_J_ij_large
    else
      write(iounit_jmolfile,'(" color $arrow",i0," [255,0,255] ")',advance='NO') i_J_ij_large
    endif
    write(iounit_jmolfile,'(" translucent ",f15.8)') &
      & 1D0-abs(J_ij_value_large_from_IT(i_J_ij_large))/MAX_J_IJ
  enddo
!   
  write(iounit_jmolfile,'("# Depiction of spin magnetic moment:")')
  write(iounit_jmolfile,'("vector scale 80 ;")')
  write(iounit_jmolfile,'("color vectors red ;")')
  write(iounit_jmolfile,'("vector on ;")')
  write(iounit_jmolfile,'("vectors 10 ;")')
! 
  write(iounit_jmolfile,'("# More readable atomic types:")')
  write(iounit_jmolfile,'("set fontsize 12")')
!
! Make the chosen atoms more easily visible:
  write(iounit_jmolfile,'("select Cu*; label ""Cu"" color orange translucent      0.80 ;")')
  write(iounit_jmolfile,'("select As*; label ""As"" color yellow translucent      0.80 ;")')
  write(iounit_jmolfile,'("select Mn*; label ""Mn"" color blue translucent      0.40 ;")')
  write(iounit_jmolfile,'("select Na*; label ""Na"" color green translucent      0.40 ;")')
  write(iounit_jmolfile,'("select Fe*; label ""Fe"" color red translucent      0.40 ;")')
  write(iounit_jmolfile,'("select Rh*; label ""Rh"" color gray translucent      0.40 ;")')
  write(iounit_jmolfile,'("select Xx*; color white translucent 0.85 ;")')

! Make the strongest interactions more visible:
  do i_J_ij_large=1,N_J_ij_large+1
    write(iounit_jmolfile,'("select atomno=",i0," ; color translucent 0.0 ;")') J_ij_endpoints(i_J_ij_large)
  enddo
!
!   do IT=1,NT
!     if( Z_ATOMIC(IT).eq.25 ) then
!       write(iounit_jmolfile,'("select manganese; label ""Mn"" color blue ")',advance='NO')
!     elseif( Z_ATOMIC(IT).eq.29 ) then
!       write(iounit_jmolfile,'("select copper; label ""Cu"" color orange ")',advance='NO')
!     elseif( Z_ATOMIC(IT).eq.33 ) then
!       write(iounit_jmolfile,'("select arsenic; label ""As"" color yellow ")',advance='NO')
!     elseif( Z_ATOMIC(IT).eq.0 ) then
! !       write(iounit_jmolfile,'("select Xx*; label ""Vc"" color white translucent 0.8 ")')
!       write(iounit_jmolfile,'("select Xx*; color white translucent 0.8 ;")')
!       cycle
!     endif
!     write(iounit_jmolfile,'("translucent ",f15.8," ;")') min( 1D0-abs(SMT(IT))/MAX_SMT,min_trans)
!     write(iounit_jmolfile,'("translucent ",f15.8," ;")') 1D0-abs(SMT(IT))/MAX_SMT
!   enddo  
! 
  write(iounit_jmolfile,'("spin axisAngle {1, 1, 1} 15 ;")')
  write(iounit_jmolfile,'("zoom 40 ;")')
  write(iounit_jmolfile,'("hide bonds ;")')
!
!   write(iounit_jmolfile,'("hide bonds ; delay 2 ; display bonds ; delay 2 ; loop ;")')
!
! Equatorial plane: location of the Neel antiphase boundary defect in CuMnAs 2D TB-KKR LIR NQ=128:
  write(iounit_jmolfile,'("# draw plane1 plane 400 (atomno=59) (atomno=58) (atomno=57) ")')
  write(iounit_jmolfile,'("# $plane1 brown translucent .4")')
! Tips:
  write(iounit_jmolfile,'("# Create animation via: CAPTURE ""/tmp/frame.png"" ")')
  write(iounit_jmolfile,'("# (issue a CAPTURE command again, to stop saving frames")')
  write(iounit_jmolfile,'("# and pack into *.mp4 format via: mencoder frame*png -ovc x264 -o test.mp4")')

! Coordinate axes, to be further customized:
  write(iounit_jmolfile,'("# reset ; zoom 50 ; rotate x 90 ; rotate y 20 ; axes unitcell ; axes offset .1")')
  
!   write(iounit_jmolfile,'("CAPTURE ""/tmp/video/frame.png"" ")')
!
! Zooming-in in 40 steps:
!   N_J_ij = 65 ! Misnomer
!   do i_J_ij_large=1,N_J_ij
!     write(iounit_jmolfile,'(" zoom ",i0," ; delay .2 ;")', advance='NO') 35+i_J_ij_large
!   enddo
!   write(iounit_jmolfile,'("")')
!
  write(iounit_jmolfile,'("console ;")')
  close(iounit_jmolfile)
!
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! 
! 
  if( .not.prepare_UppASD_input ) then
!
    goto 99999
!
  else
!
! Repeat the scan & jfile construction for the next IT in the unit cell:
!
    if( IT_CENTER.eq.NT ) goto 2000
    IT_CENTER = IT_CENTER + 1
!
    write(6,'("Repeating the scan for IT_CENTER=",i0," out of NT=",i0)') IT_CENTER,NT
!
    goto 1050
!
  endif
!
2000 continue
  write(6,'("Preparing additional ''inpsd.dat'' and ''momfile'' and ''posfile'' input for UppASD...")')
!
! Main input file: template for 'inpsd.dat', suitable values for BCC Fe example.
!
! Safeguard:
  inquire(file='inpsd.dat',exist=file_exists)
  if( file_exists ) then
    write(6,'("Pre-existing ""inpsd.dat"" present, overwrite ? (y/n)")')
    read(5,'(a)') overwrite
    if( overwrite.ne.'y' ) stop "Not overwriting. Aborting."
  endif
!
  open(iounit_uppasd,file='inpsd.dat',action='WRITE')
!
  xyzfilename = trim(potfilename)
!
  xyzfilename(length_potfilename-3:length_potfilename) = "    "
  xyzfilename = trim(xyzfilename)
!
  write(iounit_uppasd,'("simid ",a)') xyzfilename
!
  write(iounit_uppasd,'("ncell ",3i4," simulation dimensions")') 13,13,13
! 
  write(iounit_uppasd,'("BC     P   P   P        simulation boundary conditions (0=vacuum, P=periodic)")')
!
  write(iounit_uppasd,'("cell ",3f15.8)') ABAS(1:3,1)
  write(iounit_uppasd,'("     ",3f15.8)') ABAS(1:3,2)
  write(iounit_uppasd,'("     ",3f15.8)') ABAS(1:3,3)
!
  write(iounit_uppasd,'("Sym  0       Symmetry of lattice (0 for no, 1 for cubic, 2 for 2d cubic, 3 for hexagonal)")')
!
! Further input files:
  write(iounit_uppasd,'("posfile   ./posfile")')
  write(iounit_uppasd,'("momfile   ./momfile")')
  write(iounit_uppasd,'("exchange  ./jfile")')
!
! No disorder support for now:
  write(iounit_uppasd,'("do_ralloy 0")')
!
  write(iounit_uppasd,'("Mensemble 1")')
!
  write(iounit_uppasd,'("tseed 4499")')
!
!
!
  write(iounit_uppasd,'("maptype 2")')
!
  write(iounit_uppasd,'("SDEalgh  1   SDE solver: 1=midpoint, 2=heun, 3=heun3, 4=Heun_proper, 5=Depondt")')
!
  write(iounit_uppasd,'("Initmag  3   Initial config of moments (1=random, 2=cone, 3=spec. within momfile, 4=file)")')
!
  write(iounit_uppasd,'("ip_mode  M   Initial phase parameters")')
!
  write(iounit_uppasd,'("ip_mcanneal 1   --")')
!
  write(iounit_uppasd,'("10000 TEMP 1.00e-16 0.3     --")')
!
  write(iounit_uppasd,'("mode  M      S=SD, M=MC")')
!
  write(iounit_uppasd,'("temp  TEMP   Measurement phase parameters")')
!
  write(iounit_uppasd,'("mcNstep  50000")')
!
  write(iounit_uppasd,'("Nstep  50000")')
!
  write(iounit_uppasd,'("damping 0.1")')
!
  write(iounit_uppasd,'("timestep 1.0e-16")')
!
  write(iounit_uppasd,'("do_avrg   Y  Measure averages")')
!
  write(iounit_uppasd,'("do_cumu Y")')
  write(iounit_uppasd,'("cumu_step 50")')
  write(iounit_uppasd,'("cumu_buff 10")')
!
  write(iounit_uppasd,'("do_tottraj N Measure moments")')
!
  write(iounit_uppasd,'("tottraj_step 1000")')
!
  write(iounit_uppasd,'("plotenergy 1")')
!
  write(iounit_uppasd,'("do_sc C")')
  write(iounit_uppasd,'("do_ams Y")')
!
  write(iounit_uppasd,'("do_magdos Y")')
  write(iounit_uppasd,'("magdos_freq 200")')
  write(iounit_uppasd,'("magdos_sigma 30")')
!
  write(iounit_uppasd,'("qpoints C")')
!
  write(iounit_uppasd,'("do_stiffness Y")')
!
  write(iounit_uppasd,'("eta_max 12")')
  write(iounit_uppasd,'("eta_min 6")')
!
  write(iounit_uppasd,'("alat ",e15.8)') ALAT * AU_TO_ANGSTROM * 1D-12
  
  close(iounit_uppasd)
!
! Sublattice origins and occupancy:
!
  open(iounit_uppasd,file='posfile',action='WRITE')
!
  write(6,'("Neglecting for UppASD sublattices where abs(spin mag.moment) < SMT_threshold=",f15.8)') SMT_threshold
!
  JQ = 0
!
  loop_IQ_uppasd_posfile: do IQ=1,NQ
!
    JT = ITOQ(1,IQ)
!
    POS_JT(1:3) = QBASQ(1:3,IQ)
!
! Skip magnetically irrelevant atoms:
!
    if( abs(SMT(JT)).lt.SMT_threshold ) cycle loop_IQ_uppasd_posfile
!
! Take coordinates from the POTFILE enumeration, 
! even if some sublattices have been neglected:
!
    JQ = JQ + 1
! Possible alternative: use IT_ONLY_RELEVANT(...)
!
! No support substitutional disorder yet: writing out JT = JQ
!
    write(iounit_uppasd,'(i0," ",i0," ",3f15.8)') JQ, JQ, POS_JT(1:3)
!
  enddo loop_IQ_uppasd_posfile
!
  close(iounit_uppasd)
!
! Spin magnetic moments in the ground state:
!
  open(iounit_uppasd,file='momfile',action='WRITE')
!
  write(6,'("Neglecting for UppASD atoms with abs(spin mag.moment) < SMT_threshold=",f15.8)') SMT_threshold
!
  JQ = 0
!
  loop_IQ_uppasd_momfile: do IQ=1,NQ
!
    JT = ITOQ(1,IQ)
!
    MAGVECTOR(1) = sin( deg_to_rad(MTET_Q(IQ)) ) * cos( deg_to_rad(MPHI_Q(IQ)) )
    MAGVECTOR(2) = sin( deg_to_rad(MTET_Q(IQ)) ) * sin( deg_to_rad(MPHI_Q(IQ)) )
    MAGVECTOR(3) = cos( deg_to_rad(MTET_Q(IQ)) )
!
! Skip magnetically irrelevant atoms:
!
    if( abs(SMT(JT)).lt.SMT_threshold ) cycle loop_IQ_uppasd_momfile
!
    JQ = JQ + 1
! Possible alternative: use IT_ONLY_RELEVANT(...)
!
! No support substitutional disorder yet: writing out JT = JQ
!
! CHECK UppASD sign conventions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Present understanding: moment magnitude always in absolute value, 
! direction applied to the whole vector:
!
    if( SMT(JT).lt.0D0 ) MAGVECTOR(1:3) = -1D0 * MAGVECTOR(1:3)
!
    write(iounit_uppasd,'(i0," ",i0," ",f15.8,"   ",3f15.8)') JQ,1, abs(SMT(JT)),MAGVECTOR(1:3)
!
  enddo loop_IQ_uppasd_momfile
!
  close(iounit_uppasd)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
99999 continue
  deallocate( Z_ATOMIC )
  deallocate( IT_ONLY_RELEVANT )
  deallocate( TXT_T )
  deallocate( QBASQ )
  deallocate( ITOQ )
  deallocate( MTET_Q,MPHI_Q )
  deallocate( SMT )
!
  write(6,'("")')
  if( .not.prepare_UppASD_input ) then
    write(6,'("Trying to launch JMol in the background for visualization around IT_CENTER=", &
      & i0,"...")') IT_CENTER
    write(textline,'("jmol -s ",a)') trim(jmolfilename)
    call execute_command_line(textline, wait=.false.)
    write(6,'("")')
  else
    write(6,'("NOT trying to launch JMol. Do it by hand, if desired, via:")')
    do IT=1,NT
      write(6,'("jmol -s Jij-from-IT_",i0,".jmol")') IT
      write(6,'(" (needing also: Jij-from-IT_",i0,".xyz)")') IT
    enddo
    write(6,'("")')
  endif
  write(6,'("Successful end.")')
  write(6,'("")')
! 
end program recreate_mag_cluster
!
subroutine initparse_JXC_SP_SREL(jfilename,iounit_xyzfile,NT,NQ,ITOQ,IT_CENTER,Z_ATOMIC, &
  & IT_ONLY_RELEVANT, &
  & prepare_UppASD_input, &
  & QBASQ,SMT,SMT_minimal,MTET_Q,MPHI_Q,ALAT,TXT_T, &
  & J_ij_threshold,MAX_DR, &
  & NT_CLUSTER,N_J_ij_large,MAX_J_IJ, &
  & max_N_J_ij_large,J_ij_endpoints,J_ij_target_large_from_IT,J_ij_value_large_from_IT)
  implicit none
  integer,parameter :: iounit_jfile=655,iounit_jmolfile=658, &
    & iounit_uppasd=659
  real(kind(0.d0)),parameter :: AU_TO_ANGSTROM = 0.52917721067D0, &
    & SMT_rescale = 200, &
    & SMT_threshold = 0.1D0, &
    & RY_to_EV = 13.605693009D0, &
    & min_trans = 0.20D0 ! Translucent 1.0 is invisible, this prevents it.
!
  character(300), intent(IN) :: jfilename
! No support yet for substitutional CPA
  integer,intent(IN) :: iounit_xyzfile,IT_CENTER,NT,NQ, &
    & ITOQ(1,NQ),Z_ATOMIC(NT),IT_ONLY_RELEVANT(NT)
  character(7),intent(IN) :: TXT_T(NT)
  logical,intent(IN) :: prepare_UppASD_input
  real(kind(0.d0)),intent(IN) :: QBASQ(3,NQ),SMT(NT),MTET_Q(NQ),MPHI_Q(NQ),ALAT, &
    & J_ij_threshold,MAX_DR,SMT_minimal
!
  integer,intent(OUT) :: max_N_J_ij_large,NT_CLUSTER,N_J_ij_large
  integer,allocatable,intent(OUT) :: J_ij_target_large_from_IT(:), J_ij_endpoints(:)
  real(kind(0.d0)),allocatable,intent(OUT) :: J_ij_value_large_from_IT(:)
  real(kind(0.d0)),intent(OUT) :: MAX_J_IJ
!
  integer IT,IQ, JT,JQ, N1,N2,N3, iline, N_J_ij, &
    & i_J_ij_large,parsing_attempt
  real(kind(0.d0)) DX,DY,DZ, DR, J_meV,J_eV, BASQ_IT_CENTER(3), &
    & POS_IT(1:3),POS_JT(3), &
    & MAGVECTOR(3),deg_to_rad, &
    & POS_CELL_JT(1:3), DR_CELL, J_Ry
  character(300) textline,headerline
  character(3) TXT_IT,TXT_JT
  logical file_exists,clean_jfile
  DATA clean_jfile/.true./
!
  save clean_jfile
!
! Safeguard:
  inquire(file=jfilename,EXIST=file_exists)
  if( .not.file_exists ) then
    write(6,'("jfilename=",a," not found, aborting.")') trim(jfilename)
    stop "check Jij filename."
  endif
!
! Initial guess for the number of Jij's in the cluster around IT_CENTER:
  max_N_J_ij_large = 100
!
! Initial guess for the format of the file: SPRKKR 8.5.0, not Hutsepot:
!
  parsing_attempt = 1
!
149 continue
!
  if( prepare_UppASD_input ) then
!
! Start clean:
!
    if( clean_jfile) then
      close(iounit_uppasd,status='DELETE')
      open(iounit_uppasd,file='jfile',action='WRITE')
      clean_jfile = .false.
    else
      open(iounit_uppasd,file='jfile',action='WRITE',access='APPEND')
    endif
!
  endif
!
  open(iounit_jfile,file=jfilename,action='READ')
!
! Scan for the payload
!
  headerline = "        IT   IQ   JT   JQ   N1 N2 N3    DRX    DRY    DRZ    DR        J [meV]         J [eV]"
!
140 continue
!
  iline = 0
!
150 continue
!
  iline = iline + 1
!
  read(iounit_jfile,'(a)',ERR=200,END=200) textline
  if( index(textline,trim(headerline)).gt.0 ) goto 500
  goto 150
!
! Safeguard:
200 continue
  write(6,'("ERR/EOF in jfilename=",a," at iline=",i0, &
    & " at/after last read: textline=",a)') trim(jfilename),iline,textline
  write(6,'(" before finding header:",a)') headerline
!
! 2nd chance also failed: it wasn't 'jrs.omni format either.
!
  if( parsing_attempt.eq.1 ) then
!
    write(6,'("Second attempt: re-checking for ''jrs.omni'' Hutsepot file header...")')
!
! Try again, checking if it was in Hutsepot 'jrs.omni' alternative format:
!
    parsing_attempt = 2
!
! Initial part of the header:
!
    headerline = "site atom"
!
    rewind(iounit_jfile)
!
    goto 140
!
  elseif ( parsing_attempt.eq.2 ) then
!
    write(6,'("Third attempt: re-checking for SPRKKR dev. version file header...")')
!
    parsing_attempt = 3
!
    headerline = "XC-coupling constants J_ij"
!
  else
!
    write(6,'("Giving up after parsing_attempt=",i0)') parsing_attempt
!
    stop "file format problems? Aborting."
!
  endif
!
! Parse the cluster within which magnetic interactions are considered:
500 continue
!
! Find the position of the chosen atomic type, 
! at the center of the cluster with radius CLURAD (in units of ALAT):
!
  do IQ=1,NQ
    IT = ITOQ(1,IQ)
    if( IT_CENTER.eq.IT ) goto 501
  enddo
! Safeguard:
  stop "Specified IT_CENTER not found within NT. Aborting."
501 continue
!
! The chosen atom is not listed in the J_{ij} file, because it would mean
! that it interacts with itself. It is hence added at the beginning by hand,
! as atomno=1:
!
! Initial point location:
!
  BASQ_IT_CENTER(1:3) = QBASQ(1:3,IQ)
!
  open(iounit_xyzfile,status='SCRATCH',action='READWRITE')
!
  write(iounit_xyzfile,'(i0)') NQ
  write(iounit_xyzfile,'("Dump from iounit_jfile around TXT_T(IT_CENTER=", &
    & i0,")=",a," located at POS(1:3)=",3f15.8," within the unit cell, ", &
    & "with SMT_rescale=",f15.8)') IT_CENTER,TXT_T(IT_CENTER),BASQ_IT_CENTER(1:3), &
    & SMT_rescale
!
! Initial point mag.moment orientation:
!
  MAGVECTOR(1) = SMT(IT_CENTER) * sin( deg_to_rad(MTET_Q(IQ)) ) * cos( deg_to_rad(MPHI_Q(IQ)) )
  MAGVECTOR(2) = SMT(IT_CENTER) * sin( deg_to_rad(MTET_Q(IQ)) ) * sin( deg_to_rad(MPHI_Q(IQ)) )
  MAGVECTOR(3) = SMT(IT_CENTER) * cos( deg_to_rad(MTET_Q(IQ)) )
!
  MAGVECTOR(1:3) = MAGVECTOR(1:3) / SMT_rescale
!
  write(iounit_xyzfile,'(" ",i0," ",3f15.8,"   ",3f15.8," ",3f15.8)') Z_ATOMIC(IT_CENTER), &
    & ALAT*AU_TO_ANGSTROM*BASQ_IT_CENTER(1:3), &
    & MAGVECTOR(1:3)
!
! Add the other atoms belonging to the magnetic CLURAD-wide cluster:
!
  NT_CLUSTER = 0
  N_J_ij_large = 0
!
! Plus one, because of the manually added central atom in the cluster:
!
  allocate( J_ij_target_large_from_IT(max_N_J_ij_large), &
    &        J_ij_value_large_from_IT(max_N_J_ij_large), &
    &        J_ij_endpoints(max_N_J_ij_large+1) )
! For safety / for the max value search below:
  J_ij_target_large_from_IT(1:max_N_J_ij_large) = 0
  J_ij_value_large_from_IT(1:max_N_J_ij_large) = 0
  J_ij_endpoints(1:max_N_J_ij_large+1) = 0  
!
550 continue
!
  iline = iline + 1
!
! Parse as SPRKKR 8.5.0 format:
!
  if( parsing_attempt.eq.1 ) then
!
    read(iounit_jfile,'(a)',END=1000,ERR=600) textline
    read(textline,*,ERR=600) &
      & IT,IQ, JT,JQ, N1,N2,N3, &
      & DX,DY,DZ, DR, J_meV,J_eV
!
! Only write end-points from the chosen initial atomic type:
!
    if( IT.EQ.IT_CENTER ) then
!
! If the end point lies too far (but only if it is also not end-point of
! a signicant Jij's, skip it from visualization:
!
      if( (DR.GT.MAX_DR).and.(.not.(abs(J_meV).gt.J_ij_threshold)) ) goto 560
!
! Skip Jij's that end onto a too tiny spin mag.moment: 
! they would be negligible anyways.
!
      if( abs(SMT(JT)).lt.SMT_minimal ) goto 560

!
      NT_CLUSTER = NT_CLUSTER + 1
!
      POS_JT(1:3) = QBASQ(1:3,IQ)
!
! End point location:
!
      POS_JT(1) = POS_JT(1) + DX
      POS_JT(2) = POS_JT(2) + DY
      POS_JT(3) = POS_JT(3) + DZ
!
! End point mag.moment ground state orientation:
!
      MAGVECTOR(1) = sin( deg_to_rad(MTET_Q(JQ)) ) * cos( deg_to_rad(MPHI_Q(JQ)) )
      MAGVECTOR(2) = sin( deg_to_rad(MTET_Q(JQ)) ) * sin( deg_to_rad(MPHI_Q(JQ)) )
      MAGVECTOR(3) = cos( deg_to_rad(MTET_Q(JQ)) )
!
      write(iounit_xyzfile,'(" ",i0," ",3f15.8,"   ",3f15.8)') &
        & Z_ATOMIC(JT), &
        & ALAT*AU_TO_ANGSTROM*POS_JT(1:3), &
        & MAGVECTOR(1:3) * SMT(JT) / SMT_rescale
!
! maptype = 1:
!
!       if( prepare_UppASD_input ) write(iounit_uppasd,'(i0," ",i0," ",3f15.8, &
!         & " ",f15.8," ",f15.8)') IT_CENTER,JT,POS_JT(1:3),J_meV / RY_to_EV, DR
!
! maptype = 2:
!
      if( prepare_UppASD_input ) write(iounit_uppasd,'(i0," ",i0," ",3i7, &
        & " ",f15.8," ",f15.8)') IT_ONLY_RELEVANT(IT_CENTER),IT_ONLY_RELEVANT(JT), &
        & N1,N2,N3, J_meV / RY_to_EV, DR
!
! Store signicantly large Jij's, magnitude, and end-point within
! the cluster around IT_CENTER:
!
      if( abs(J_meV).gt.J_ij_threshold ) then
!
        N_J_ij_large = N_J_ij_large + 1
!
        if( N_J_ij_large.gt.max_N_J_ij_large ) then
          write(6,'("max_N_J_ij_large=",i0, &
            & " was too small for N_J_ij_large. Reallocating & re-parsing...")') max_N_J_ij_large
!  
          deallocate( J_ij_target_large_from_IT,J_ij_value_large_from_IT, &
            & J_ij_endpoints )
!
          close(iounit_xyzfile)
          close(iounit_jfile)
!
          if( prepare_UppASD_input ) close(iounit_uppasd)
!
          max_N_J_ij_large = max_N_J_ij_large * 2
!
          goto 149
        endif
!
! Adding +1, because of the first atom at the origin added by hand above:
!
        J_ij_target_large_from_IT(N_J_ij_large) = NT_CLUSTER + 1
        J_ij_value_large_from_IT(N_J_ij_large) = J_meV
!
        goto 560
!
      endif
!
560   continue
!
    endif
!
  elseif( parsing_attempt.eq.2 ) then
!
    write(6,'("Assuming jfilename=",a," is in ''jrs.omni'' Hutsepot format")') trim(jfilename)
!     1   Na        0.0000000000        3.9674232552        8.6447720000    5   As        0.0000000000        3.9674232552        2.8771300000           0.0000000000        0.0000000000        0.0000000000        0.0000000000    0.00000000000000000000E+00
!
    read(iounit_jfile,'(a)',END=1000,ERR=600) textline
    read(textline,*,ERR=600) &
      & IT,TXT_IT,POS_IT(1:3), &
      & JT,TXT_JT,POS_JT(1:3), &
      & POS_CELL_JT(1:3), DR_CELL, J_Ry
!
    write(6,'("DEBUG ''jrs.omni'':")')
    write(6,*) IT,TXT_IT,POS_IT(1:3), &
      & JT,TXT_JT,POS_JT(1:3), &
      & POS_CELL_JT(1:3), DR_CELL, J_Ry
!
! Only write end-points from the chosen initial atomic type:
!
    if( IT.EQ.IT_CENTER ) then
! !
! ! If the end point lies too far (but only if it is also not end-point of
! ! a signicant Jij's, skip it from visualization:
! !
!       if( (DR.GT.MAX_DR).and.(.not.(abs(J_meV).gt.J_ij_threshold)) ) goto 560
! !
!       NT_CLUSTER = NT_CLUSTER + 1
! !
!       POS_JT(1:3) = QBASQ(1:3,IQ)
! !
! ! End point location:
! !
!       POS_JT(1) = POS_JT(1) + DX
!       POS_JT(2) = POS_JT(2) + DY
!       POS_JT(3) = POS_JT(3) + DZ
! !
! ! End point mag.moment ground state orientation:
! !
!       MAGVECTOR(1) = sin( deg_to_rad(MTET_Q(JQ)) ) * cos( deg_to_rad(MPHI_Q(JQ)) )
!       MAGVECTOR(2) = sin( deg_to_rad(MTET_Q(JQ)) ) * sin( deg_to_rad(MPHI_Q(JQ)) )
!       MAGVECTOR(3) = cos( deg_to_rad(MTET_Q(JQ)) )
! !
!       write(iounit_xyzfile,'(" ",i0," ",3f15.8,"   ",3f15.8)') &
!         & Z_ATOMIC(JT), &
!         & ALAT*AU_TO_ANGSTROM*POS_JT(1:3), &
!         & MAGVECTOR(1:3) * SMT(JT) / SMT_rescale
! !
!       if( prepare_UppASD_input ) write(iounit_uppasd,'(i0," ",i0," ",3f15.8, &
!         & " ",f15.8," ",f15.8)') IT_CENTER,JT,POS_JT(1:3),J_meV / RY_to_EV, DR
! !
! ! Store signicantly large Jij's, magnitude, and end-point within
! ! the cluster around IT_CENTER:
! !
!       if( abs(J_meV).gt.J_ij_threshold ) then
! !
!         N_J_ij_large = N_J_ij_large + 1
! !
!         if( N_J_ij_large.gt.max_N_J_ij_large ) then
!           write(6,'("max_N_J_ij_large=",i0, &
!             & " was too small for N_J_ij_large. Reallocating & re-parsing...")') max_N_J_ij_large
! !  
!           deallocate( J_ij_target_large_from_IT,J_ij_value_large_from_IT, &
!             & J_ij_endpoints )
! !
!           close(iounit_xyzfile)
!           close(iounit_jfile)
! !
!           if( prepare_UppASD_input ) close(iounit_uppasd)
! !
!           max_N_J_ij_large = max_N_J_ij_large * 2
! !
!           goto 149
!         endif
! !
! ! Adding +1, because of the first atom at the origin added by hand above:
! !
!         J_ij_target_large_from_IT(N_J_ij_large) = NT_CLUSTER + 1
!         J_ij_value_large_from_IT(N_J_ij_large) = J_meV
! !
!         goto 560
! !
!       endif
! !
! 560   continue
! !
    endif
! !
    stop "not implemented yet. Aborting."
!
  elseif( parsing_attempt.eq.3 ) then
    write(6,'("Assuming jfilename=",a," is in SPRKKR dev.version format")') trim(jfilename)
    stop "not implemented yet. Aborting."
  else
    write(6,'("parsing_attempt=",i0)') parsing_attempt
    stop "not implemented yet. Aborting."
  endif
!
! Continue the Jij's scan:
!
  goto 550
!
! Safeguard:
600 continue
  write(6,'("ERR in jfilename=",a," on iline=",i0, &
    & " parsing textline=",a)') trim(jfilename),iline,textline
  stop "file format problems? Aborting."
!
1000 continue
  N_J_ij = iline - 1
  write(6,'("Whole jfilename=",a," read-in, up to iline=",i0)') trim(jfilename),N_J_ij
!
  write(6,'("Depicting with arrows Jij from TXT_T(IT_CENTER=",i0,")=",a, &
    & " i.e. cluster atomno=",i0," found larger than J_ij_threshold=", &
    & f15.8," and lying within MAX_DR=",f15.8," (ALAT)")') IT_CENTER,TXT_T(IT_CENTER), &
    & 1, J_ij_threshold,MAX_DR
  write(6,'("J_{i=",i0,",j=",i0,"}=",f15.8)') (1, &
    & J_ij_target_large_from_IT(i_J_ij_large), &
    &  J_ij_value_large_from_IT(i_J_ij_large), i_J_ij_large=1,N_J_ij_large)
!
  MAX_J_IJ = maxval( abs(J_ij_value_large_from_IT) )
  write(6,'("Strongest exchange coupling shown (in meV, absolute value):", &
    & f15.8)') MAX_J_IJ
!
!   allocate( Jij_list(6,N_J_ij) )
!   rewind(iounit_jfile)
!   goto 149
! 1100 continue
!
  close(iounit_jfile)
!
  if( prepare_UppASD_input ) close(iounit_uppasd)
!
end subroutine initparse_JXC_SP_SREL

real(kind(0.d0)) function deg_to_rad(angle)
  implicit none
  real(kind(0.d0)),parameter :: pi=3.141592653589793D0
  real(kind(0.d0)),intent(IN) :: angle
  deg_to_rad = angle * pi / 180D0
end function deg_to_rad

character(len=*) function string_to_lowercase(string)
  implicit none
  character(len=*),intent(IN) :: string
  integer i,ich
!
! Convert to lower case ASCII characters A-Z:
!
  do i = 1, len_trim(string)
    ich = ichar(string(i:i)) ! Get position of string(i:i) in ASCII set
! range A to Z:
    if( (ich.ge.65).and.(ich.le.90) ) then
! range a to z:
      string_to_lowercase(i:i) = char(ich+32)
!
    else
      string_to_lowercase(i:i) = string(i:i)
    endif
  enddo
!
end function string_to_lowercase
!
subroutine create_gnuplot_T_critical_fit()
  implicit none
!   set xlabel 'Temperature (K)'
! set ylabel 'spin mag.moment (\mu_B)'
! set title 'CuMnAs AF from SPRKKR SP-SREL'
! 
! set grid
! 
! p 'thermal.dat' u 1:2 w lp lw 4
! 
! m(x) = norm * (1 - x/Tc)**beta
! 
! # Initial Ansatz for the fit:
! Tc = 1000
! beta = .2
! fit m(x) "thermal.dat" u 1:2 via Tc,beta,norm

end subroutine create_gnuplot_T_critical_fit
