!
! Program to depict the magnetic ground state and show the Heisenberg 
! exchange interactions around a chosen atomic type
!
! AND/OR
!
! prepare this output in suitable form as input for Atomistic Spin 
! Dynamics pacakges like UppASD.
!
! AND
!
! estimate the Heisenberg energy and local effective magnetic field.
!
! A.Marmodoro, FZU-Praha, April 2020: initial version.
!                           May 2020: initial support for UppASD.
!                          June 2020: initial support for CPA potential, Heisenberg energy and effective mag.field.
!
program recreate_mag_cluster
  implicit none
  integer,parameter :: iounit_potfile=654,iounit_xyzfile=656, &
    & iounit_xyzfile2=657,iounit_jmolfile=658,iounit_uppasd=659, &
    & iounit_history=660,iounit_Heisenberg_energy=661,iounit_Bfield=662, &
! Hard-coded max. number of CPA alternative on any sublattice IQ:
    & nalts_max=10
  real(kind(0.d0)),parameter :: AU_TO_ANGSTROM = 0.52917721067D0, &
    & SMT_rescale = 400, &
    & SMT_threshold = 0.01D0, & ! Skipped magnetically irrelevant atoms within 'posfile' etc. for UppASD.
    & min_trans = 0.20D0 ! Translucent 1.0 is invisible, this prevents it.
!
  integer IT,IQ, JT,JQ, iline, NCORT,NVALT,NSEMCORSHLT, IA,NT,NQ, &
    & IT_CENTER, IQ_CENTER, NT_CLUSTER, N_J_ij_large,i_J_ij_large,length_filename, &
    & IREFQ,IMQ,iarg,nargs,max_N_J_ij_large,IT_RELEVANT,NT_RELEVANT
  integer,allocatable :: Z_ATOMIC(:),J_ij_target_large_from_IT(:),J_ij_endpoints(:), &
    & ITOQ(:,:),IT_ONLY_RELEVANT(:),NOQ(:)
  real(kind(0.d0)) ALAT,ABAS(3,3),BASQ_IT_CENTER(3), &
    & QEL,NOS,OMT,HFF,MAX_SMT,MAX_J_IJ,J_ij_threshold,POS_JT(3), MAX_DR, &
    & MAGVECTOR(3),deg_to_rad, &
    & Heisenberg_energy_from_IT_CENTER_to_every_JT, Heisenberg_energy_global_sum, &
    & Bfield_on_IT_CENTER_from_every_JT(1:3),CONC_TEST
  real(kind(0.d0)),allocatable :: QBASQ(:,:),SMT(:),J_ij_value_large_from_IT(:), &
    & MTET_Q(:),MPHI_Q(:),CONC(:)
  character overwrite
  character(300) textline,scratchtext,headerline,potfilename,jfileKKRname, &
    & xyzfilename,jmolfilename,Heisenberg_energy_filename,Bfield_filename
  character(300) string_to_lowercase
  character(7),allocatable :: TXT_T(:)
  logical file_exists,prepare_UppASD_input,SPRKKR_supported_format,chosen_xyz_cell_format
! Obsolete, non-standard:
!   external iargc
!
  interface initparse_JXC_SP_SREL
    subroutine initparse_JXC_SP_SREL( jfileKKRname, &
      & iounit_xyzfile,NT,NQ, &
      & NOQ,ITOQ,IT_CENTER,Z_ATOMIC, &
      & IT_ONLY_RELEVANT, &
      & prepare_UppASD_input, &
      & QBASQ,SMT,SMT_threshold,MTET_Q,MPHI_Q,ALAT,TXT_T, &
      & J_ij_threshold,MAX_DR, &
      & NT_CLUSTER,N_J_ij_large,MAX_J_IJ, &
      & max_N_J_ij_large,J_ij_endpoints,J_ij_target_large_from_IT, &
      & J_ij_value_large_from_IT, &
      & CONC, &
      & Heisenberg_energy_from_IT_CENTER_to_every_JT, &
      & Bfield_on_IT_CENTER_from_every_JT )
    implicit none
    character(300), intent(IN) :: jfileKKRname
! Hard-coded max. number of CPA alternative on any sublattice IQ:
    integer,parameter :: nalts_max=10
    integer,intent(IN) :: iounit_xyzfile,IT_CENTER,NT,NQ, &
      & NOQ(NQ),ITOQ(nalts_max,NQ),Z_ATOMIC(NT),IT_ONLY_RELEVANT(NT)
    character(7),intent(IN) :: TXT_T(NT)
    logical,intent(IN) :: prepare_UppASD_input
    real(kind(0.d0)),intent(IN) :: QBASQ(3,NQ),SMT(NT),MTET_Q(NQ),MPHI_Q(NQ),ALAT, &
      & J_ij_threshold,MAX_DR,SMT_threshold,CONC(NT)
!
    integer,intent(OUT) :: max_N_J_ij_large,NT_CLUSTER,N_J_ij_large
    integer,allocatable,intent(OUT) :: J_ij_target_large_from_IT(:), J_ij_endpoints(:)
    real(kind(0.d0)),allocatable,intent(OUT) :: J_ij_value_large_from_IT(:)
    real(kind(0.d0)),intent(OUT) :: MAX_J_IJ,Heisenberg_energy_from_IT_CENTER_to_every_JT, &
      & Bfield_on_IT_CENTER_from_every_JT(1:3)
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
    write(6,'("Expected arguments: CONVERGED potential (without SYMMETRY) filename;")')
    write(6,'("                    JXC *.dat (for IREL=2) filename;")')
    write(6,'("                    IT_CENTER i.e. magnetic IT for which to depict Jij''s")')
    write(6,'("                    minimal J_ij_threshold (meV) of Jij''s to be depicted")')
    write(6,'(" (optional) maximal cluster radius DR shown DR (in multiples of ALAT)")')
    write(6,'(" (optional) UppASD [case-insensitive]: input for Uppsala Atomistic Spin Dynamics package")')
    stop "Start again with 4 or 5 or 6 arguments"
!
  else
!
! Obsolete, non-standard Fortran:
!
!     call getarg(1,potfilename)
    call get_command_argument(1,potfilename)
!
    call get_command_argument(2,jfileKKRname)
!
    call get_command_argument(3,textline)
    read(textline,*,ERR=5) IT_CENTER
    goto 6
5   continue
    stop "IT_CENTER must be an integer \in [1,NT]. Aborting."
6   continue
!
    call get_command_argument(4,textline)
    read(textline,*,ERR=7) J_ij_threshold
    goto 8
7   continue
    stop "4th argument must be J_ij_threshold (meV). Aborting."
8   continue
!
! Safe defaults:
!
    MAX_DR = 99999D0
    prepare_UppASD_input = .FALSE.
!
! Safe defaults:
!
    Heisenberg_energy_global_sum = 0D0
!
    if( nargs.ge.5 ) then
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
        write(6,'("Ignoring input IT_CENTER, and scanning through all IT...")')
!
      else
!
        prepare_UppASD_input = .FALSE.
        write(6,'("requested Atomistic Spin Dynamics package :",a," not supported. Continuing...")') trim(textline)
!
      endif
!
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
! Keep track of how the tool was last used:
!
  open(iounit_history,file='last_command_line.dat')
  write(iounit_history,'("Last command line used:")')
  do iarg=1,nargs
    call get_command_argument(iarg,textline)
    write(iounit_history,'(a," ")',advance='NO') trim(textline)
  enddo
  write(iounit_history,'("")')
  close(iounit_history)
!
  write(6,'("Reading unit cell geometry, occupancy and spin mag.moment from potfilename: ",a)') trim(potfilename)
  write(6,'("Reading Jij geometry and coupling energy from jfileKKRname: ",a)') trim(jfileKKRname)
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
    write(6,'("potfilename: ",a," not found, aborting.")') trim(potfilename)
    stop "check converged potential filename."
  endif
!
! single file with all Jij's energies relative to each IT, starting empty:
!
    write(Heisenberg_energy_filename,'("Heisenberg_energy_from_each_IT.dat")')
!
    open(iounit_Heisenberg_energy,file=Heisenberg_energy_filename)
    close(iounit_Heisenberg_energy,status='DELETE')
!
    open(iounit_Heisenberg_energy,file=Heisenberg_energy_filename,action='WRITE',access='APPEND')
!
! Header:
!
    write(iounit_Heisenberg_energy,'("# IT_CENTER: TXT_T: CONC: IQ_CENTER:  Heisenberg_energy_from_IT_CENTER_to_every_JT (meV):")')
!
! single file with all effective Bfield from every JT acting on each IT_CENTER, starting empty:
!
    write(Bfield_filename,'("Bfield_on_IT_CENTER_from_every_JT.dat")')
!
    open(iounit_Bfield,file=Bfield_filename)
    close(iounit_Bfield,status='DELETE')
!
    open(iounit_Bfield,file=Bfield_filename,action='WRITE',access='APPEND')
!
    write(iounit_Bfield,'("# potfilename: ",a,", jfileKKRname: ",a)') trim(potfilename),trim(jfileKKRname)
    write(iounit_Bfield,'("# IT_CENTER:  IQ_CENTER: spin mag.moment on IT_CENTER (m_x,m_y,m_z):   ", &
      & "Bfield_on_IT_CENTER_from_every_JT (B_x,B_y,B_z):")')
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
!
! Safeguard SPRKKR potential file format version:
!
!     write(6,'("textline=",a)') textline
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
  write(6,'("ERR/EOF at iline=",i0," after parsing textline: ",a, &
    & " of potfilename: ",a)') iline,textline,trim(potfilename)
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
  allocate( QBASQ(3,NQ) )
!
! Scan for the payload:
!
  rewind(iounit_potfile)
  headerline="        IQ        QBAS(X)               QBAS(Y)               QBAS(Z)"
  iline = 0
81 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=82,END=82) textline
  if( index(textline,headerline).gt.0 ) goto 83
  goto 81
82 continue
  write(6,'("ERR/EOF in potfilename: ",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "QBAS not parsed. Aborting"
83 continue
  do JQ=1,NQ
    read(iounit_potfile,*) IQ,QBASQ(1:3,IQ)
  enddo
  write(6,'("Parsed unit cell sublattices origins.")')
!   write(6,'("QBASQ(IQ=",i0,")=",3f15.8)') (IQ,QBASQ(1:3,IQ),IQ=1,NQ)
!
  allocate( ITOQ(nalts_max,NQ),NOQ(NQ) )
  allocate( CONC(NT) )
!
! Scan for the payload:
!
  rewind(iounit_potfile)
  headerline = "        IQ     IREFQ       IMQ       NOQ  ITOQ  CONC"
  iline = 0
84 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=85,END=85) textline
  if( index(textline,headerline).gt.0 ) goto 86
  goto 84
85 continue
  write(6,'("ERR/EOF in potfilename: ",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "ITOQ not parsed. Aborting"
86 continue
!
  loop_JQ_parse_occupancy: do JQ=1,NQ
!
    read(iounit_potfile,*) IQ,IREFQ,IMQ,NOQ(JQ),(ITOQ(IA,IQ),CONC(ITOQ(IA,IQ)), &
      & IA=1,NOQ(JQ))
!
! Safeguard:
!
    if( NOQ(JQ).le.0 ) &
      & stop "NOQ < 0 alternatives per site does not make sense. Aborting"
!
    CONC_TEST = 0D0
!
    loop_IA_check_conc: do IA=1,NOQ(JQ)
!
! Safeguard:
!
      IT = ITOQ(IA,IQ)
      CONC_TEST = CONC_TEST + CONC(IT)
!
      if( (IT.lt.1).or.(IT.gt.NT) ) then
        stop "inconsistent IT < 1 or > NT. Aborting."
!
      elseif( (CONC(IT).lt.0D0).or.(CONC(IT).gt.1D0) ) then
        write(6,'("CONC(IT=",i0,")=",f15.8)') IT,CONC(IT)
        stop "CONC(IT) < 0 or > 1. Aborting."
      endif
!
    enddo loop_IA_check_conc
!
! Safeguard:
!
    if( abs(CONC_TEST-1D0).gt.0.000001D0 ) then
      write(6,'("CONC_TEST=",f15.8," on IQ=",i0)') CONC_TEST,IQ
      stop "CONC not adding up to 100%. Aborting."
    endif
!
  enddo loop_JQ_parse_occupancy
!
! Safeguard:
!
  do IQ=1,NQ
    do IA=1,NOQ(IQ)
      if( (ITOQ(IA,IQ).le.0).or.(ITOQ(IA,IQ).gt.NT) ) then
        write(6,'("ITOQ(IA=",i0,",IQ=",i0,")=",i0)') IA,IQ,ITOQ(IA,IQ)
        stop "ITOQ cannot be < 0 or > NT on any IQ. Aborting."
      endif
    enddo
  enddo
!  
  allocate( Z_ATOMIC(NT) )
  allocate( TXT_T(NT) )
  allocate( IT_ONLY_RELEVANT(NT) )
!
! Scan for the payload:
!
  rewind(iounit_potfile)
  headerline="   IT     TXT_T            ZT     NCORT     NVALT    NSEMCORSHLT"
  iline = 0
90 continue
  iline = iline + 1
  read(iounit_potfile,'(a)',ERR=91,END=91) textline
  if( index(textline,headerline).gt.0 ) goto 95
  goto 90
91 continue
  write(6,'("ERR/EOF in potfilename: ",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "Z(IT) not parsed. Aborting"
!
95 continue
  do JT=1,NT
    read(iounit_potfile,*) IT,TXT_T(IT),Z_ATOMIC(IT), NCORT,NVALT,NSEMCORSHLT
    TXT_T(IT) = TRIM(TXT_T(IT))
  enddo
!
!   write(6,'("Z(IT=",i0,")=",i0)') (IT,Z_ATOMIC(IT),IT=1,NT)
  write(6,'("Parsed unit cell occupation:")')
  do IQ=1,NQ
    do IA=1,NOQ(IQ)
      IT = ITOQ(IA,IQ)
      if( NOQ(IQ).gt.1 ) write(6,'("   ")',advance='NO')
      write(6,'("ITOQ(IA=",i0,",IQ=",i0,")=",i0," with Z_ATOMIC(IT=",i0,")=",i0, &
        & " label: ",a," CONC=",f15.8)') IA,IQ,IT,IT,Z_ATOMIC(IT),TXT_T(IT),CONC(IT)
    enddo
  enddo
!
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
  write(6,'("ERR/EOF in potfilename: ",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "MPHI_Q, MTET_Q not parsed (older POTFILE version?). Aborting."
118 continue
  do JQ=1,NQ
    read(iounit_potfile,*) IQ, MTET_Q(IQ),MPHI_Q(IQ)
  enddo
! Safeguard:
  if( IQ.NE.NQ ) stop "Inconsistency in parsing angles. Aborting."
!
  write(6,'("Parsed MTET_Q, MPHI_Q:")')
!
  do IQ=1,NQ
    do IA=1,NOQ(IQ)
!
      IT = ITOQ(IA,IQ)
!
      if( NOQ(IQ).gt.1 ) write(6,'("   ")',advance='NO')
      write(6,'("TXT_T(IT=",i0,")=",a, &
        & ": MTET_Q(IQ=",i0,")=",f15.8," MPHI_Q(IQ=",i0,")=",f15.8)') &
        & IT,TXT_T(IT), IQ,MTET_Q(IQ), IQ,MPHI_Q(IQ)
    
    enddo
  enddo
  write(6,'("")')
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
  write(6,'("ERR/EOF in potfilename: ",a," before finding headerline=", &
    & a)') trim(potfilename),headerline
  stop "SMT not parsed. Aborting."
148 continue
  do IT=1,NT
    read(iounit_potfile,'(a)') textline 
    read(iounit_potfile,*) QEL,NOS,SMT(IT),OMT,HFF
    read(iounit_potfile,'(a)') textline
  enddo
!
  close(iounit_potfile)
!
  write(6,'("Parsed SMT:")')
  write(6,'("TXT_T(IT=",i0,")=",a," SMT(IT=",i0,")=",f15.8)') (IT, &
    & TXT_T(IT),IT,SMT(IT),IT=1,NT)
  write(6,'("")')
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
    if( abs(SMT(IT)).ge.SMT_threshold ) then
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
  NT_RELEVANT = IT_RELEVANT
!
  write(6,'("")')
  write(6,'("Due to SMT_threshold=",f15.8,", NT_RELEVANT=",i0)') SMT_threshold,NT_RELEVANT
  write(6,'(" and will use the following mapping down to relevant IT only:")')
  write(6,'("IT_ONLY_RELEVANT(IT=",i0,")=",i0)') (IT,IT_ONLY_RELEVANT(IT),IT=1,NT)
  write(6,'("")')
!
1050 continue
! 
! Safeguard:
!
  if( abs(SMT(IT_CENTER)).lt.SMT_threshold ) then
!
    write(6,'("For TXT_T(IT_CENTER=",i0,")=",a," SMT(IT_CENTER=",i0,")=", &
      & f15.8," is too tiny, there would be no Jij to visualize.")') &
      &  IT_CENTER,TXT_T(IT_CENTER),IT_CENTER,SMT(IT_CENTER)
    write(6,'("Lower SMT_threshold if this IT_CENTER was really desired.")')
!
    if( prepare_UppASD_input ) then
!
      write(6,'("Skipping IT_CENTER=",i0," from UppASD jfile: SMT.lt.SMT_threshold=",f15.8)') IT_CENTER, &
        & SMT_threshold
      write(6,*) ''
!
      IT_CENTER = IT_CENTER + 1
      if( IT_CENTER.ge.NT ) goto 2000
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
  length_filename = len(trim(potfilename))
  xyzfilename = trim(potfilename)
!
  xyzfilename(length_filename-3:length_filename) = ".xsf"
!
! As (periodically repeating i.e. 3D periodic) XSF format:
!
  chosen_xyz_cell_format = .FALSE.
  call write_cell(chosen_xyz_cell_format,xyzfilename,ALAT,ABAS,NQ,NT,MTET_Q,MPHI_Q,SMT,SMT_rescale,Z_ATOMIC, &
    & ITOQ,NOQ,CONC,QBASQ)
!
! As (non-repeating i.e. 0D periodic) XYZ format:
!
  xyzfilename(length_filename-3:length_filename) = ".xyz"
!
  chosen_xyz_cell_format = .TRUE.
  call write_cell(chosen_xyz_cell_format,xyzfilename,ALAT,ABAS,NQ,NT,MTET_Q,MPHI_Q,SMT,SMT_rescale,Z_ATOMIC, &
    & ITOQ,NOQ,CONC,QBASQ)
!
!   open(iounit_xyzfile,file=xyzfilename,action='WRITE')
! !
!   write(iounit_xyzfile,'(i0)') NQ
! !
! ! XCrysDen defaults for forces is too long for depicting SMT in (mu_B): rescale by a fixed factor:
! !
!   write(iounit_xyzfile,'("Dump from potfilename: ",a," with SMT_rescale=",f15.8)') trim(potfilename),SMT_rescale
!   do IQ=1,NQ
!   
!     IT = ITOQ(1,IQ)
! !
!     MAGVECTOR(1) = SMT(IT) * sin( deg_to_rad(MTET_Q(IQ)) ) * cos( deg_to_rad(MPHI_Q(IQ)) )
!     MAGVECTOR(2) = SMT(IT) * sin( deg_to_rad(MTET_Q(IQ)) ) * sin( deg_to_rad(MPHI_Q(IQ)) )
!     MAGVECTOR(3) = SMT(IT) * cos( deg_to_rad(MTET_Q(IQ)) )
! !
! ! XCrysDen defaults for forces is too long for depicting SMT in (mu_B): rescale by a fixed factor:
! !
!     MAGVECTOR(1:3) = MAGVECTOR(1:3) / SMT_rescale
! !    
!     write(iounit_xyzfile,'(" ",i0," ",3f15.8," ",3f15.8)') Z_ATOMIC(IT),ALAT*AU_TO_ANGSTROM*QBASQ(1:3,IQ), &
!       & MAGVECTOR(1:3)
!   enddo
! !
!   close(iounit_xyzfile)
! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! ! !
!
  write(6,*) ''
  write(6,'("Calling initparse_JXC_SP_SREL for IT_CENTER=",i0)') IT_CENTER
!
  call initparse_JXC_SP_SREL(jfileKKRname,iounit_xyzfile,NT,NQ, &
    & NOQ,ITOQ,IT_CENTER,Z_ATOMIC, &
    & IT_ONLY_RELEVANT, &
    & prepare_UppASD_input, &
    & QBASQ,SMT,SMT_threshold,MTET_Q,MPHI_Q,ALAT,TXT_T, &
    & J_ij_threshold,MAX_DR, &
    & NT_CLUSTER,N_J_ij_large,MAX_J_IJ, &
    & max_N_J_ij_large,J_ij_endpoints,J_ij_target_large_from_IT, &
    & J_ij_value_large_from_IT, &
    & CONC, &
    & Heisenberg_energy_from_IT_CENTER_to_every_JT, &
    & Bfield_on_IT_CENTER_from_every_JT)
!
! Re-find the position of the chosen atomic type:
!
  do IQ=1,NQ
    do IA=1,NOQ(IQ)
      IT = ITOQ(IA,IQ)
      if( IT_CENTER.eq.IT ) goto 2501
    enddo
  enddo
!
! Safeguard:
  stop "Specified IT_CENTER not found within NT. Aborting."
2501 continue
!
! Use the found position:
!
  IQ_CENTER = IQ
!
  MAGVECTOR(1) = SMT(IT_CENTER) * sin( deg_to_rad(MTET_Q(IQ_CENTER)) ) * cos( deg_to_rad(MPHI_Q(IQ_CENTER)) )
  MAGVECTOR(2) = SMT(IT_CENTER) * sin( deg_to_rad(MTET_Q(IQ_CENTER)) ) * sin( deg_to_rad(MPHI_Q(IQ_CENTER)) )
  MAGVECTOR(3) = SMT(IT_CENTER) * cos( deg_to_rad(MTET_Q(IQ_CENTER)) )
!
  write(6,'("From IT_CENTER=",i0," on IQ_CENTER=",i0," with MAGVECTOR=",3f15.8," and CONC=",f15.8, &
    & ", Heisenberg_energy_from_IT_CENTER_to_every_JT=",f15.8, &
    & " upon adding up to MAX_DR=",f15.8)') IT_CENTER,IQ_CENTER, MAGVECTOR(1:3), &
    & CONC(IT_CENTER), &
    & Heisenberg_energy_from_IT_CENTER_to_every_JT, MAX_DR
!
  write(6,'("From IT_CENTER=",i0," on IQ_CENTER=",i0," with MAGVECTOR=",3f15.8," and CONC=",f15.8, &
    & ", Bfield_on_IT_CENTER_from_every_JT=",3f15.8, &
    & " upon adding up to MAX_DR=",f15.8)') IT_CENTER, IQ_CENTER, MAGVECTOR(1:3), &
    & CONC(IT_CENTER), &
    & Bfield_on_IT_CENTER_from_every_JT(1:3), MAX_DR
!
  write(iounit_Bfield,'(i0," ",a," CONC=",f15.8," ",i0," ",3f15.8," ",3f15.8)') IT_CENTER,TXT_T(IT_CENTER),CONC(IT_CENTER), &
    & IQ_CENTER, MAGVECTOR(1:3), &
    & Bfield_on_IT_CENTER_from_every_JT(1:3)
    
!
! Store for plotting into a single file:
  write(Heisenberg_energy_filename,'("Heisenberg_energy_from_each_IT.dat")')
!
! Payload for this IT_CENTER. If it was a non-magnetic atom, it could be almost zero
! and just bring 'noise' to the plotting. Write it out, but masked for gnuplot:
!
  if( abs(Heisenberg_energy_from_IT_CENTER_to_every_JT).gt.0.1D0 ) then
    write(iounit_Heisenberg_energy,'(i0," ",a," CONC=",f15.8," ",i0," ",f15.8)') IT_CENTER,TXT_T(IT_CENTER), &
      & CONC(IT_CENTER), IQ_CENTER, Heisenberg_energy_from_IT_CENTER_to_every_JT
  else
    write(iounit_Heisenberg_energy,'("# ",i0," ",a," CONC=",f15.8," ",i0," ",f15.8)') IT_CENTER,TXT_T(IT_CENTER), &
      & CONC(IT_CENTER), IQ_CENTER, Heisenberg_energy_from_IT_CENTER_to_every_JT
  endif
!
! Same as performing at the end:
!
! cat all_J_from_each_IT.dat | awk '{a=a+$3; print a}'
!
! PLUS concentration-weighted averaging:
!
  Heisenberg_energy_global_sum = Heisenberg_energy_global_sum + Heisenberg_energy_from_IT_CENTER_to_every_JT * CONC(IT_CENTER)
!
! Write out also the total at the end:
!
  if( IT_CENTER.eq.NT ) then
    write(iounit_Heisenberg_energy,'("# Heisenberg_energy_global_sum=",f15.8," (meV)")') &
      & Heisenberg_energy_global_sum
  endif
!
! Copy with correct header:
!
  write(xyzfilename,'("Jij-from-IT_",i0,".xyz")') IT_CENTER
!
  open(iounit_xyzfile2,file=xyzfilename,action='WRITE')
!
! The atom at the center has been added by hand:
  write(iounit_xyzfile2,'(i0)') NT_CLUSTER + 1
!
  write(iounit_xyzfile2,'("J_ij dump from jfileKKRname: ",a, &
    & ", around TXT_T(IT_CENTER=",i0,")=",a," with CONC=",f15.8,", ONLY above SMT_threshold=", &
    & f15.8," MAX_DR=",f15.8)') trim(jfileKKRname), &
    & IT_CENTER,TXT_T(IT_CENTER),CONC(IT_CENTER),MAX_SMT,MAX_DR
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
    & i0,")=",a," appearing with CONC=",f15.8)') IT_CENTER,TXT_T(IT_CENTER),CONC(IT_CENTER)
  write(iounit_jmolfile,'("# located within the unit cell at POS(1:3)=",f15.8)') BASQ_IT_CENTER(1:3)
  write(iounit_jmolfile,'("# and only depicted if above J_ij_threshold=",f15.8," (meV)")') J_ij_threshold
!
  write(iounit_jmolfile,'("# ")')
  write(iounit_jmolfile,'("# Interactively customize: open ''File'' menu, then ''Console''")')
  write(iounit_jmolfile,'("# Documentation from: https://chemapps.stolaf.edu/jmol/docs/")')
!
! Keep track of how the tool was last used:
!
  write(iounit_jmolfile,'("# Last command line used to create this file:")')
  write(iounit_jmolfile,'("# ")',advance='NO') 
  do iarg=1,nargs
    call get_command_argument(iarg,textline)
    write(iounit_jmolfile,'(a," ")',advance='NO') trim(textline)
  enddo
  write(iounit_jmolfile,'("")')
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
  write(iounit_jmolfile,'("select Ti*; label ""Ti"" color brown translucent      0.40 ;")')
  write(iounit_jmolfile,'("select Al*; label ""Al"" color pink translucent      0.40 ;")')
  write(iounit_jmolfile,'("select Cr*; label ""Cr"" color gray translucent      0.40 ;")')
  write(iounit_jmolfile,'("select V*; label ""V"" color white translucent      0.40 ;")')
!
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
  write(iounit_jmolfile,'("hide bonds ;")')
  write(iounit_jmolfile,'("# spin axisAngle {1, 1, 1} 15 ;")')
!   write(iounit_jmolfile,'("zoom 40 ;")')
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
  write(iounit_jmolfile,'("reset ; zoom 70 ; rotate z -90 ; rotate y 90 ; rotate x 4")')
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
  xyzfilename(length_filename-3:length_filename) = "    "
  xyzfilename = trim(xyzfilename)
!
  write(iounit_uppasd,'("simid ",a)') xyzfilename
!
  write(iounit_uppasd,'("ncell ",3i4," simulation dimensions")') 19,19,19
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
!
  write(iounit_uppasd,'("do_ralloy 0")')
  do IQ=1,NQ
!     if( NOQ(IQ).gt.1 ) stop "1 CPA addition ongoing. Aborting"
  enddo
!
  write(iounit_uppasd,'("Mensemble 1")')
!
  write(iounit_uppasd,'("tseed 4499")')
!
! Ensure C(artesian) rather than D(irect) coordinates:
!
  write(iounit_uppasd,'("posfiletype C")')
!
  write(iounit_uppasd,'("maptype 2")')
!
  write(iounit_uppasd,'("SDEalgh  1   SDE solver: 1=midpoint, 2=heun, 3=heun3, 4=Heun_proper, 5=Depondt")')
!
  write(iounit_uppasd,'("Initmag  3   Initial config of moments (1=random, 2=cone, 3=spec. within momfile, 4=file)")')
!
  write(iounit_uppasd,'("ip_mode  M   Initial phase parameters ")',advance='NO')
  write(iounit_uppasd,'("(S=SD, M=Monte Carlo, H=Heat bath Monte Carlo,N=none)")')
!
  write(iounit_uppasd,'("ip_mcanneal 1   --")')
!
  write(iounit_uppasd,'("10000 TEMP 1.00e-16 0.3     --")')
!
  write(iounit_uppasd,'("mode  M      S=SD, M=MC")')
!
  write(iounit_uppasd,'("temp  TEMP   Measurement phase parameters")')
!
  write(iounit_uppasd,'("mcNstep  150000")')
!
  write(iounit_uppasd,'("Nstep  150000")')
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
! Resolved output:
!
  write(iounit_uppasd,'("do_proj_avrg   Y     Magnetizations of sublattices")')
  write(iounit_uppasd,'("do_cumu_proj   Y")')
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
    loop_IA_uppasd_posfile: do IA=1,NOQ(IQ)
!
      JT = ITOQ(IA,IQ)
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
!       if( IA.gt.1 ) stop "2 CPA addition ongoing. Aborting" ! Safeguard.
!
      write(iounit_uppasd,'(i0," ",i0," ",3f15.8)') JQ, JQ, POS_JT(1:3)
!
    enddo loop_IA_uppasd_posfile
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
    loop_IA_uppasd_momfile: do IA=1,NOQ(IQ)
!
      JT = ITOQ(IA,IQ)
!
!       if( IA.gt.1 ) stop "3 CPA addition ongoing. Aborting" ! Safeguard.
!
! Skip magnetically irrelevant alternative occupants:
!
      if( abs(SMT(JT)).lt.SMT_threshold ) cycle loop_IA_uppasd_momfile
      
!       if( IA.gt.1 ) stop "4 CPA addition ongoing. Read UppASD example first. Aborting"
!
! IQ-only dependent, but repeated with tiny overhead for each IA, for sign handling
! reasons according to SMT(JT) below:
!
      MAGVECTOR(1) = sin( deg_to_rad(MTET_Q(IQ)) ) * cos( deg_to_rad(MPHI_Q(IQ)) )
      MAGVECTOR(2) = sin( deg_to_rad(MTET_Q(IQ)) ) * sin( deg_to_rad(MPHI_Q(IQ)) )
      MAGVECTOR(3) = cos( deg_to_rad(MTET_Q(IQ)) )
! CHECK UppASD sign conventions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! Present understanding: moment magnitude always in absolute value, 
! direction applied to the whole vector:
!
      if( SMT(JT).lt.0D0 ) MAGVECTOR(1:3) = -1D0 * MAGVECTOR(1:3)
!
!
! No support substitutional disorder yet: writing out JT = JQ
!
!     if( NOQ(IQ).gt.1 ) stop "5 CPA addition ongoing. Aborting" ! Safeguard.
!
! Skip magnetically irrelevant positions:
! Possible alternative: use IT_ONLY_RELEVANT(...)
! ! ! !     if( maxval(abs(SMT(1:NT))).ge.SMT_threshold ) then
!
      JQ = JQ + 1
!
      write(iounit_uppasd,'(i0," ",i0," ",f15.8,"   ",3f15.8)') JQ,1, abs(SMT(JT)),MAGVECTOR(1:3)
! 4-11-2020 revision for AF:
!       write(iounit_uppasd,'(i0," ",i0," ",f15.8,"   ",3f15.8)') JQ,1, SMT(JT),abs(MAGVECTOR(1:3))
!
! ! ! !     endif
    enddo loop_IA_uppasd_momfile
!
  enddo loop_IQ_uppasd_momfile
!
  close(iounit_uppasd)
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
99999 continue
!
  call get_critical_temperature_mean_field(jfileKKRname,NT,NQ,NOQ,ITOQ,CONC)
!
  deallocate( Z_ATOMIC )
  deallocate( QBASQ )
!
  deallocate( ITOQ,NOQ )
!
  deallocate( MTET_Q,MPHI_Q )
!
  close(iounit_Heisenberg_energy)
!
! Support GNUPLOT script to plot the Heisenberg energy:
!
  length_filename = len(trim(Heisenberg_energy_filename))
!
  Heisenberg_energy_filename(length_filename-3:) = ".gnuplot"

  open(iounit_Heisenberg_energy,file=Heisenberg_energy_filename,action='WRITE')

  write(iounit_Heisenberg_energy,'("set xlabel ''atom index IQ''; set xtics 7,4; set xrange [21:77]")')
! set xtics 3,4
! set x2tics 3,4
  write(iounit_Heisenberg_energy,'("set ylabel ''E^{Heisenberg} on site IQ (meV)''")')
!   write(iounit_Heisenberg_energy,'("set y2label ''e_z on site IQ''")')
!   write(iounit_Heisenberg_energy,'("unset y2tics")')
!
  write(iounit_Heisenberg_energy,'("set key bottom center")')
  write(iounit_Heisenberg_energy,'("set zeroaxis")')
!
!   write(iounit_Heisenberg_energy,'("set title ''stability assessment: mag.moment (circles) ", &
!     & "parallel or antiparallel to effective Bfield from neighbors?''")')

  write(iounit_Heisenberg_energy,'("set ytics nomirror")')
  write(iounit_Heisenberg_energy,'("set y2label ''E^{Heisenberg} sum up to IQ (meV)'' ; set y2tics nomirror")')
! 
  write(iounit_Heisenberg_energy,'("")')
  write(iounit_Heisenberg_energy,'("cumulative_sum(x,partial)=(partial = partial + x, partial)")')
  write(iounit_Heisenberg_energy,'("")')
  write(iounit_Heisenberg_energy,'("IQBOT=",i0," ; IQTOP=",i0)') 27,71 ! 1,NQ
  write(iounit_Heisenberg_energy,'("")')
!   write(iounit_Heisenberg_energy,'("p ''Heisenberg_energy_from_each_IT.dat'' u 5:6 w lp lw 3, \")')
  write(iounit_Heisenberg_energy,'("partial=0 ; p ''Heisenberg_energy_from_each_IT.dat'' u 5:6 w lp lw 3, \")')
!   write(iounit_Heisenberg_energy,'("'''' u 5:(cumulative_sum($6,partial)) w lp lw 3 axes x1y2")')
  write(iounit_Heisenberg_energy,'("'''' u 5:( ( ($5 >= IQBOT) && ($5 <= IQTOP)) ? \")')
  write(iounit_Heisenberg_energy,'("   cumulative_sum($6,partial) : NaN) w lp lw 3 axes x1y2")')
  write(iounit_Heisenberg_energy,'("")')
!
  write(iounit_Heisenberg_energy,'("print ''Total between IQBOT='',IQBOT,'' \")')
  write(iounit_Heisenberg_energy,'("   and IQTOP='',IQTOP,'' : '',partial,'' (meV)''")')
!
  close(iounit_Heisenberg_energy)
!
  close(iounit_Bfield)
!
  if( prepare_UppASD_input ) then
!
    write(6,'("Preparing queue submission script and trial temperature listing for UppASD...")')
!
! Submission scripts and trial temperatures range:
!
    call create_T_scan_scripts()
!
! Support GNUPLOT script to fit and plot the magnetization vs. temperature 
! from Monte Carlo UppASD scan.
!
! Helper parser:
!
    call parse_IT_resolved_MC_averages(NT_RELEVANT)
!
! Gnuplot for numerical fit & results ploting:
!
    call create_gnuplot_T_critical_fit(maxval(abs(SMT(1:NT))))
!
  endif
!
  deallocate( SMT )
!
! Support GNUPLOT script to plot the effective mag.field:
!
  length_filename = len(trim(Bfield_filename))
!
  Bfield_filename(length_filename-3:) = ".gnuplot"

  open(iounit_Bfield,file=Bfield_filename,action='WRITE')

  write(iounit_Bfield,'("set xlabel ''atom index IQ''; set xtics 7,4; set xrange [21:77]")')
! set xtics 3,4
! set x2tics 3,4
  write(iounit_Bfield,'("set ylabel ''B_z^{eff.} perceived on site IQ''")')
  write(iounit_Bfield,'("set y2label ''e_z on site IQ''")')
  write(iounit_Bfield,'("unset y2tics")')
!
  write(iounit_Bfield,'("set key bottom center")')
!
  write(iounit_Bfield,'("set title ''stability assessment: mag.moment (circles) ", &
    & "parallel or antiparallel to effective Bfield from neighbors?''")')
!
! set ytics 17
! 
  write(iounit_Bfield,'("set grid")')
!
  write(iounit_Bfield,'("p ''Bfield_on_IT_CENTER_from_every_JT.dat'' u 5:11 w lp lw 3 pt 5 ps 1 lc 0", &
    & " t ''eff. Bfield from all Jij'', \")')
  write(iounit_Bfield,'(" '''' u 5:(abs($8) >= 1 ? ($8 >= ",f15.8)',advance='NO') SMT_threshold
  write(iounit_Bfield,'(" ? 1 : -1) : NaN)  axes x1y2 w lp ps 3 pt 7", &
    & " lc rgb ''brown'' t ''Mn atoms spin mag.mom. sign''")')
!
  close(iounit_Bfield)
!
  write(6,'("IT_CENTER=",i0," i.e. ",a," with CONC=",f15.8," on IQ_CENTER=",i0, &
    & ": Heisenberg_energy_global_sum=",f15.8)') IT_CENTER, &
    & TXT_T(IT_CENTER),CONC(IT_CENTER), IQ_CENTER, &
    & Heisenberg_energy_global_sum
!
  write(6,'("   Bfield_on_IT_CENTER_from_every_JT=",3f15.8)') Bfield_on_IT_CENTER_from_every_JT(1:3)
!
  write(6,'("")')
!
  deallocate( TXT_T )
  deallocate( CONC )
!
  if( .not.prepare_UppASD_input ) then
!
    write(6,'("Trying to launch JMol in the background for visualization around IT_CENTER=", &
      & i0,"...")') IT_CENTER
!
    write(textline,'("jmol -s ",a)') trim(jmolfilename)
    call execute_command_line(textline, wait=.false.)
!
    write(6,'("")')
  else
!
    write(6,'("NOT trying to launch JMol. Do it by hand, if desired, via:")')
!
    do IT=1,NT
      if( IT_ONLY_RELEVANT(IT).gt.0 ) then
        write(6,'("jmol -s Jij-from-IT_",i0,".jmol")') IT
        write(6,'(" (needing also: Jij-from-IT_",i0,".xyz)")') IT
      endif
    enddo
    write(6,'("")')
  endif
!
  deallocate( IT_ONLY_RELEVANT )
  write(6,'("Successful end.")')
  write(6,'("")')
! 
end program recreate_mag_cluster
!
subroutine initparse_JXC_SP_SREL(jfileKKRname,iounit_xyzfile,NT,NQ, &
  & NOQ,ITOQ,IT_CENTER,Z_ATOMIC, &
  & IT_ONLY_RELEVANT, &
  & prepare_UppASD_input, &
  & QBASQ,SMT,SMT_threshold,MTET_Q,MPHI_Q,ALAT,TXT_T, &
  & J_ij_threshold,MAX_DR, &
  & NT_CLUSTER,N_J_ij_large,MAX_J_IJ, &
  & max_N_J_ij_large,J_ij_endpoints,J_ij_target_large_from_IT,J_ij_value_large_from_IT, &
  & CONC, &
  & Heisenberg_energy_from_IT_CENTER_to_every_JT, Bfield_on_IT_CENTER_from_every_JT )
  implicit none
  integer,parameter :: iounit_KKR_jfile=655,iounit_jmolfile=658, &
    & iounit_uppasd=659, &
! Hard-coded max. number of CPA alternative on any sublattice IQ:
    & nalts_max=10
  real(kind(0.d0)),parameter :: AU_TO_ANGSTROM = 0.52917721067D0, &
    & SMT_rescale = 100, &
    & RY_to_EV = 13.605693009D0, &
    & min_trans = 0.20D0 ! Translucent 1.0 is invisible, this prevents it.
!
  character(300), intent(IN) :: jfileKKRname
! No support yet for substitutional CPA
  integer,intent(IN) :: iounit_xyzfile,IT_CENTER,NT,NQ, &
    & ITOQ(nalts_max,NQ),NOQ(NQ),Z_ATOMIC(NT),IT_ONLY_RELEVANT(NT)
  character(7),intent(IN) :: TXT_T(NT)
  logical,intent(IN) :: prepare_UppASD_input
  real(kind(0.d0)),intent(IN) :: QBASQ(3,NQ),SMT(NT),MTET_Q(NQ),MPHI_Q(NQ),ALAT, &
    & J_ij_threshold,MAX_DR,SMT_threshold,CONC(NT)
!
  integer,intent(OUT) :: max_N_J_ij_large,NT_CLUSTER,N_J_ij_large
  integer,allocatable,intent(OUT) :: J_ij_target_large_from_IT(:), J_ij_endpoints(:)
  real(kind(0.d0)),allocatable,intent(OUT) :: J_ij_value_large_from_IT(:)
  real(kind(0.d0)),intent(OUT) :: MAX_J_IJ,Heisenberg_energy_from_IT_CENTER_to_every_JT, &
    & Bfield_on_IT_CENTER_from_every_JT(1:3)
!
  integer IT,IA,IQ, JT,JQ, N1,N2,N3, iline, N_J_ij, &
    & i_J_ij_large,parsing_attempt,IQ_CENTER
  real(kind(0.d0)) DX,DY,DZ, DR, J_meV,J_eV, BASQ_IT_CENTER(3), &
    & POS_IT(1:3),POS_JT(3), &
    & MAGVECTOR_JT(3),MAGVECTOR_IT_CENTER(3), &
    & deg_to_rad, &
    & J_Ry
  character(300) textline,headerline
  character(3) TXT_IT,TXT_JT
  logical file_exists,clean_jfile,firstmessage
  DATA clean_jfile/.true./
  save clean_jfile
!
!   clean_jfile = .true.
!
! Safeguard:
  inquire(file=trim(jfileKKRname),EXIST=file_exists)
  if( .not.file_exists ) then
    write(6,'("jfileKKRname: ",a," not found, aborting.")') trim(jfileKKRname)
    stop "check Jij filename."
  endif
!
  write(6,'("Parsing for IT_CENTER=",i0," with CONC=",f15.8)') IT_CENTER,CONC(IT_CENTER)
!
  firstmessage = .true.
!
! Initial guess for the number of Jij's in the cluster around IT_CENTER:
  max_N_J_ij_large = 1000
!
! Initial guess for the format of the file: SPRKKR 8.5.0, not Hutsepot:
!
  parsing_attempt = 1
!
! Allow to re-start with enlarged arrays:
!
149 continue
!
  Heisenberg_energy_from_IT_CENTER_to_every_JT = 0D0
!
  Bfield_on_IT_CENTER_from_every_JT(1:3) = 0D0
!
  if( prepare_UppASD_input ) then
!
! Start clean:
!
    if( clean_jfile ) then
      close(iounit_uppasd,status='DELETE')
      open(iounit_uppasd,file='jfile',action='WRITE')
      clean_jfile = .false.
    else
      open(iounit_uppasd,file='jfile',action='WRITE',access='APPEND')
    endif
!
  endif
!
  open(iounit_KKR_jfile,file=jfileKKRname,action='READ',ERR=20)
! Safeguard
  goto 30
20 stop "Error opening jfileKKRname. Aborting."
30 continue
!
! Scan for the payload:
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
  read(iounit_KKR_jfile,'(a)',ERR=200,END=200) textline
  if( index(textline,trim(headerline)).gt.0 ) goto 500
  goto 150
!
! Safeguard:
200 continue
  write(6,'("ERR/EOF in jfileKKRname: ",a," at iline=",i0, &
    & " at/after last read: textline=",a)') trim(jfileKKRname),iline,textline
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
    rewind(iounit_KKR_jfile)
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
!
500 continue
!
! Find the position of the chosen atomic type, 
! at the center of the cluster with radius CLURAD (in units of ALAT):
!
  do IQ=1,NQ
    do IA=1,NOQ(IQ)
      IT = ITOQ(IA,IQ)
      if( IT_CENTER.eq.IT ) goto 501
    enddo
  enddo
!
! Safeguard:
  stop "Specified IT_CENTER not found within NT. Aborting."
501 continue
!
! Location of IT_CENTER from the scan above:
  IQ_CENTER = IQ
!
  BASQ_IT_CENTER(1:3) = QBASQ(1:3,IQ_CENTER)
!
! The chosen atom is not listed in the J_{ij} file, because it would mean
! that it interacts with itself. It is hence added at the beginning by hand,
! as atomno=1:
!
  open(iounit_xyzfile,status='SCRATCH',action='READWRITE')
!
  write(iounit_xyzfile,'(i0)') NQ
  write(iounit_xyzfile,'("Dump from iounit_KKR_jfile around TXT_T(IT_CENTER=", &
    & i0,")=",a," located at POS(1:3)=",3f15.8," as IQ=",i0," within the unit cell, ", &
    & "with SMT_rescale=",f15.8)') IT_CENTER,TXT_T(IT_CENTER),BASQ_IT_CENTER(1:3), &
    & IQ_CENTER, SMT_rescale
!Debug:
!   write(6,'("Dump from iounit_KKR_jfile around TXT_T(IT_CENTER=", &
!     & i0,")=",a," located at POS(1:3)=",3f15.8," as IQ=",i0," within the unit cell, ", &
!     & "with SMT_rescale=",f15.8)') IT_CENTER,TXT_T(IT_CENTER),BASQ_IT_CENTER(1:3), IQ, &
!     & SMT_rescale
!
! Initial point mag.moment orientation:
!
  MAGVECTOR_IT_CENTER(1) = SMT(IT_CENTER) * sin( deg_to_rad(MTET_Q(IQ_CENTER)) ) * cos( deg_to_rad(MPHI_Q(IQ_CENTER)) )
  MAGVECTOR_IT_CENTER(2) = SMT(IT_CENTER) * sin( deg_to_rad(MTET_Q(IQ_CENTER)) ) * sin( deg_to_rad(MPHI_Q(IQ_CENTER)) )
  MAGVECTOR_IT_CENTER(3) = SMT(IT_CENTER) * cos( deg_to_rad(MTET_Q(IQ_CENTER)) )
!
  MAGVECTOR_IT_CENTER(1:3) = MAGVECTOR_IT_CENTER(1:3) / SMT_rescale
!
! Format: index, cartesian coordinates (in \AA), spin mag.moment components (in rescaled \mu_B)
!
  write(iounit_xyzfile,'(" ",i0," ",3f15.8,"   ",3f15.8," ",3f15.8)') Z_ATOMIC(IT_CENTER), &
    & ALAT*AU_TO_ANGSTROM*BASQ_IT_CENTER(1:3), &
    & MAGVECTOR_IT_CENTER(1:3)
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
    if( firstmessage ) then
      firstmessage = .false.
      write(6,'(" Assuming jfileKKRname: ",a," is in SPRKKR 8.5.0 format")') trim(jfileKKRname)
    endif
!
    read(iounit_KKR_jfile,'(a)',END=1000,ERR=600) textline
    read(textline,*,ERR=600) &
      & IT,IQ, JT,JQ, N1,N2,N3, &
      & DX,DY,DZ, DR, J_meV,J_eV
! Debug:
!     write(6,*) textline  
!
! In SPRKKR, the sign of Heisenberg exhange energy 
! is always referring to a FM coupling. Swap the sign
! if the two atoms are instead in a AF ground state:
!
    J_meV = sign(1D0,SMT(IT))*sign(1D0,SMT(JT)) * J_meV
    J_eV = sign(1D0,SMT(IT))*sign(1D0,SMT(JT)) * J_eV
!
! Only write end-points from the chosen initial atomic type:
!
    if( IT.EQ.IT_CENTER ) then
!
! If the contribution is smaller than the chosen cutoff, skip it:
!
      if( abs(J_meV).lt.J_ij_threshold ) goto 560
!
! If the end point lies too far, skip it from visualization:
!
      if( DR.GT.MAX_DR ) goto 560
!
! If the end point lies too far (but only if it is also not end-point of
! a signicant Jij's, skip it from visualization:
!       if( (DR.GT.MAX_DR).and.(.not.(abs(J_meV).gt.J_ij_threshold)) ) goto 560
!
! Skip Jij's that end onto a too tiny spin mag.moment: 
! they would be negligible anyways.
!
      if( abs(SMT(JT)).lt.SMT_threshold ) goto 560
!
      NT_CLUSTER = NT_CLUSTER + 1
!
! Safeguard:
!
      if( IQ.ne.IQ_CENTER ) stop "Starting atom at IQ <> IQ_CENTER. Aborting."
!
! End atom location, given as offset from the starting atom:
!
      POS_JT(1) = QBASQ(1,IQ) + DX
      POS_JT(2) = QBASQ(2,IQ) + DY
      POS_JT(3) = QBASQ(3,IQ) + DZ
!
! End point mag.moment ground state orientation on JT. In substitutional CPA cases, 
! the orientation is here given position-wise, i.e. indexed by the common JQ, not JT.
!
      MAGVECTOR_JT(1) = sin( deg_to_rad(MTET_Q(JQ)) ) * cos( deg_to_rad(MPHI_Q(JQ)) )
      MAGVECTOR_JT(2) = sin( deg_to_rad(MTET_Q(JQ)) ) * sin( deg_to_rad(MPHI_Q(JQ)) )
      MAGVECTOR_JT(3) = cos( deg_to_rad(MTET_Q(JQ)) )
!
! Accumulate the effective magnetic field perceived by IT_CENTER, 
! as generated by all the other mag.moments JT located at JQ:
!
      Bfield_on_IT_CENTER_from_every_JT(1:3) = Bfield_on_IT_CENTER_from_every_JT(1:3) + &
        & J_meV * MAGVECTOR_JT(1:3) * sign(1D0,SMT(JT)) * CONC(JT)
!
! Additional sign manipulation to get the energy, with Kamil on 27.5.2020:
!
! H(IT_CENTER) = - \sum_{JT=1}^{NT} J_ij * e_i * e_j, with e_i, e_j versors pointing
! up or down:
!       Heisenberg_energy_from_IT_CENTER_to_every_JT = Heisenberg_energy_from_IT_CENTER_to_every_JT - &
!         & J_meV * sign(1D0,SMT(IT))*sign(1D0,SMT(JT))
!
! In substitutional CPA cases, the orientation is here given position-wise, i.e. indexed 
! by the common IQ, not IT.
!
      MAGVECTOR_IT_CENTER(1) = sin( deg_to_rad(MTET_Q(IQ_CENTER)) ) * cos( deg_to_rad(MPHI_Q(IQ_CENTER)) )
      MAGVECTOR_IT_CENTER(2) = sin( deg_to_rad(MTET_Q(IQ_CENTER)) ) * sin( deg_to_rad(MPHI_Q(IQ_CENTER)) )
      MAGVECTOR_IT_CENTER(3) = cos( deg_to_rad(MTET_Q(IQ_CENTER)) )
!
      Heisenberg_energy_from_IT_CENTER_to_every_JT = Heisenberg_energy_from_IT_CENTER_to_every_JT - &
        & J_meV * &
        & dot_product( MAGVECTOR_IT_CENTER(1:3) * CONC(IT_CENTER) * sign(1D0,SMT(IT_CENTER)) / norm2(MAGVECTOR_IT_CENTER(1:3)), &
        &                     MAGVECTOR_JT(1:3) * CONC(JT) *        sign(1D0,SMT(JT))        / norm2(MAGVECTOR_JT(1:3)) )
!
      write(iounit_xyzfile,'(" ",i0," ",3f15.8,"   ",3f15.8)') &
        & Z_ATOMIC(JT), &
        & ALAT*AU_TO_ANGSTROM*POS_JT(1:3), &
        & MAGVECTOR_JT(1:3) * SMT(JT) / SMT_rescale
!
! maptype = 1:
!
!       if( prepare_UppASD_input ) write(iounit_uppasd,'(i0," ",i0," ",3f15.8, &
!         & " ",f15.8," ",f15.8)') IT_CENTER,JT,POS_JT(1:3),J_meV / RY_to_EV, DR
!
! maptype = 2:
!
! Convertion into (mRy) on-the-fly:
!
      if( prepare_UppASD_input ) write(iounit_uppasd,'(i0," ",i0," ",3i7, &
        & " ",f15.8," ",f15.8)') IT_ONLY_RELEVANT(IT_CENTER),IT_ONLY_RELEVANT(JT), &
        & N1,N2,N3, J_meV / RY_to_EV, DR
!
! Store significantly large Jij's, magnitude, and end-point within
! the cluster around IT_CENTER:
!
      if( abs(J_meV).gt.J_ij_threshold ) then
!
        N_J_ij_large = N_J_ij_large + 1
!
        if( N_J_ij_large.gt.max_N_J_ij_large ) then
!
          write(6,'("max_N_J_ij_large=",i0, &
            & " was too small for N_J_ij_large. Reallocating & re-parsing...")') max_N_J_ij_large
!  
          deallocate( J_ij_target_large_from_IT,J_ij_value_large_from_IT, &
            & J_ij_endpoints )
!
          close(iounit_xyzfile)
          close(iounit_KKR_jfile)
!
          if( prepare_UppASD_input ) close(iounit_uppasd)
!
          max_N_J_ij_large = max_N_J_ij_large * 2
!
          goto 149
!
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
    endif ! IT.eq.IT_CENTER
!
  elseif( parsing_attempt.eq.2 ) then
    if( firstmessage ) then
      firstmessage = .false.
!
      write(6,'("Assuming jfileKKRname: ",a," is in ''jrs.omni'' Hutsepot format")') trim(jfileKKRname)
    endif  
!     1   Na        0.0000000000        3.9674232552        8.6447720000    5   As        0.0000000000        3.9674232552        2.8771300000           0.0000000000        0.0000000000        0.0000000000        0.0000000000    0.00000000000000000000E+00
!
! Safeguard:
    if( (abs(MTET_Q(IT_CENTER)).gt.0.00001D0).or.(abs(MPHI_Q(IT_CENTER)).gt.0.00001D0) ) &
      & stop "not globally colinear ground state not supported (?) yet from Hutsepot. Aborting."
    
    if( NOQ(IT_CENTER).ne.1 ) stop "No support yet for CPA substitutions. Aborting."
!
    read(iounit_KKR_jfile,'(a)',END=1000,ERR=600) textline
    read(textline,*,ERR=600) &
      & IT,TXT_IT,POS_IT(1:3), &
      & JT,TXT_JT,POS_JT(1:3), &
      & DX,DY,DZ, DR, J_Ry
!
! No support yet for CPA substitutional disorder:
!
    JQ = JT
!     
    J_meV = J_Ry * RY_to_EV * 1000D0
!
! Only write end-points from the chosen initial atomic type:
!
    if( IT.EQ.IT_CENTER ) then
!
! If the end point lies too far (but only if it is also not end-point of
! a signicant Jij's, skip it from visualization:
!
      if( (DR.GT.MAX_DR).and.(.not.(abs(J_meV).gt.J_ij_threshold)) ) goto 580
!
      NT_CLUSTER = NT_CLUSTER + 1
!
! ! Different conventions, not needed.      POS_JT(1:3) = QBASQ(1:3,IQ)
!
! End point location:
!
      POS_JT(1) = POS_IT(1) + DX
      POS_JT(2) = POS_IT(2) + DY
      POS_JT(3) = POS_IT(3) + DZ
!
! End point mag.moment ground state orientation:
!
      MAGVECTOR_JT(1) = sin( deg_to_rad(MTET_Q(JQ)) ) * cos( deg_to_rad(MPHI_Q(JQ)) )
      MAGVECTOR_JT(2) = sin( deg_to_rad(MTET_Q(JQ)) ) * sin( deg_to_rad(MPHI_Q(JQ)) )
      MAGVECTOR_JT(3) = cos( deg_to_rad(MTET_Q(JQ)) )
!
      write(iounit_xyzfile,'(" ",i0," ",3f15.8,"   ",3f15.8)') &
        & Z_ATOMIC(JT), &
        & AU_TO_ANGSTROM*POS_JT(1:3), &
        & MAGVECTOR_JT(1:3) * SMT(JT) / SMT_rescale
!
      if( prepare_UppASD_input ) write(iounit_uppasd,'(i0," ",i0," ",3i7, &
        & " ",f15.8," ",f15.8)') IT_ONLY_RELEVANT(IT_CENTER),IT_ONLY_RELEVANT(JT), &
        & NINT(POS_JT(1:3) / ALAT), &
        & J_Ry*1000, SQRT( POS_JT(1)**2 + POS_JT(2)**2 + POS_JT(3)**2 ) / ALAT
!
! Store signicantly large Jij's, magnitude, and end-point within
! the cluster around IT_CENTER:
!
      if( abs(J_meV).gt.J_ij_threshold ) then
!
        N_J_ij_large = N_J_ij_large + 1
!
        if( N_J_ij_large.gt.max_N_J_ij_large ) then
!
          write(6,'("max_N_J_ij_large=",i0, &
            & " was too small for N_J_ij_large. Reallocating & re-parsing...")') max_N_J_ij_large
!  
          deallocate( J_ij_target_large_from_IT,J_ij_value_large_from_IT, &
            & J_ij_endpoints )
!
          close(iounit_xyzfile)
          close(iounit_KKR_jfile)
!
          if( prepare_UppASD_input ) close(iounit_uppasd)
!
          max_N_J_ij_large = max_N_J_ij_large * 2
!
          goto 149
!
        endif
!
! Adding +1, because of the first atom at the origin added by hand above:
!
        J_ij_target_large_from_IT(N_J_ij_large) = NT_CLUSTER + 1
        J_ij_value_large_from_IT(N_J_ij_large) = J_meV
!
        goto 580
!
      endif
!
580   continue
!
    endif
!
!     stop "not implemented yet. Aborting."
!
  elseif( parsing_attempt.eq.3 ) then
    write(6,'("Assuming jfileKKRname: ",a," is in SPRKKR dev.version format")') trim(jfileKKRname)
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
  write(6,'("ERR in jfileKKRname: ",a," on iline=",i0, &
    & " parsing textline: ",a)') trim(jfileKKRname),iline,textline
  stop "file format problems? Aborting."
!
1000 continue
  N_J_ij = iline - 1
  write(6,'("Whole jfileKKRname: ",a," read-in, up to iline=",i0)') trim(jfileKKRname),N_J_ij
!
  write(6,'("Depicting with arrows Jij from TXT_T(IT_CENTER=",i0,")=",a, &
    & " on IQ_CENTER=",i0, &
    & " i.e. cluster atomno=",i0," found larger than J_ij_threshold=", &
    & f15.8," and lying within MAX_DR=",f15.8," (ALAT)")') IT_CENTER,TXT_T(IT_CENTER), &
    & IQ_CENTER, &
    & 1, J_ij_threshold,MAX_DR
!
  if( N_J_ij_large.ge.1 ) then
    write(6,'("J_{i=",i0,",j=",i0,"}=",f15.8)') (1, &
      & J_ij_target_large_from_IT(i_J_ij_large), &
      &  J_ij_value_large_from_IT(i_J_ij_large), i_J_ij_large=1,N_J_ij_large)
  endif
!
  MAX_J_IJ = maxval( abs(J_ij_value_large_from_IT) )
  write(6,'("Strongest exchange coupling shown from IT_CENTER=",i0," (in meV, absolute value):", &
    & f15.8)') IT_CENTER,MAX_J_IJ
!
!   allocate( Jij_list(6,N_J_ij) )
!   rewind(iounit_KKR_jfile)
!   goto 149
! 1100 continue
!
! Save & close output files:
!
  close(iounit_KKR_jfile)
!
  if( prepare_UppASD_input ) close(iounit_uppasd)
!
end subroutine initparse_JXC_SP_SREL

pure real(kind(0.d0)) function deg_to_rad(angle)
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
!
    ich = ichar(string(i:i)) ! Get position of string(i:i) in ASCII set
!
! Range A to Z:
    if( (ich.ge.65).and.(ich.le.90) ) then
!
! Range a to z:
      string_to_lowercase(i:i) = char(ich+32)
!
    else
!
      string_to_lowercase(i:i) = string(i:i)
!
    endif
!
  enddo
!
end function string_to_lowercase
!
subroutine create_T_scan_scripts()
  implicit none
  integer,parameter :: iounit_T_scan=665
!
!
! List of trial temperatures in Kelvin, one by line, to be edited by hand:
  open(iounit_T_scan,file='sampled_temperatures.dat',action='WRITE')
!
  write(iounit_T_scan,'("1")')
  write(iounit_T_scan,'("10")')
  write(iounit_T_scan,'("20")')
  write(iounit_T_scan,'("30")')
!
  close(iounit_T_scan)
!
! Master submission script, exploring the range of trial temperatures in parallal:
  open(iounit_T_scan,file='viper_submit_range.sh',action='WRITE')
!
  write(iounit_T_scan,'("#!/bin/bash")')
  write(iounit_T_scan,'("for Temp in `cat sampled_temperatures.dat`; do \")')
  write(iounit_T_scan,'(" echo ""Submitting for temp=""$Temp ; qsub viper_UppASD-commandline_temp.job $Temp ; \")')
  write(iounit_T_scan,'("done")')
!
  close(iounit_T_scan)
!
! Submission script of each trial temperature for the queuing system:
  open(iounit_T_scan,file='viper_UppASD-commandline_temp.job',action='WRITE')
!
  write(iounit_T_scan,'("#!/bin/bash")')
  write(iounit_T_scan,'("#$ -cwd")')
  write(iounit_T_scan,'("#$ -N hello_world")')
!
  write(iounit_T_scan,'("#$ -q sprkkr.q")')
! Viper entirely single-node queue. Risk of overlap with split-node queue, hence commented out.
  write(iounit_T_scan,'("######$ -pe slots 20")')
!
! Number of OpenMP threads for each execution:
  write(iounit_T_scan,'("#$ -pe mpi 40")')
  write(iounit_T_scan,'("#$ -V")')
  write(iounit_T_scan,'("echo ''Running with NSLOTS=''$NSLOTS")')
!
!   write(iounit_T_scan,'("#! /bin/bash")')
!
! Trial temperature as first invocation argument:
  write(iounit_T_scan,'("Temp=$1")')

  write(iounit_T_scan,'("echo ''Simulating Temp=''$Temp ''...''")')
  write(iounit_T_scan,'("echo '' '' ")')
  write(iounit_T_scan,'("mkdir Temp_$Temp/ 2>/dev/null")')

! Enter each prepared directory, and customize its trial temperature from the list:
  write(iounit_T_scan,'("cp inpsd.dat jfile posfile momfile Temp_$Temp/")')
  write(iounit_T_scan,'("cd Temp_$Temp/")')
  write(iounit_T_scan,'("sed -i ""s/TEMP/$Temp/g"" inpsd.dat")')

! Execute:
  write(iounit_T_scan,'("TMPDIR=/tmp OMP_NUM_THREADS=$NSLOTS ~/UppASD/source/sd > output.log")')

  close(iounit_T_scan)
end subroutine create_T_scan_scripts
!
subroutine create_gnuplot_T_critical_fit(largest_SMT)
  implicit none
  integer,parameter :: iounit_Tc_fit=663
  real(kind(0.d0)),intent(IN) :: largest_SMT
!
  open(iounit_Tc_fit,file='Tc_fit.gnuplot',action='WRITE')
!
! Parse the Monte Carlo outcome, type-resolved (to account for compensated cases like AF):
  write(iounit_Tc_fit,'("# Parse results from UppASD output in each directory i.e. sampled temperature:")')
  write(iounit_Tc_fit,'("! . printM-type_resolved.sh")')

  write(iounit_Tc_fit,'("set xlabel ''Temperature (K)'' ; set grid ")')
!
  write(iounit_Tc_fit,'("set ytics nomirror ; set y2tics nomirror")')
  write(iounit_Tc_fit,'("set ylabel ''Type-resolved mag.moment'' tc rgb ''black''")')
  write(iounit_Tc_fit,'("set y2label ''Type-resolved susceptibility'' tc rgb ''red''")')
! 
! set ylabel 'spin mag.moment (\mu_B)'
! set title 'CuMnAs AF from SPRKKR SP-SREL'
! 
  write(iounit_Tc_fit,'("m(x) = norm * (1 - x/Tc)**beta")')
! 
  write(iounit_Tc_fit,'("# Initial Ansatz of critical temperature fit parameters (adjust as needed):")')
!
  write(iounit_Tc_fit,'("Tc = 500 ; beta = .2 ; norm=",f15.8)') abs(largest_SMT)
! 
!   write(iounit_Tc_fit,'("fit m(x) ''thermal.dat'' u 1:2 via Tc,beta,norm")')
! Always present first element:
  write(iounit_Tc_fit,'("  fit m(x) ''thermal_type_1.dat'' u 1:3 via Tc,beta")') 
  ! Skipping fit of the maximal value, which may be best left to the DFT result.
!   write(iounit_Tc_fit,'("  fit m(x) ''thermal_type_1.dat'' u 1:3 via Tc,beta,norm")')
!
  write(iounit_Tc_fit,'("print ''Fit results for type-resolved magnetization dependence ")',advance='NO')
  write(iounit_Tc_fit,'(" m(x) = norm * (1 - x/Tc)**beta : Tc='',Tc, '' beta='',beta, '' norm='',norm")')
!
! Always present first element:
  write(iounit_Tc_fit,'("p ''thermal_type_1.dat'' u 1:3  w p pt 5 ps 2 lw 4 lc rgb ''black'', ")',advance='NO')
  write(iounit_Tc_fit,'(" m(x) w l lw 3 lc rgb ''black'', \")')
  write(iounit_Tc_fit,'(" '''' u 1:5 axes x1y2 w lp lw 3 lc rgb ''red''")')
!
  close(iounit_Tc_fit)
!
end subroutine create_gnuplot_T_critical_fit
!
subroutine parse_IT_resolved_MC_averages(NT_RELEVANT)
! Auxiliary script to extract input for create_gnuplot_T_critical_fit
  implicit none
  integer,intent(IN) :: NT_RELEVANT
  integer,parameter :: iounit_parsescript=664
!
  open(iounit_parsescript,file='printM-type_resolved.sh',action='WRITE')
!
  write(iounit_parsescript,'("#! /bin/bash")')
  write(iounit_parsescript,'("echo ''Parsing results from trial temperatures listed in file: sampled_temperatures.dat''")')
!
  write(iounit_parsescript,'("m_max=`tail -n1 Temp_*/cumulants.*.out | grep E | awk ''{ print $2 }'' | sort -gr | head -1 `")')
! #cv_max=`tail -n1 Temp_*/cumulants.*.out | grep E | awk '{ print $7 }' | sort -gr | head -1 `
  write(iounit_parsescript,'("cv_max=1")')
  
  write(iounit_parsescript,'("x_max=`tail -n1 T*/cumulants.*.out | grep E | awk ''{ print $6 }'' | sort -gr | head -1 `")')
! Headers:
  write(iounit_parsescript,'("echo ''# Temp.   Mavg     UBinder    Susc.      Cv'' > thermal.dat")')
  write(iounit_parsescript,'("echo ''# Temp.   Mavg     UBinder    Susc.      Cv'' > thermal.norm.dat")')
!
  write(iounit_parsescript,'("for Temp in `cat sampled_temperatures.dat`")')
  write(iounit_parsescript,'("   do")')
  write(iounit_parsescript,'("     echo ''Parsing AVERAGE from calculation performed at Temp=''$Temp")')
!
  write(iounit_parsescript,'("     tail -1 Temp_$Temp/cumulants.*.out | ")',advance='NO')
  write(iounit_parsescript,'(" awk ''{ printf ""%6i  %f  %f  %f  %f\n"", temp, $2,$5,$6,$7 }'' temp=$Temp >> thermal.dat")')
!
  write(iounit_parsescript,'("     tail -1 Temp_$Temp/cumulants.*.out | ")',advance='NO')
  write(iounit_parsescript,'(" awk ''{ printf ""%6i  %f  %f  %f  %f\n"", temp, $2/m_max,$5,$6/x_max,$7/cv_max }'' ")',advance='NO')
  write(iounit_parsescript,'(" m_max=$m_max x_max=$x_max cv_max=$cv_max temp=$Temp >> thermal.norm.dat")')
!
  write(iounit_parsescript,'("   done")')

  write(iounit_parsescript,'("# Antiferromagnetic case, with 2 magnetic sublatties:")')
!
  write(iounit_parsescript,'("NTYPES=",i0)') NT_RELEVANT
!
  write(iounit_parsescript,'("echo ''Parsing TYPE RESOLVED...''")')
!
  write(iounit_parsescript,'("for type in `seq 1 $NTYPES`")')
  write(iounit_parsescript,'("   do ")')
  write(iounit_parsescript,'("     echo ''# Temp. Type  Mavg     UBinder    Susc.'' > ''thermal_type_''$type''.dat''")')
  write(iounit_parsescript,'("     for Temp in `cat sampled_temperatures.dat`")')
  write(iounit_parsescript,'("        do")')
  write(iounit_parsescript,'("          echo ''Parsing Type=''$type '' from calculation performed at Temp=''$Temp")')
!   
  write(iounit_parsescript,'("          tail -$NTYPES Temp_$Temp/projcumulants.*.out | ")',advance='NO')
  write(iounit_parsescript,'(" awk -v type=$type ''($2 == type) ")',advance='NO')
  write(iounit_parsescript,'(" { printf ""%6i %2i %f  %f  %f\n"", temp, type, $3,$6,$7 }'' temp=$Temp ")',advance='NO')
  write(iounit_parsescript,'("  >> ''thermal_type_''$type''.dat''")')
!   
  write(iounit_parsescript,'("        done")')
  write(iounit_parsescript,'("  done")')
!
  close(iounit_parsescript)
!
end subroutine parse_IT_resolved_MC_averages
!
subroutine write_cell(chosen_xyz_cell_format,xyzfilename,ALAT,ABAS,NQ,NT,MTET_Q,MPHI_Q,SMT,SMT_rescale,Z_ATOMIC, &
    & ITOQ,NOQ,CONC,QBASQ)
  implicit none
!
  real(kind(0.d0)),parameter :: AU_TO_ANGSTROM = 0.52917721067D0
  integer,parameter :: iounit_xyzfile=656, &
! Hard-coded max. number of CPA alternative on any sublattice IQ:
    & nalts_max=10
!
  character(len=*),intent(IN) :: xyzfilename
  integer,intent(IN) :: NQ,NT,Z_ATOMIC(NT),ITOQ(nalts_max,NQ),NOQ(NQ)
  real(kind(0.d0)),intent(IN) :: SMT_rescale, SMT(NT), &
    & QBASQ(3,NQ),MTET_Q(NQ),MPHI_Q(NQ),ALAT,CONC(NT), ABAS(3,3)
  logical,intent(IN) :: chosen_xyz_cell_format
!
  real(kind(0.d0)) MAGVECTOR(3),deg_to_rad
  integer IQ,IT,IA
!
  if( chosen_xyz_cell_format ) then
!
    write(6,'("Storing 0D-periodic unit cell in XYZ format as xyzfilename: ",a)') trim(xyzfilename)
!
  else
!
    write(6,'("Storing 3D-periodic unit cell in XSF format as xyzfilename: ",a)') trim(xyzfilename)
    write(6,'(" with 3D lattice basis vectors: ABAS(1:3,1)=",3f15.8)') ABAS(1:3,1)
    write(6,'("                                ABAS(1:3,2)=",3f15.8)') ABAS(1:3,2)
    write(6,'("                                ABAS(1:3,3)=",3f15.8)') ABAS(1:3,3)
!
  endif
!
! Safeguard:
!
  if( SMT_rescale.le.0D0 ) &
    stop "SMT_rescale cannot be 0 or negative"
!
  if( ALAT.le.0D0 ) &
    stop "ALAT cannot be 0 or negative"    
!
  open(iounit_xyzfile,file=trim(xyzfilename),action='WRITE')
!
  if( chosen_xyz_cell_format ) then
    write(iounit_xyzfile,'(i0)') NQ
!
! XCrysDen defaults for forces is too long for depicting SMT in (mu_B): rescale by a fixed factor:
!
    write(iounit_xyzfile,'("Dump from SPRKKR potential file, with SMT_rescale=",f15.8)') SMT_rescale
!
  else
!
    write(iounit_xyzfile,'("  CRYSTAL")')
    write(iounit_xyzfile,'("PRIMVEC")')
!
    write(iounit_xyzfile,'(3f15.8)') ABAS(1:3,1)*ALAT*AU_TO_ANGSTROM
    write(iounit_xyzfile,'(3f15.8)') ABAS(1:3,2)*ALAT*AU_TO_ANGSTROM
    write(iounit_xyzfile,'(3f15.8)') ABAS(1:3,3)*ALAT*AU_TO_ANGSTROM
!
    write(iounit_xyzfile,'("CONVVEC")')
    write(iounit_xyzfile,'(3f15.8)') ABAS(1:3,1)*ALAT*AU_TO_ANGSTROM
    write(iounit_xyzfile,'(3f15.8)') ABAS(1:3,2)*ALAT*AU_TO_ANGSTROM
    write(iounit_xyzfile,'(3f15.8)') ABAS(1:3,3)*ALAT*AU_TO_ANGSTROM
!
    write(iounit_xyzfile,'("PRIMCOORD")')
    write(iounit_xyzfile,'(i0," 1")') NQ
!
  endif
!
  do IQ=1,NQ
    do IA=1,NOQ(IQ)
!
      IT = ITOQ(IA,IQ)
!
      MAGVECTOR(1) = SMT(IT) * sin( deg_to_rad(MTET_Q(IQ)) ) * cos( deg_to_rad(MPHI_Q(IQ)) )
      MAGVECTOR(2) = SMT(IT) * sin( deg_to_rad(MTET_Q(IQ)) ) * sin( deg_to_rad(MPHI_Q(IQ)) )
      MAGVECTOR(3) = SMT(IT) * cos( deg_to_rad(MTET_Q(IQ)) )
!
! XCrysDen defaults for forces is too long for depicting SMT in (mu_B): rescale by a fixed factor:
!
      MAGVECTOR(1:3) = MAGVECTOR(1:3) / SMT_rescale
!
! XYZ format does not specify anything for CPA substitutions. 
!
! XCrysDen will just ignore further columns: these are used for the purpose of CPA alternatives.
      if( IA.eq.1 ) then
        write(iounit_xyzfile,'(" ",i0," ",3f15.8," ",3f15.8)',advance='NO') Z_ATOMIC(IT), &
          & ALAT*AU_TO_ANGSTROM*QBASQ(1:3,IQ), MAGVECTOR(1:3)
      else
        write(iounit_xyzfile,'(" Z(ITOQ(IA=",i0,",IQ=",i0,"))=",i0," with CONC=",f15.8,", ")',advance='NO') &
          & IA,IQ, Z_ATOMIC(IT), CONC(IT)
      endif
!
    enddo
    write(iounit_xyzfile,'("")')
  enddo
!
  close(iounit_xyzfile)
!
end subroutine write_cell

subroutine get_critical_temperature_mean_field(jfileKKRname,NT,NQ,NOQ,ITOQ,CONC)
!
! From xcpljij.f, for cross-checking the parsing etc.:
!
  implicit none
  integer,parameter :: nalts_max = 10, &
    & iounit_KKR_jfile=655
! Boltzmann constant:
  real(kind(0.d0)),parameter :: KB_SI = 1.38064852D-23, &
    & EV_J = 1.6021766208D-19
!
  character(300), intent(IN) :: jfileKKRname
  integer,intent(IN) :: NT,NQ,NOQ(NQ),ITOQ(nalts_max,NQ)
  real(kind(0.d0)),intent(IN) :: CONC(NT)
!
  real(kind(0.d0)),allocatable :: W4(:)
  integer LWMAX,IQ,IA,IT, JQ,IO,JT, LWK, BLAS_INFO,iline, N1,N2,N3, N_J_ij
!   read_info,
  real(kind(0.d0)) T_C(NT,NT), &
    & DX,DY,DZ, DR, J_meV,J_eV
!     POS_IT(3),POS_JT(3), 
!   character(7) TXT_T(NT)
  character(300) textline,headerline
!
  real(kind(0.d0)) TCURIE,WI(NQ),WR(NQ),WW1(NQ,NQ),WW2(NQ,NQ)
!
  write(6,'("")')
  write(6,'("Re-opening jfileKKRname: ",a," to estimate mean-field T_c...")') trim(jfileKKRname)
  open(iounit_KKR_jfile,file=trim(jfileKKRname),action='READ',ERR=20)
  goto 30
!
20 stop "Error opening jfileKKRname. Aborting."
30 continue
!
  T_C(1:NT,1:NT) = 0.D0
!
! Scan for the payload. Only SPRKKR 8.5.0 IREL=2 version supported for now.
!
  headerline = "        IT   IQ   JT   JQ   N1 N2 N3    DRX    DRY    DRZ    DR        J [meV]         J [eV]"
40 continue
!
  iline = 0
!
50 continue
!
  iline = iline + 1
!
  read(iounit_KKR_jfile,'(a)',ERR=200,END=200) textline
  if( index(textline,trim(headerline)).gt.0 ) goto 500
  goto 50
!
! Safeguard:
200 continue
  write(6,'("ERR/EOF in jfileKKRname: ",a," at iline=",i0, &
    & " at/after last read: textline=",a)') trim(jfileKKRname),iline,textline
  write(6,'(" before finding header:",a)') headerline
!
500 continue
  iline = iline + 1
!
! Parse as SPRKKR 8.5.0 format:
!
  read(iounit_KKR_jfile,'(a)',END=1000,ERR=600) textline
  read(textline,*,ERR=600) &
    & IT,IQ, JT,JQ, N1,N2,N3, &
    & DX,DY,DZ, DR, J_meV,J_eV
! Debug:
!   write(6,*) textline  

! In SPRKKR, the sign of Heisenberg exhange energy 
! is always referring to a FM coupling. Swap the sign
! if the two atoms are instead in a AF ground state:
!
!     J_meV = sign(1D0,SMT(IT))*sign(1D0,SMT(JT)) * J_meV
!     J_eV = sign(1D0,SMT(IT))*sign(1D0,SMT(JT)) * J_eV
!
  T_C(IT,JT) = T_C(IT,JT) + CONC(IT)*CONC(JT)*J_eV
!
  goto 500
!
600 continue
  write(6,'("ERR in jfileKKRname: ",a," on iline=",i0, &
    & " parsing textline: ",a)') trim(jfileKKRname),iline,textline
  stop "file format problems? Aborting."
!
1000 continue
  N_J_ij = iline - 1
  write(6,'("Whole jfileKKRname: ",a," read-in, up to iline=",i0)') trim(jfileKKRname),N_J_ij
  
!   print *,T_C(1:NT,1:NT)
!
! Continue the Jij's scan:
!
!
  LWMAX = 20*NT*NT
  allocate( W4(LWMAX) )
!
! Safe blanking out, although they should be written-in anyways:
!
!   T_C(1:NT,1:NT) = 0.D0
!
! Format of the file: SPRKKR 8.5.0, not Hutsepot:
!
! Skip the header:
! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! 
! load into memory the Jij's:
!
!   read(textline,*,ERR=200,IOSTAT=read_info) &
!     & IT,TXT_T(IT),POS_IT(1:3), &
!     & JT,TXT_T(JT),POS_JT(1:3), &
!     & DX,DY,DZ, DR, J_Ry
! !
!   if( read_info.lt.0 ) goto 400  
! !
! 200 stop "Error while re-parsing JXC_J_ij.dat: aborting."
! !
! 400 continue
  loop_IQ: DO IQ = 1,NQ
!
    loop_IA: DO IA = 1,NOQ(IQ)
!
      IT = ITOQ(IA,IQ)
!
      loop_JQ: DO JQ = 1,NQ
!
!         JXC0(IT) = 0D0
!
        loop_IO: DO IO = 1,NOQ(JQ)
!
          JT = ITOQ(IO,JQ)
!
! meV conversion:
!
! !           JXC0(IT) = CONC(JT)*JXCITJT_0(IT,JT)*1.D-3 
!           T_C(IT,JT) = CONC(IT)*JXC0(IT)
!
        ENDDO loop_IO
!
      ENDDO loop_JQ
!
    ENDDO loop_IA
!
  ENDDO loop_IQ
!
!   NTC = NT
!
  CALL DGEEV('N','N',NT,T_C,NT,WR,WI,WW1,NT,WW2,NT,W4,-1,BLAS_INFO)

  LWK = MIN(LWMAX,INT(W4(1)))
  CALL DGEEV('N','N',NT,T_C,NT,WR,WI,WW1,NT,WW2,NT,W4,LWK,BLAS_INFO)
!
  IF( BLAS_INFO.EQ.0 ) THEN
!
    TCURIE = 0D0
!
    DO IT = 1,NT
!
      TCURIE = MAX(TCURIE,T_C(IT,IT))
!
    ENDDO
!
    write(6,'("Mean field critical temperature estimate: ",f15.8)') TCURIE*(2D0/3D0)*EV_J/KB_SI
    write(6,'("")')
!
  ELSE
!
    write(6,'("BLAS_INFO=",i0)') BLAS_INFO
    stop "WARNING:  <SGEGV> INFO <> 0"
!
  ENDIF
!
  close(iounit_KKR_jfile)
  deallocate( W4 )
!
end subroutine get_critical_temperature_mean_field
