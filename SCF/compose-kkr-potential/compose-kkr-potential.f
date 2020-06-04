      PROGRAM compose_kkr_potential
c     =============================

c Program composes a SPR-KKR potential in the Munchen format by modifying
c   an existing Munchen potential file.
c   Modifying potentials can be either Munchen or Julich format potentials.


      implicit none

      integer ntmax, irdm
      parameter ( ntmax=1001, irdm=1001 )

      integer iargc, IT, NMODPOT, I, NARG, N1(99), IEND, IBEG, IWRI, 
     $     JP_MODIF, JP_ORIG, JRWSMX, NR(NTMAX), IR, POTFMTINP,
     $     NT, JRWS(NTMAX), ITDAT
      double precision RAD(IRDM), VT(IRDM), BT(IRDM), EF, EXPDX, 
     $     RADMOD(IRDM,NTMAX), VTMOD(IRDM,NTMAX), BTMOD(IRDM,NTMAX),
     $     EF_MODIF(NTMAX), R1(NTMAX), EF_ORIG, RWS_MODIF, DXCALC,
     $     RWS_ORIG,  RWS(NTMAX), ef_expl, oldfr(ntmax), newfr(ntmax),
     $     sum
      character TEXT*256, INFILE*128, BASFILE*128, NEWFILE*128,
     $     FMT08*40
      logical MODIF(NTMAX)


      if ( iargc() .ne. 2 .and. iargc() .ne. 3 ) then
        write (*,'(a,/,9(a))') 'Two or three arguments needed: ',
     $       '   <Old-KKR-pot-file> (in)    <New-KKR-pot-file> (out)',
     $       '    [E_Fermi]'
        stop 'Program ended BADLY !'
      endif  

      open(55,status='scratch' )

      do it = 1, ntmax
        modif(it) = .false.

        RWS(it) = 0.0d0                     ! Non-sense default value

        R1(it) = 1.0d-6                     ! Hard-coded defaults
        jrws(it) = 721
      enddo 


c Browses modifying potential files one by one and extracts the potentials

      nmodpot = 0

 10   continue 
        do i = 1, len(text)
          text(i:i) = ' '
        enddo 
        write (6,'(19(a))') '<Modif_Pot_File>    ',
     $       '<Ipot-modif (from)>  <RWS-modif (from)>    ',
     $       '<Ipot-orig (to)>  <RWS-orig (to)>    (or EOF) !'
        read (5,'(a)',end=90) text

        if ( text(1:1) .eq. '#'  .or.  text(1:1) .eq. '!' ) goto 10

        call countarg2( text, narg, n1 )
        call findsize( text, ibeg, iend )

        if ( narg .ne. 5 ) then
          write (*,*) 'Illegal format of the input line:'
          write (*,*) 'Input =''', text(ibeg:iend), ''''
          stop 'Program ended BADLY !'
        endif 

        infile = text(n1(1):n1(2)-2)
        open( 25, file=infile, status='old', action='read')

        write (55,'(a)') text(n1(2)-1:iend+1)
        backspace(55)
        read (55,*) JP_MODIF, rws_modif, JP_ORIG, rws_orig

        if ( JP_MODIF .le. 0 .or. JP_MODIF .gt. 2001 .or.
     $       JP_ORIG .le. 0 .or. JP_ORIG .gt. 2001 ) then
          write (*,*) 'Illegal JP_MODIF, JP_orig !'
          write (*,*) 'Input line: ''', text(n1(2)-1:iend+1), ''''
          stop 'Program ended BADLY !'
        endif 
        if ( RWS_MODIF .le. 0.1 .or. RWS_MODIF .gt. 10.0 .or.
     $       RWS_ORIG .le. 0.1 .or. RWS_ORIG .gt. 10.0 ) then
          write (*,*) 'Illegal RWS_MODIF, RWS_orig !'
          write (*,*) 'Input line: ''', text(n1(2)-1:iend+1), ''''
          stop 'Program ended BADLY !'
        endif 

c Generates radial grid for potential type JP_ORIG

        RWS(JP_ORIG) = RWS_ORIG
        

c Potential from the modifying potential file is read for JRWS mesh points
        
        call readpot( 25, JP_MODIF, JRWSMX, RWS_MODIF, RAD, VT, BT, EF )

        close(25)


c Potential to be later inserted is stored together with its radial grid

        NR(JP_ORIG) = JRWSMX

        do IR = 1, NR(JP_ORIG)
          RADMOD(IR,JP_ORIG) = RAD(IR)
          VTMOD(IR,JP_ORIG) = VT(IR)
          BTMOD(IR,JP_ORIG) = BT(IR)
        enddo  

        EF_MODIF(JP_ORIG) = EF

        modif(jp_orig) = .true.

        nmodpot = nmodpot + 1
      goto 10

 90   write (*,'(i4,a)') nmodpot,
     $     ' modifying potentials identified and read'


c Opens the original potential file to be modified

      call getarg( 1, basfile )
      open (15, file=basfile, status='old',action='read')

      call getarg( 2, newfile )
      open (16, file=newfile, status='unknown',action='write')


      iwri = 0

c Sets POTFMTINP to 7 by default

      CALL READKWINT(15,'FORMAT    ',POTFMTINP,7,iwri,0)      


      CALL READKWINT(15,'NT        ',NT,0,IWRI,0)
      if ( NT .eq. 0 ) then
        CALL READKWINT(15,'NTCLU     ',NT,0,IWRI,1)
      endif 
      if ( nt.eq.0 ) stop 'original potential file nt.eq.0 !'


c Reads and copies all lines BEFORE the line with Fermi energy

      call copy_before_fermi( 15, 16 )

c Reads the Fermi energy of the original potental file

      CALL READKWREAL(15,'EF        ',EF_ORIG,0.0D0,IWRI,1)

      if ( iargc().eq.3 ) then
        do i = 1, len(text)
          text(i:i) = ' '
        enddo 
        call getarg( 3, text )
        write (55,'(9(a))') '  ', text, '  '
        backspace(55)
        read (55,*) EF_EXPL
      else 
        EF_EXPL = EF_ORIG
      endif 

      FMT08 = '(A10,3(1P,E22.14))'
      WRITE (16,FMT=FMT08) 'EF        ', EF_EXPL

c Reads and copies further lines BEFORE the part containing the potentials

      call copy_between_fermi_and_potential( 15, 16 )


c Reads potentials and substitutes them with interpolated modifying potentials 
c   if needed 

      do IT = 1, NT
        call readnextpot( 15, POTFMTINP, ITDAT, JRWS, VT, BT )
        
        if ( .not.modif(itdat) ) then

c Unless it is a bulk vacuum sphere, adjust the original potential 
c    to the new Fermi energy

          sum = 0.0d0
          do IR = 1, JRWS(ITDAT)
            sum = sum + abs(VT(IR))
          enddo  
          
          if ( sum .gt. 0.0001 ) then           
            do IR = 1, JRWS(ITDAT)
              VT(IR) = VT(IR) + EF_EXPL - EF_ORIG
            enddo  
          endif 

          call writepot( 16, ITDAT, POTFMTINP, JRWS(ITDAT), VT, BT )
          
        else 

c Generates radial grid as in the original potential file
          
          RAD(1) = R1(ITDAT)
          DXCALC = log(RWS(ITDAT)/R1(ITDAT)) / dble(JRWS(JP_ORIG)-1)
          EXPDX = dexp(DXCALC)
          
          DO IR = 2, JRWS(ITDAT)
            RAD(IR) = RAD(IR-1) * EXPDX
          END DO
          

c Interpolates modifying potential into the radial grid of the original file
c    For numerical reasons, interpolation of potential is done after R multiplication

          do ir = 1, NR(ITDAT)
            oldfr(ir) = RADMOD(ir,ITDAT) * VTMOD(ir,ITDAT)
          enddo 

          call qint2( oldfr, RADMOD(1,ITDAT), NR(ITDAT),
     $                newfr, RAD, JRWS(ITDAT) )

          do ir = 1, JRWS(ITDAT)
            vt(ir) = newfr(ir) / RAD(IR)
          enddo 


          call qint2( BTMOD(1,ITDAT), RADMOD(1,ITDAT), NR(ITDAT),
     $                BT, RAD, JRWS(ITDAT) )

c Unless it is a bulk vacuum sphere, adjust the modifying potential 
c    to the Fermi energy of the original potential

          sum = 0.0d0
          do IR = 1, JRWS(ITDAT)
            sum = sum + abs(VT(IR))
          enddo  
          
          if ( sum .gt. 0.0001 ) then
            do IR = 1, JRWS(ITDAT)
              VT(IR) = VT(IR) + EF_EXPL - EF_MODIF(ITDAT)
            enddo  
          endif 
          

c Writes the modifying potential in place of the original potential

          call writepot( 16, ITDAT, POTFMTINP, JRWS(ITDAT), VT, BT )

        endif 

      enddo 


c Reads and copies all lines AFTER the potential part

      call copy_after_potential( 15, 16 )


      stop
      end




      subroutine READPOT( IOF, JP, NR, RWS, RAD, VT, BT, EF ) 
c     ===================

      implicit none

      integer IOF, JP, NR
      double precision RAD(*), VT(*), BT(*), EF, RWS

      logical form_julich, form_munchen


      form_julich = .false.
      form_munchen = .false.

c Identifies whether the potential is in Munchen or Julich format

      call findformat( iof, form_julich, form_munchen )

      rewind( IOF )

      if ( FORM_JULICH ) then

        call read_julich_pot( IOF, JP, NR, RAD, VT, BT, EF )

      else if ( FORM_MUNCHEN ) then

        call read_munchen_pot( IOF, JP, NR, RWS, RAD, VT, BT, EF ) 

      else 
        write (*,*)
     $       'Neither Julich nor Munchen potential file identified !'
        stop 'Program ended BADLY !'
      endif 

      return 
      end




      SUBROUTINE READ_MUNCHEN_POT( IOF, JP, NR, RWS, RAD, VT, BT, EF ) 
c     ===========================

c Reads potential type JP in Munchen format

      implicit none

      integer IOF, JP, NR
      double precision RAD(*), VT(*), BT(*), EF, RWS


      integer ntmax
      parameter ( ntmax=1001 )

      integer iend, iwri, nt, ir, i, itype, POTFMTINP
      double precision R1, DXCALC, expdx
      character FMT07*16, text*128


      iwri = 0

c Generates the radial grid (hard-coded, determined by RWS)

      NR = 721
      R1 = 1.0d-6
      RAD(1) = R1
      DXCALC = log(RWS/R1) / dble(NR-1)
      EXPDX = dexp(DXCALC)
      
      DO IR = 2, NR
        RAD(IR) = RAD(IR-1) * EXPDX
      END DO




      CALL READKWINT(iof,'FORMAT    ',POTFMTINP,7,IWRI,0)

c Reads Fermi energy 

      CALL READKWREAL(IOF,'EF        ',EF,0.0D0,IWRI,1)

c --------------------------------------------------------------------
c Reads the "occupation" section containing details about radial grids

      CALL READKWINT(IOF,'NT        ',NT,0,IWRI,0)
      if ( NT .eq. 0 ) then
        CALL READKWINT(IOF,'NTCLU     ',NT,0,IWRI,1)
      endif 
      if ( nt.eq.0 ) stop 'READ_MUNCHEN_POT:   nt.eq.0 !'

      if ( nt .gt. ntmax ) then
        stop ' nt .gt. ntmax  !'
      endif 


c -------------------------------
c Reads the potential for type JP

      rewind(iof)

      IF ( POTFMTINP.EQ.6 ) THEN
        FMT07 = '(1P,5E16.9)'
      else 
        FMT07 = '(1P,5E22.14)'
      endif 


 10   continue 
        do i = 1, len(text)
          text(i:i) = ' '
        enddo 

        read (iof,'(a)',end=661) text

        if ( text(1:9).eq.'POTENTIAL' ) goto 190

      goto 10 
 190  continue 

 20   continue 
        do i = 1, len(text)
          text(i:i) = ' '
        enddo 

        read (iof,'(a)',end=662) text

        if ( text(1:6).eq.'TYPE  ' ) then
          call findend ( text, iend)
          write (55,'(9(a))') '  ', text(7:iend), '  '
          backspace(55)
          read (55,*,end=663,err=663) itype
          if ( itype .eq. JP ) goto 290
        endif 

      goto 20 
 290  continue 

      READ (iof,FMT=FMT07) (VT(I),I=1,NR)
      READ (iof,FMT=FMT07) (BT(I),I=1,NR)

      return 

 661  stop 'READ_MUNCHEN_POT:  EOF while scanning POTENTIAL target !'
 662  stop 'READ_MUNCHEN_POT:  EOF while scanning TYPE target !'
 663  stop 'READ_MUNCHEN_POT:  EOF or ERR while reading ITYPE !'

      end




      subroutine read_julich_pot( IOF, ITYP, JWS, RAD, VT, BT, EFNEW )
c     ==========================

c Reads the Julich code ASA potential and transforms it to Munchen form

      implicit none

      integer IOF, ITYP, JWS
      double precision vt(*), bt(*), EFNEW, RAD(*)


      integer irdm
      parameter ( irdm=401 )

      integer IRWS, NCORE, INEW, IH, IS, IR, ICORE, NSPINS, jjr
      double precision RMT, ALAT, Z, RWS, A, B,  RMTNEW, VBC,
     $     VM2Z(irdm,2), vpdn(irdm), vpup(irdm), dummy
      character*128 TITLE, CORELINE
      

      nspins = 2

c
c---> read title of potential card
c
      IH = 0

 10   IH = IH + 1
        if ( IH .gt. ITYP ) goto 90

        do IS = 1, NSPINS
          READ (IOF, '(a)', end=61 ) TITLE

c---  >read muffin-tin radius , lattice constant and new muffin radius
c      (new mt radius is adapted to the given radial mesh)

          READ (IOF,FMT=*) RMT, ALAT, RMTNEW

c---> read nuclear charge , lmax of the core states ,
c     wigner seitz radius , fermi energy and energy difference
c     between electrostatic zero and muffin tin zero

          READ (IOF,FMT=*) Z
          READ (IOF,FMT=*) RWS, EFNEW, VBC

c---> read : number of radial mesh points

          READ (IOF,FMT=*) IRWS
          READ (IOF,FMT=*) A, B
          READ (IOF,FMT=*) NCORE, INEW
          
c
c---> read the different core states : l and energy
c
          IF ( NCORE .ge. 1 ) THEN
            DO ICORE = 1, NCORE
              READ (IOF,'(a)') CORELINE
            END DO
          END IF

c--->  read the input potential without the nuclear pot.
c      Omit first point because it corresponds to zero coordinate     

          jws = irws - 1
          
          READ (IOF,FMT=*) dummy, ( VM2Z(jjr,IS), jjr = 1, jws )

        enddo
      goto 10

 90   continue 


      DO IR = 2, IRWS
        JJR = IR - 1
        RAD(JJR) = B * ( EXP( A * DBLE(IR-1) ) - 1.0D0 )
      ENDDO 


      if ( nspins .eq. 2 ) then
        do jjr = 1, jws
          vpdn(jjr) = VM2Z(jjr,1)
          vpup(jjr) = VM2Z(jjr,2)
        enddo 
      else 
        do jjr = 1, jws
          vpdn(jjr) = VM2Z(jjr,1)
          vpup(jjr) = VM2Z(jjr,1)
        enddo 
      endif  


c Potentials defitions:
c    V_DN = VT - BT + 2Z/R           is=1
c    V_UP = VT + BT + 2Z/R           is=2

c Isolates spin-independent and spin-dependent part and adds Coulombic part

      do jjr = 1, jws
        bt(jjr) = 0.50d0 * ( vpup(jjr) - vpdn(jjr) )
      enddo 

      do jjr = 1, jws
        vt(jjr) = 0.50d0 * ( vpup(jjr) + vpdn(jjr) ) - 2.0d0*Z/rad(jjr)
      enddo 

      return 

 61   stop 'EOF before reaching appropriate Julich potential section'

Cabc 9030 FORMAT (3f12.8)
Cabc 9040 FORMAT (f10.5,/,f10.5,2f15.10)
Cabc 9050 FORMAT (i3,/,2d15.8,/,2i2)

      end





      subroutine FINDFORMAT( iof, FORM_JULICH, FORM_MUNCHEN )
c     =====================

c Finds out whether unit IOF contains potential in Munchen or Julich format

      implicit none

      integer iof
      logical FORM_JULICH, FORM_MUNCHEN 

      integer I, iprint
      double precision rdum(9)
      logical found

      data iprint / 0 /

Cabc      iprint = 1

      FORM_JULICH = .false.
      FORM_MUNCHEN  = .false.

c Zkousim, jestli je to Mnichovsky format

      FORM_MUNCHEN = .true.

      call locateline( iof, 'EF  ', 4, found )
      if ( .not.found ) FORM_MUNCHEN  = .false.
      if ( iprint .ge. 1 ) then
        write (*,*) 'Munchen test:   EF --> ', found
      endif 

      call locateline( iof, 'SITES', 5, found )
      if ( .not.found ) FORM_MUNCHEN  = .false.
      if ( iprint .ge. 1 ) then
        write (*,*) 'Munchen test:   SITES --> ', found
      endif 

      call locateline( iof, 'TYPES', 5, found )
      if ( .not.found ) FORM_MUNCHEN  = .false.
      if ( iprint .ge. 1 ) then
        write (*,*) 'Munchen test:   TYPES --> ', found
      endif 

      call locateline( iof, 'POTENTIAL', 9, found )
      if ( .not.found ) FORM_MUNCHEN  = .false.
      if ( iprint .ge. 1 ) then
        write (*,*) 'Munchen test:   POTENTIAL --> ', found
      endif 


      if ( FORM_MUNCHEN ) return 


c Tries whether it is a Julich format
      
      FORM_JULICH = .true.

      call locateline( iof, ' DOWN ', 6, found )
      if ( .not.found ) FORM_JULICH  = .false.
      if ( iprint .ge. 1 ) then
        write (*,*) 'Julich test:   DOWN --> ', found
      endif 
Cabc      write (*,*) 'After DOWN:  FORM_JULICH  = ', FORM_JULICH

      call locateline( iof, ' UP ', 4, found )
      if ( .not.found ) FORM_JULICH  = .false.
      if ( iprint .ge. 1 ) then
        write (*,*) 'Julich test:   UP --> ', found
      endif 
Cabc      write (*,*) 'After UP:  FORM_JULICH  = ', FORM_JULICH

      rewind(iof)

      read (iof,*,err=66,end=66)

      read (iof,*,err=66,end=66) rdum(1), rdum(2), rdum(3)
      do i = 1, 3
        if ( rdum(i).le.0.0 .or. rdum(i).ge.1000.0 ) then
          FORM_JULICH  = .false.
          if ( iprint .ge. 1 ) then
            if ( i.eq.1 ) then
              write (*,*) 'Julich test:    -->  RMT failed'
            else if ( i.eq.2 ) then
              write (*,*) 'Julich test:    -->  ALAT failed'
            else if ( i.eq.3 ) then
              write (*,*) 'Julich test:    -->  RMTNEW failed'
            endif 
          endif 
        endif 
      enddo 
Cabc      write (*,*) 'After RMT, ALAT, RMTNEW:  FORM_JULICH  = ',
Cabc     $     FORM_JULICH

      read (iof,*,err=66,end=66)
      read (iof,*,err=66,end=66)
      read (iof,*,err=66,end=66)

      read (iof,*,err=66,end=66) rdum(1), rdum(2)
      do i = 1, 2
        if ( rdum(i).le.0.0 .or. rdum(i).ge.1000.0 ) then
          FORM_JULICH  = .false.
          if ( iprint .ge. 1 ) then
            if ( i.eq.1 ) then
              write (*,*) 'Julich test:    -->  A failed'
            else if ( i.eq.2 ) then
              write (*,*) 'Julich test:    -->  B failed'
            endif 
          endif 
        endif 
      enddo 
Cabc      write (*,*) 'After A, B:  FORM_JULICH  = ', FORM_JULICH

      return 


 66   FORM_JULICH  = .false.
      return 

      end



      subroutine copy_before_fermi( ifil_old, ifil_new )
c     ============================

C Reads and copies all lines BEFORE the part containing the potentials

      implicit none

      integer ifil_old, ifil_new 

      integer i, iend
      character text*128

      rewind( ifil_old )


 10   do i = 1, len(text)
        text(i:i) = ' '
      enddo 

        read(ifil_old,'(a)',end=61,err=62) text
        if ( text(1:4).eq.'EF  ' ) goto 90

        call findend( text, iend)
        write (ifil_new,'(a)') text(1:iend)


      goto 10
 90   continue 


      return 
      
 61   stop 'EOF before EF section identified !'
 62   stop 'I/O ERR while scanning for EF target !'

      end




      subroutine copy_between_fermi_and_potential( ifil_old, ifil_new )
c     ===========================================

C Reads and copies all lines BEFORE the part containing the potentials

      implicit none

      integer ifil_old, ifil_new 

      integer i, iend
      character text*128


      rewind( ifil_old )


c Finds EF target

 10   do i = 1, len(text)
        text(i:i) = ' '
      enddo 

        read(ifil_old,'(a)',end=61,err=62) text
        if ( text(1:4).eq.'EF  ' ) goto 90

      goto 10
 90   continue 


c Copies lines up to 'POTENTIAL'

 110  do i = 1, len(text)
        text(i:i) = ' '
      enddo 

         read(ifil_old,'(a)',end=161,err=162) text
         call findend( text, iend)
         write (ifil_new,'(a)') text(1:iend)

         if ( text(1:9).eq.'POTENTIAL' ) goto 190

      goto 110
 190  continue 


      return 
      
 61   stop 'EOF while scanning for AFTER EF target !'
 62   stop 'I/O ERR while scanning for AFTER EF target !'      
 161  stop 'EOF before POTENTIAL section identified !'
 162  stop 'I/O ERR while scanning for POTENTIAL target !'
      
      end




      subroutine copy_after_potential( ifil_old, ifil_new )
c     ===============================

C Reads and copies all lines BEFORE the part containing the potentials

      implicit none

      integer ifil_old, ifil_new 

      integer i, iend
      character text*128


      read (ifil_old,*)              ! Skips one '========' line


 10   do i = 1, len(text)
        text(i:i) = ' '
      enddo 

        read(ifil_old,'(a)',end=90,err=62) text
        call findend( text, iend)
        write (ifil_new,'(a)') text(1:iend)

      goto 10
 90   continue 


      return 
      
 62   stop 'I/O ERR while scanning for POTENTIAL target !'

      end




      subroutine readnextpot( iof, POTFMTINP, ITDAT, JTOP, VT, BT )
c     ======================

c Reads SPR-KKR potential

      integer ITDAT, JTOP(*), iof, POTFMTINP
      double precision vt(*), bt(*)

      integer I, iend
      character text*128, FMT07*16


 10   do i = 1, len(text)
        text(i:i) = ' '
      enddo 

      read (iof,'(a)',end=61,err=62) text
      call findend( text, iend)
      
      if ( text(1:6) .eq. 'TYPE  ' ) then
        write (55,'(2x,a,2x)') text(7:iend)
        backspace(55)
        read (55,*) ITDAT
      else 
        goto 10
      endif 


      IF ( POTFMTINP.EQ.6 ) THEN
        FMT07 = '(1P,5E16.9)'
      else 
        FMT07 = '(1P,5E22.14)'
      endif 

      READ (IOF,FMT=FMT07) ( VT(I), I = 1, JTOP(ITDAT) )
      READ (IOF,FMT=FMT07) ( BT(I), I = 1, JTOP(ITDAT) )

      
      return 

 61   stop 'EOF while scanning for TYPE !'
 62   stop 'I/O ERR while scanning for TYPE !'

      end




      subroutine writepot( IOF, IT, POTFMTINP, JTOP, VT, BT )
c     ===================

c Reads SPR-KKR potential

      integer IT, JTOP, iof, POTFMTINP
      double precision VT(*), BT(*)

      integer i
      character FMT07*16


      IF ( POTFMTINP.EQ.6 ) THEN
        FMT07 = '(1P,5E16.9)'
      else 
        FMT07 = '(1P,5E22.14)'
      endif 


      WRITE (IOF,99017) 'TYPE      ',IT
      WRITE (IOF,FMT=FMT07) (VT(I),I=1,JTOP)
      WRITE (IOF,FMT=FMT07) (BT(I),I=1,JTOP)
      
      WRITE (IOF,99011)


      return 
99017 FORMAT (A10,6I10)
99011 FORMAT (79('='))
      end




      subroutine FINDEND( TEXT, IEND )
c     ==================
c
      implicit none
c
c Finds the "logical end" of a character string, i.e. the position
c   of the last non-space or non-blank character.
c
      integer IEND
      character*(*) TEXT
c
      integer len, DIMTEXT, I
c
      DIMTEXT = len( TEXT )
c
      do I = DIMTEXT, 1, -1
        if ( TEXT(I:I) .ne. ' ' ) then
          IEND = I
          goto 1
        endif
      enddo
      IEND = 0
 1    continue
c
      return
      end




      subroutine COUNTARG2( TEXT, NARG, N1 )
c     ====================
c
c Pocita pocet cisel ("argumentu") NARG na zadane radce.
c   V retezci N1(i) jsou pozice zacatku jednotlivych argumentu.
c
c   Called by:  READCTRL
c
      implicit none
c
      integer NARG, N1(*)
      character*(*) TEXT
c
      integer NLEN, I, NN
c
      NLEN = len( TEXT )
c
      NARG = 0
      do I = NLEN-1, 1, -1
        if ( TEXT(I+1:I+1) .eq. ' '  .and.  TEXT(I:I) .ne. ' ' ) then
          NARG = NARG + 1
        endif
      enddo
c
      NN = 0
      if ( TEXT(1:1) .ne. ' ' ) then
        NN = NN + 1
        N1(NN) = 1
      endif
c
      do I = 2, NLEN
        if ( TEXT(I-1:I-1) .eq. ' '  .and.  TEXT(I:I) .ne. ' ' ) then
          NN = NN + 1
          N1(NN) = I
        endif
      enddo
c
      if ( NARG .ne. NN ) then
        write(6,*) ' Source code error in  COUNTARG2  !!'
        write(6,*) '    NARG =', NARG, ',   NN =', NN
        write(6,*) '    TEXT =''', TEXT(1:NLEN), ''''
        stop
      endif
c
      return
      end




      subroutine FINDSIZE( TEXT, IBEG, IEND )
c     ===================
c
      implicit none
c
c Finds the "logical beginning" and the "logical end" of a character string,
c   i.e. the position of the first and of the last non-space or non-blank
c   character.
c If there are no non-blank characters in the string, then IBEG=DimTEXT+1
c   and IEND=0.
c
      integer IBEG, IEND
      character*(*) TEXT
c
      integer len, DIMTEXT, I
c
      DIMTEXT = len( TEXT )
c
      do I = 1, DIMTEXT
        if ( TEXT(I:I) .ne. ' ' ) then
          IBEG = I
          goto 1
        endif
      enddo
      IBEG = DIMTEXT + 1
 1    continue
c
      do I = DIMTEXT, 1, -1
        if ( TEXT(I:I) .ne. ' ' ) then
          IEND = I
          goto 2
        endif
      enddo
      IEND = 0
 2    continue
c
      return
      end
      



      subroutine QINT2( FX, X, NX, FY, Y, NY )
c     ================
c
      implicit none
c      
      integer IY, IX, NX, NY
      double precision   X(*), Y(*), YY, A1, A2, A3, A4
      double precision  FX(*), FY(*), A12, A13, A23, A24, F1, F2
c
c Piecewise quadratic interpolation from grid X to grid Y, by Aitken's 
c   divided difference scheme.  Note NX, NY are array dimensions.
c
c The output grid Y(NY) is unchanged in any circumstance.
c
c Interpolace se provadi dvakrat (ze dvou predchazejicich a jednoho
c   nasledujiciho bodu a z jednoho predchazejiciho a dvou nasledujicich
c   bodu) a z toho se dela prumer.
c
      if ( NX .lt. 3 ) then
        write(6,'(//,a,i5,/,a,/)')
     7       'QINT2 input error:   Too small NX =', NX,
     &       '   -->  Correct calling subroutine and re-run !'
        stop 'Program ended BADLY !'
      endif

c The X-grid has to be increasing otherwise it crashes...

      do IX = 2, NX
        if ( X(IX) .le. X(IX-1) ) stop 'QINT2:  Decreasing X-grid'
      enddo  

      IY=1
      do 2 IX=3,NX-1
1     YY=Y(IY)
      if ( YY .gt. X(IX) ) goto 2
      A1=X(IX-2)-YY
      A2=X(IX-1)-YY
      A3=X(IX)-YY
      A4=X(IX+1)-YY
      A12=(FX(IX-2)*A2-FX(IX-1)*A1)/(X(IX-1)-X(IX-2))
      A13=(FX(IX-2)*A3-FX(IX)*A1)/(X(IX)-X(IX-2))
      A23=(FX(IX-1)*A3-FX(IX)*A2)/(X(IX)-X(IX-1))
      A24=(FX(IX-1)*A4-FX(IX+1)*A2)/(X(IX+1)-X(IX-1))
      F1=(A12*A3-A13*A2)/(X(IX)-X(IX-1))
      if ( YY .ge. X(2) ) then
        F2=(A23*A4-A24*A3)/(X(IX+1)-X(IX))
      else
        F2=F1
      endif
      FY(IY)=(F1+F2)/2.
      IY=IY+1
      if ( IY .gt. NY ) goto 3
      goto 1
   2  continue
c
      IX=NX-1
  10  YY=Y(IY)
      A2=X(IX-1)-YY
      A3=X(IX)-YY
      A4=X(IX+1)-YY
      A23=(FX(IX-1)*A3-FX(IX)*A2)/(X(IX)-X(IX-1))
      A24=(FX(IX-1)*A4-FX(IX+1)*A2)/(X(IX+1)-X(IX-1))
      FY(IY)=(A23*A4-A24*A3)/(X(IX+1)-X(IX))
      IY=IY+1
      if ( IY .gt. NY ) goto 3
      goto 10
c
   3  continue
c
      return
      end
      



      SUBROUTINE READKWINT(IFIL,KW,I,DEF,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFIL                 *
C   *  at the beginning of a line and read following  INTEGER I        *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER DEF,I,IFIL,IFLAG,IWRI
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      INTEGER LL
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFIL
 100  CONTINUE
      READ (IFIL,'(A)',ERR=200,END=200) LINE
      CALL REPLACE_TAB(LINE)
      LL = 80
      CALL SKPLBLK(LINE,LL)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,I
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFIL,STR10,I
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFIL
         STOP
      ELSE
         I = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFIL,STR10,I
         REWIND IFIL
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFIL,KW,LINE
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' INTEGER  input: ',A,I5)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFIL',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,A,/,1X,79('*'),/)
      END




      SUBROUTINE READKWREAL(IFIL,KW,A,DEF,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFIL                 *
C   *  at the beginning of a line and read following  REAL    A        *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      REAL*8 A,DEF
      INTEGER IFIL,IFLAG,IWRI
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      INTEGER LL
      CHARACTER*10 STR10
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFIL
 100  CONTINUE
      READ (IFIL,'(A)',ERR=200,END=200) LINE
      CALL REPLACE_TAB(LINE)
      LL = 80
      CALL SKPLBLK(LINE,LL)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,A
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFIL,STR10,A
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFIL
         STOP
      ELSE
         A = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFIL,STR10,A
         REWIND IFIL
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFIL,KW
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' REAL     input: ',A,E15.8)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFIL',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,1X,79('*'),/)
      END




      SUBROUTINE REPLACE_TAB(STRING)
C   ********************************************************************
C   *                                                                  *
C   *  replace horizontal TABs  by blanks                              *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,L
C
C*** End of declarations rewritten by SPAG
C
      L = LEN_TRIM(STRING)
C
      DO I = 1,L
         IF ( IACHAR(STRING(I:I)).EQ.9 ) STRING(I:I) = ' '
      END DO
      END



      SUBROUTINE SKPLBLK(STRING,LL)
C   ********************************************************************
C   *                                                                  *
C   *  skip leading blanks in   STRING                                 *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER LL
      CHARACTER*(*) STRING
C
C Local variables
C
      INTEGER I,J,K
C
C*** End of declarations rewritten by SPAG
C
      DO I = 1,LL
         IF ( STRING(I:I).NE.' ' ) THEN
            IF ( I.EQ.1 ) RETURN
            K = 0
            DO J = I,LL
               K = K + 1
               STRING(K:K) = STRING(J:J)
            END DO
            LL = LL - I + 1
            RETURN
         END IF
      END DO
      END




      SUBROUTINE POSFIL(IFIL,STR,IPOS,ISTOP)
C   ********************************************************************
C   *                                                                  *
C   *   position file  IFIL  to line starting with string  STR         *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
      INTEGER IFIL,IPOS,ISTOP
C
      CHARACTER*10 STRIN,STR
 
      REWIND IFIL
 
 100  CONTINUE
      READ (IFIL,FMT='(A)',END=200) STRIN
      IF ( STR.NE.STRIN ) GOTO 100
      IPOS = 1
      RETURN
 
 200  CONTINUE
      IPOS = 0
      IF ( ISTOP.NE.0 ) THEN
         WRITE (6,*) ' STOP IN <POSFIL>'
         WRITE (6,*) ' STRING ',STR,' NOT FOUND IN FILE ',IFIL
         STOP
      ELSE
         REWIND IFIL
         RETURN
      END IF
      END





      subroutine locateline( iof, text, lentext, found )
c     =====================

c Finds out whether unit IOF contain string TEXT

      implicit none

      integer iof, lentext
      character text*(*)
      logical found

      character line*128

      found = .false.

      rewind( iof )
      
 10   read (iof,'(a)',end=90) line

      if ( index( line, text(1:lentext) ) .gt. 0 ) then
        found = .true.
        return 
      endif 

      goto 10

 90   return 

      end




      SUBROUTINE READKWSTR(IFIL,KW,STR,DEF,L,IWRI,IFLAG)
C   ********************************************************************
C   *                                                                  *
C   *  find the keyword   KW (STR*10)  in  file   IFIL                 *
C   *  at the beginning of a line and read following  STRING STR*(L)   *
C   *  DEF        default value                                        *
C   *  IWRI  = 1  write result to chanel IWRI                          *
C   *          0  no output                                            *
C   *  IFLAG = 1  stop if KW not found                                 *
C   *          0  continue even if KW not found and set to default     *
C   *                                                                  *
C   ********************************************************************
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER L
      CHARACTER*(L) DEF,STR
      INTEGER IFIL,IFLAG,IWRI
      CHARACTER*10 KW
C
C Local variables
C
      CHARACTER*80 LINE
      INTEGER LL
      CHARACTER*10 STR10
      CHARACTER*70 STR70
C
C*** End of declarations rewritten by SPAG
C
      REWIND IFIL
 100  CONTINUE
      READ (IFIL,'(A)',ERR=200,END=200) LINE
      CALL REPLACE_TAB(LINE)
      LL = 80
      CALL SKPLBLK(LINE,LL)
      IF ( INDEX(LINE,KW).NE.1 ) GOTO 100
C
      READ (LINE,*,ERR=300,END=300) STR10,STR70
      STR = STR70(1:MIN(L,70))
      IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFIL,STR10,STR
      RETURN
C
 200  CONTINUE
      IF ( IFLAG.NE.0 ) THEN
         WRITE (6,99002) KW,IFIL
         STOP
      ELSE
         STR = DEF
         IF ( IWRI.NE.0 ) WRITE (IWRI,99001) KW,IFIL,STR10,STR
         REWIND IFIL
         RETURN
      END IF
C
 300  CONTINUE
      WRITE (6,99003) IFIL,KW
      STOP
99001 FORMAT (2X,'keyword ',A,' file',I4,' STRING   input: ',2A)
99002 FORMAT (//,1X,79('*'),/,2X,'keyword  ',A,
     &        ' not found in file IFIL',I4,/,2X,
     &        'the STOP flag was set indicating input obligatory',/,
     &        79('*'),/)
99003 FORMAT (//,1X,79('*'),/,2X,'error reading file',I4,
     &        ' for keyword ',A,/,1X,79('*'),/)
      END
