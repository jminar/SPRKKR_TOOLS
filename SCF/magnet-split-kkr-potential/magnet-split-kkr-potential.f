      program magnet_split_kkr_potential
c     ==================================

      implicit none

      integer narg, n1(1001), iend, i, iargc, iprint, jp, itc(99),
     $     nput, jws, itx, jpx, ibeg  
      double precision  VT(1001), BT(1001), fac(99), add(99)
      logical poten, nasel
      character*128 inpfile, outfile
      character*128 text, textuc
      CHARACTER*40 FMT07

      data  JWS / 721 /

c -------------------------------------------------------------------------
c Input potential format
c -------------------------------------------------------------------------

      FMT07 = '(1P,5E22.14)'



      iprint = 1

      if ( iargc() .ne. 2 ) then
        write (*,*) 'Two command-line arguments needed:'
        write (*,*) '-->   <Old_potfile>   <New_potfile>'
        stop 'Program ended BADLY !'
      endif  


c Opens the files 

      call getarg( 1, inpfile )
      open (15, file=inpfile, status='old',action='read')

      call getarg( 2, outfile )
      open (16, file=outfile, status='unknown',action='write')
      
      open(55,status='scratch' )
 
c Browses given sites and reads magnetic splitting factor

      jp = 0

 10   continue 
        do i = 1, len(text)
          text(i:i) = ' '
        enddo 
        write (6,'(9(a))') 'Enter:   ',
     $       '<IT>      [ +<addmag> | *<facmag>]    (or EOF) !'
        read (5,'(a)',end=90) text

        if ( text(1:1) .eq. '#'  .or.  text(1:1) .eq. '!' ) goto 10

        call countarg2( text, narg, n1 )
        call findend( text, iend )

        
        if ( narg.eq.1 ) then
          call FINDSIZE(text,ibeg,iend)
          textuc(1:iend) = text(1:iend)
          call CNVTUC( textuc, 1, iend )
          if ( index(textuc(ibeg:iend),'EOF') .gt. 0 .or.
     $         index(textuc(ibeg:iend),'END') .gt. 0 ) then
            goto 90
          endif 
        endif 

        if ( narg.ne.2 ) then
          write (*,*) 'Illegal format of the input line:  NARG =', narg
          write (*,*) 'Input =''', text(1:iend), ''''
          stop 'Program ended BADLY !'
        endif 

        jp = jp + 1

        write (55,'(a)') text(1:n1(2)-1)
        backspace(55)
        read (55,*) itc(jp)
        if ( itc(jp) .le. 0 .or. itc(jp) .gt. 2001 ) then
          write (*,*) 'Illegal ITC  !'
          write (*,*) 'Input line: ''', text(1:iend), ''''
          stop 'Program ended BADLY !'
        endif 

        if ( text(n1(2):n1(2)) .eq. '+' ) then
          fac(jp) = 1.0d0
          write (55,'(9(a))') '  ', text(n1(2)+1:iend), '  '
          backspace(55)
          read (55,*) add(jp)

        else if ( text(n1(2):n1(2)) .eq. '*' ) then
          add(jp) = 0.0d0
          write (55,'(9(a))') '  ', text(n1(2)+1:iend), '  ' 
          backspace(55)
          read (55,*) fac(jp)

        else 
          write (*,*) 'Illegal FAC or ADD format  !'
          write (*,*) 'Input line: ''', text(1:iend), ''''
          stop 'Program ended BADLY !'
        endif 

      goto 10

 90   nput = jp

      
c Initial setting of special flags

      poten = .false.

c Reads host input potential line by line 

 100  do I = 1, len(text)
        text(i:i) = ' '
      enddo 

      read(15,'(a)',end=900) text

      call FINDEND(text,iend)

c -------------------------------------------------------------------------
c Looks for potential initialization targets

      if ( text(1:9) .eq. 'POTENTIAL' ) then
        poten = .true.        
        write(16,'(a)') text(1:iend)
        goto 100
      endif 

c No search for atomic types in CHARGE, MOMENTS and CORE sections
      
      if ( poten .and. text(1:6) .eq. 'CHARGE' ) then
        poten = .false.
        write(16,'(a)') text(1:iend)
        goto 100
      endif 
      if ( text(1:7) .eq. 'MOMENTS' ) then
        poten = .false.
        write(16,'(a)') text(1:iend)
        goto 100
      endif       
      if ( text(1:4) .eq. 'CORE' ) then
        poten = .false.
        write(16,'(a)') text(1:iend)
        goto 100
      endif 

c -------------------------------------------------------------------------
c Looks if it specifies one of the potentials to be kicked

      if ( poten ) then
        if ( text(1:5) .eq. 'TYPE ' ) then
          write(55,'(a)') text(5:len(text))
          backspace(55)
          read(55,*) itx
          nasel = .false.

          do jp = 1, nput
            if ( itx .eq. itc(jp) ) then
              
c Je to typ, ktery ma byt vymenen

              READ (15,FMT=FMT07) ( vt(i), I = 1, JWS )
              READ (15,FMT=FMT07) ( bt(i), I = 1, JWS )

              jpx = jp
              nasel = .true.
              goto 200
            endif 
          enddo 

c Vypise radek s cislem typu atomu 

 200      call findend(text,iend)
          write(16,'(a)') text(1:iend)
          
c Neni to jeden z typu ktere budou vymeneny
          
          if ( .not.nasel ) then
            goto 100            
          endif 

c Alters the magnetic part of the potential

          do i = 1, jws
            BT(I) = BT(I) + add(jpx)
            BT(I) = fac(jpx) * BT(I)
          enddo 
                
          write (16,FMT=FMT07) ( VT(I), I = 1, JWS )
          write (16,FMT=FMT07) ( BT(I), I = 1, JWS )      

          if ( iprint .ge. 1 ) then
            write (*,'(a,i4,a,f8.4,a,f8.4)')
     $           'Potential kicked for type', itx,
     $           '   MagADD =', add(jpx), '   MagFAC =', fac(jpx)
          endif 

          goto 100
          
c Neni to radek s oznacenim typu atomu

        else           
          call findend(text,iend)
          write(16,'(a)') text(1:iend)
          goto 100
        endif 

      endif         

c ---------------------------------------------------------
c Copies the line into new file 
      
      call findend(text,iend)
      write(16,'(a)') text(1:iend)
      goto 100

c The "old" or host input file has been finished

 900  continue 

      stop
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



      SUBROUTINE CNVTUC(STRING,I1,I2)
C   ********************************************************************
C   *                                                                  *
C   *  convert characters to upper case in  STRING(I1:I2)              *
C   *                                                                  *
C   ********************************************************************
C
      IMPLICIT NONE
C
C*** Start of declarations rewritten by SPAG
C
C Dummy arguments
C
      INTEGER I1,I2
      CHARACTER*(*) STRING
C
C Local variables
C
      CHARACTER CHAR
      INTEGER I,IA,INC,IZ,J
      INTEGER ICHAR
C
C*** End of declarations rewritten by SPAG
C
      IA = ICHAR('a')
      IZ = ICHAR('z')
      INC = ICHAR('A') - IA
C
      DO I = I1,I2
         J = ICHAR(STRING(I:I))
         IF ( J.GE.IA .AND. J.LE.IZ ) STRING(I:I) = CHAR(J+INC)
      END DO
      END
