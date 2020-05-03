c
c Read-in SPRKKR x-ray absorption rate output (*.rat files),
C and deploy sum rules originally from: "X-ray circular dichroism 
C as a probe of orbital magnetization" by B.T.Thole, P.Carra, F.Sette, G.van der Laan,
C PRL68,1943 (1992), https://doi.org/10.1103/PhysRevLett.68.1943 
C or:
C "Synchrotron Radiation Techniques Based on X-ray Magnetic Circular Dichroism" 
C by G.Schütz, E.Goering, H.Stoll, Handbook of Magnetism and Advanced Magnetic 
C Materials (2007), https://doi.org/10.1002/9780470022184.hmm304
C
C ... in order to provide estimates of the effective (S_z together with T_z) spin
C and orbital magnetic moment
C
c Originally by Ondrej Sipr, Division of Solid State Physics, Department of 
C Structure Analysis, Czech Academy of Science Praha: sipr@fzu.cz
c
c Minor revisions: A.Marmodoro, LMU Muenchen, July 2018 (automatic arguments checks)
C                                         February 2019 (padding and effective EF
C                                                        shifting for non-SCF NEGF
C                                                        usage)
C                                         August 2019 (preliminary addition of M_{4,5}
C                                                      absorption edge; array sizes)
C

C
C Hard-coded maximal problem size. Modify and recompile if ever needed:
C
#define _NTMAX_  151
#define _NQMAX_  151
#define _NEMAX_  6001

C
C Collecting auxiliary procedures into modules allows the compiler
C to ensure each subroutine / function gets the expected arguments:
C
      MODULE text_utilities
C
C Convenience subroutines, to parse ASCII text from plain *.rat files:
C
      IMPLICIT NONE
C
      CONTAINS
C
c **********************************************************************
      SUBROUTINE FINDEND( TEXT, IEND )
c **********************************************************************
c
      IMPLICIT NONE
c
c Finds the "LOGICAL END" of a CHARACTER string, 
C i.e. the position of the last non-space or non-blank character.
c
      INTEGER,INTENT(OUT) :: IEND
      CHARACTER*(*),INTENT(IN) :: TEXT
c
      INTEGER LEN, DIMTEXT, I
c
      DIMTEXT = LEN( TEXT )
C
c Loop from the end of the string:
C
      DO I = DIMTEXT, 1, -1
         IF ( TEXT(I:I) .NE. ' ' ) THEN
            IEND = I
            GOTO 1
         ENDIF
      ENDDO
      IEND = 0
 1    CONTINUE
c
      END SUBROUTINE FINDEND
C
      SUBROUTINE TXTTEXT( TXTT,LTXTT, NT, NTMAX )
c **********************************************************************
c * Add an extension '_a', '_b' etc. to the element name TXTT, in case *
C * there are more types with the same atomic number Z                 *
c **********************************************************************
      INTEGER,INTENT(IN) :: NT, NTMAX
      CHARACTER(LEN=10),INTENT(INOUT) :: TXTT(NTMAX)
      INTEGER,INTENT(OUT) :: LTXTT(NTMAX)
C
      INTEGER ICAL, ICAU, ICZL, ICZU, IT, I1, I2, iext, jt, ic2,
     &   KDONE(NTMAX)
C
      ICAL = ICHAR('a')
      ICAU = ICHAR('A')
      ICZL = ICHAR('z')
      ICZU = ICHAR('Z')

      DO IT=1,NT
         TXTT(IT) = TXTT(IT)(1:2)//'        '
         IC2 = ICHAR(TXTT(IT)(2:2))

         IF ( ((IC2.LT.ICAL) .OR. (IC2.GT.ICZL)) .AND.
     &       ((IC2.LT.ICAU) .OR. (IC2.GT.ICZU)) ) THEN
            TXTT(IT)(2:2) = ' '
            LTXTT(IT) = 1
         ELSE
            LTXTT(IT) = 2
         ENDIF  
         KDONE(IT) = 0
      ENDDO

      DO IT = 1,NT
         I1=LTXTT(IT)+1
         I2=I1+1
         IEXT = 0
         IF ( KDONE(IT) .NE. 1 ) THEN
            DO JT=1,NT
               IF ( (IT.NE.JT) .AND. 
     &             (TXTT(IT)(1:2).EQ.TXTT(JT)(1:2)) ) THEN
                  TXTT(IT)(I1:I2) = '_a'
                  IEXT = IEXT + 1
                  TXTT(JT)(I1:I2) = '_'//CHAR( ICHAR('a') + IEXT )
                  KDONE(JT) = 1
                  LTXTT(JT)=LTXTT(JT)+2
               ENDIF
            ENDDO
         ENDIF
         IF ( IEXT .NE. 0 ) LTXTT(IT)=LTXTT(IT)+2
         KDONE(IT) = 1
C
      ENDDO
C
      END SUBROUTINE TXTTEXT
C
      END MODULE text_utilities
C
C =========================================================
C
      MODULE math_utilities
C
C Case-specific interpolation scheme, and trapezoidal
C quadrature rule for numerical integration.
C
      IMPLICIT NONE
      CONTAINS
C
      SUBROUTINE QINT2( FX, X, NX, FY, Y_single, NY )
c
c ******************************************************************
c * Piecewise QUADRATIC interpolation from grid X to grid Y, using *
C * Aitken's divided difference scheme.  Note NX, NY are the array *
C * dimensions.                                                    *
c *                                                                *
c * The output grid Y(NY) is unchanged in any circumstance.        *
c *                                                                *
c * Czech: Interpolace se provadi dvakrat (ze dvou predchazejicich *
c * a jednoho nasledujiciho bodu a z jednoho predchazejiciho       *
C * a dvou nasledujicich bodu) a z toho se dela prumer.            *
c *                                                                *
c * Interpolation is performed as an average of two quadratic      *
C interpolations, that proceed left-to-right and right-to-left.    *
c ******************************************************************
C
      IMPLICIT NONE
C
C Portable double precision definition, used everywhere:
C
      INTEGER,PARAMETER :: DP = selected_REAL_kind(15, 307)
C      
      INTEGER,INTENT(IN) :: NX, NY
C
      REAL(DP),INTENT(IN) :: X(NX), FX(NX), Y_single
      REAL(DP),INTENT(OUT) ::  FY(NY)
C
      REAL(DP) YY, A1, A2, A3, A4,
     &             A12, A13, A23, A24, F1, F2
      INTEGER IY, IX
C
C Safeguard:
C
      IF ( NX .LT. 3 ) THEN
         WRITE(6,'(//,a,i5,/,a,/)')
     7       'QINT2 input error:   Too small NX =', NX,
     &       '   -->  Correct calling SUBROUTINE and re-run !'
         STOP 'PROGRAM ended BADLY !'
      ELSEIF ( NY.NE.1 ) THEN
         STOP 'QINT2 assuming NY=1, aborting.'
      ENDIF
C
c The X-grid values have to be monotonically increasing, 
C otherwise it crashes...
C
      DO IX = 2, NX
         IF ( X(IX) .LE. X(IX-1) ) THEN
            WRITE(6,*) 'IX-1 =', IX-1, '      X(IX-1) =', X(IX-1)
            WRITE(6,*) 'IX   =', IX,   '      X(IX)   =', X(IX)
            STOP 'QINT2:  Decreasing X-grid values, aborting'
         ENDIF
      ENDDO
C
      IY = 1
      loop_IX : DO IX=3,NX-1
1       CONTINUE 
C
        YY = Y_single
C
        IF ( YY .GT. X(IX) ) GOTO 2
C
        A1 = X(IX-2) - YY
        A2 = X(IX-1) - YY
        A3 = X(IX) - YY
        A4 = X(IX+1) - YY
C
        A12 = ( FX(IX-2)*A2-FX(IX-1)*A1 ) / ( X(IX-1)-X(IX-2) )
        A13 = ( FX(IX-2)*A3-FX(IX)*A1 ) / ( X(IX)-X(IX-2) )
        A23 = ( FX(IX-1)*A3-FX(IX)*A2 ) / ( X(IX)-X(IX-1) )
        A24 = ( FX(IX-1)*A4-FX(IX+1)*A2 ) / ( X(IX+1)-X(IX-1) )
C
        F1 = ( A12*A3-A13*A2 ) / ( X(IX)-X(IX-1) )
    
        IF ( YY .GE. X(2) ) THEN
          F2 = ( A23*A4 - A24*A3 ) / ( X(IX+1)-X(IX) )
        ELSE
          F2 = F1
        ENDIF
C
        FY(IY) = (F1+F2) / 2.0_DP
C
        IY = IY+1
        IF ( IY .GT. NY ) GOTO 3
        GOTO 1
2       CONTINUE
      ENDDO loop_IX
c
      IX = NX-1
  10  CONTINUE
C
      YY = Y_single
C
      A2 = X(IX-1)-YY
      A3 = X(IX)-YY
      A4 = X(IX+1)-YY
      A23 = ( FX(IX-1)*A3-FX(IX)*A2 ) / ( X(IX)-X(IX-1) )
      A24 = ( FX(IX-1)*A4-FX(IX+1)*A2 ) / ( X(IX+1)-X(IX-1) )
C
      FY(IY) = ( A23*A4-A24*A3 ) / ( X(IX+1)-X(IX) )
C
      IY = IY+1
      IF ( IY .GT. NY ) GOTO 3
      GOTO 10
c
   3  CONTINUE
c
      END SUBROUTINE QINT2
C      
      PURE FUNCTION TRAPEZ( X0, Y0, N_samples, X, Y )
c *********************************************************
c * Trapezoidal rule integration of function Y=F(X)       *
c *********************************************************
c
      IMPLICIT NONE
c
C Portable double precision definition, used everywhere:
      INTEGER,PARAMETER :: DP = selected_REAL_kind(15, 307)
      REAL(DP),PARAMETER :: half = 0.5_DP
C
      INTEGER,INTENT(IN) :: N_samples
      REAL(DP),INTENT(IN) :: x0, y0, x(N_samples), y(N_samples)
      REAL(DP) trapez
C
      INTEGER i
C
      TRAPEZ = half * ( X(1) - X0 ) * ( Y(1) + Y0 )
      DO i = 2, N_samples
         TRAPEZ = TRAPEZ + half * (Y(i) + Y(i-1)) * (X(i) - X(i-1))
      ENDDO
C
      END FUNCTION trapez
C
      END MODULE math_utilities
C
C ==================================================================
C
      MODULE io_utilities
C
      IMPLICIT NONE
C
      CONTAINS
C
c **********************************************************************
c Reads the data from a SPRKKR x-ray absorption rate '*.rat' output file
C and apply ADDITIONAL padding etc.
c **********************************************************************
C
      SUBROUTINE READ_RXAS( ratfile_io,EMESH_FROM_RATFILE, 
     &   XMU_PL_L3, XMU_MN_L3, 
     &   XMU_Z0_L3,XMU_PL_L2, XMU_MN_L2, XMU_Z0_L2, ECORE,
     &   NE_FROM_RATFILE, NTXG, NCST, NQT, CONC, TXTT, LTXTT,
     &   ecoreav, NKTYP, EBAND,
     &   set_to_zero_below_EF, add_padded_samples_below_EF,
     &   kedge,l23edge,m45edge)
C
      USE text_utilities,ONLY: FINDEND,TXTTEXT
C
      IMPLICIT NONE
C
      INTEGER,PARAMETER :: IP = 0,
     &   NE_MAX_PADDED = 50 
C
C Portable double precision definition, used everywhere:
      INTEGER,PARAMETER :: DP = selected_REAL_kind(15, 307)
C
      REAL(DP),PARAMETER :: RYOVEV = 13.605826_DP,
     &                      RZERO = 0.0_DP
C
      CHARACTER(LEN=2),PARAMETER :: EUNIT = 'eV'
C
      INTEGER,INTENT(IN) :: ratfile_io
      REAL(DP),INTENT(IN) :: EBAND
      LOGICAL,INTENT(IN) :: set_to_zero_below_EF,
     &                    add_padded_samples_below_EF
      REAL(DP),INTENT(OUT) :: EMESH_FROM_RATFILE(_NEMAX_),CONC(_NTMAX_),
     &   XMU_PL_L3(_NEMAX_), XMU_MN_L3(_NEMAX_), XMU_Z0_L3(_NEMAX_),
     &   XMU_PL_L2(_NEMAX_), XMU_MN_L2(_NEMAX_), XMU_Z0_L2(_NEMAX_),
     &   ecoreav(9),
     &   ECORE(_NTMAX_,10) ! 10 for the 5 'd' initial states of M_{4,5} edge;
                           !  6 for the 3 'p' initial states of L_{2,3} edge.
      INTEGER,INTENT(OUT) :: NE_FROM_RATFILE,NKTYP,NCST,NTXG,
     &   NQT(_NTMAX_),LTXTT(_NTMAX_)
      CHARACTER(LEN=10),INTENT(OUT) :: TXTT(_NTMAX_)
      LOGICAL,INTENT(OUT) :: kedge,l23edge,m45edge
C
      INTEGER NQ, NT, IREL, IQP, IQ, NLQ(_NQMAX_), ITP,
     &   IQOFT(_NQMAX_,_NTMAX_),
     &   NLT(_NTMAX_), ITXG, ITXRSGRP(_NTMAX_), NCXRAY(_NTMAX_), 
     &   LCXRAY(_NTMAX_), IDIPOL, NPOL, I, ICST, ITXGP,
     &   NE0, IFMT, KUNIT, IT, IE, 
     &   IE_closest_to_EF, NE_PADDED, 
     &   JE, II, I1, I2, IQT, 
     &   LHEADER, K,LSYSTEM, IPOL, JCST, ktyp(9), 
     &   NKPCOR(10), !  6 for the 3 'p' initial states of L_{2,3} edge;
                     ! 10 for the 5 'd' initial states of M_{4,5} edge.
     &   KTYPCOR(10), KAPCOR(10), MM05COR(10), IKMCOR(10,2)
      REAL(DP) RD(_NEMAX_,_NTMAX_,10,3), !  6 for the 3 'p' initial states of L_{2,3} edge;
                                         ! 10 for the 5 'd' initial states of M_{4,5} edge.
     &  EREFER, dist, wgt, 
     &     XNORM(6), VUC, SUM, DELE, SUMK(9),EF
      COMPLEX(DP) ERYD, ABSRTA, ABSRTB, ABSRATE
      CHARACTER(LEN=32) KEYWORD
      CHARACTER(LEN=80) HEADER, INFO, SYSTEM, FMT710, FMT740 ! LINE, 
      LOGICAL fneg
C
C Safeguard: Check that the *.rat input really is a SPRKKR x-ray absorption rate file:
C
      READ(ratfile_io,'(10X,A)') KEYWORD
      IF ( KEYWORD(1:4) .NE. 'RXAS' ) THEN
         WRITE(6,*) 'KEYWORD =''', KEYWORD, ''''
         STOP 'Keyword RXAS not found ! Valid *.rat file?'
      ENDIF
C
c Basic dimensions and specifiers:
C
      READ(ratfile_io,'(10X,A80)') HEADER
      READ(ratfile_io,'(10X,A80)') SYSTEM
      READ(ratfile_io,'(10X,I10)') NQ
      READ(ratfile_io,'(10X,I10)') NT
      READ(ratfile_io,*) 
      READ(ratfile_io,'(10X,I10)') NE_FROM_RATFILE
      READ(ratfile_io,'(10X,I10)') IREL
      READ(ratfile_io,'(10X,F10.5)') EF
      READ(ratfile_io,'(10X,A80)') INFO
      READ(ratfile_io,*)
C
c Check the dimensions for the statically allocated arrays or die:
C
      IF ( NQ.GT._NQMAX_ ) THEN
         WRITE(6,'("Parsed NQ=",i0," but _NQMAX_=",i0,
     &      ": too small: increase it within sum_rules_theory.f")')
     &   NQ,_NQMAX_
C
      ELSEIF (NT.GT._NTMAX_ ) THEN
         WRITE(6,'("Parsed NT=",i0," but _NTMAX_=",i0,
     &      ": too small: increase it within sum_rules_theory.f")')
     &   NT,_NTMAX_
C
      ELSEIF (NE_FROM_RATFILE.GT._NEMAX_ ) THEN
         WRITE(6,'("Parsed NE_FROM_RATFILE=",i0," but _NEMAX_=",i0,
     &      ": too small, increase it within sum_rules_theory.f")')
     &   NE_FROM_RATFILE,_NEMAX_
C
      ELSE
C All OK, proceed further:
         GOTO 10
C
      ENDIF 
      STOP 'PROGRAM ended BADLY! Arrays size problem?'
C
10    CONTINUE
! Debug:
      WRITE(6,'("From ratfile: NE_FROM_RATFILE=",i0)') NE_FROM_RATFILE
      WRITE(6,'("From ratfile: EF=",F15.8," [Ry] relative to VMTZ")') 
     &   EF
C ----------------------------------------------------------------------
c Reads atomic sites (IQ) -resolved data:
C
      READ(ratfile_io,*)
C Angular momentum expansion, site IQ -resolved:
      DO IQP = 1,NQ
!          READ(ratfile_io,'(A)') LINE
!          WRITE(88,'(A)') LINE(2:80)
!          BACKSPACE(88)
!          READ(88,*) IQ, NLQ(IQ)
!          READ(ratfile_io,'(A)') LINE(2:80)
!          WRITE(88,'(A)') LINE(2:80)
!          BACKSPACE(88)
! Simpler:
         READ(ratfile_io,*) IQ, NLQ(IQ)
      ENDDO
C Debug:
      WRITE(6,'("From ratfile: NLQ(IQ=",i0,")=",i0)') (IQ,NLQ(IQ),
     &   IQ=1,NQ)
C
Cabc      READ(ratfile_io,*)
Cabc      DO ITP = 1, NT
Cabc         READ(ratfile_io,'(1X,I4,1X,A10,F10.5,I5,10I3,:,(41X,10I3))') IT,
Cabc     &       TXTT(IT), CONC(IT), NQT(IT),
Cabc     &       ( IQOFT(IQT,IT), IQT = 1, NQT(IT) )
Cabc         NLT(IT) = NLQ(IQOFT(1,IT))
Cabc      ENDDO
C
C ----------------------------------------------------------------------
c Reads atomic types (IT) -resolved data:
      READ(ratfile_io,'(10X,A80)')
C
C Name, Concentration, angular momentum expansion, occupied site:
      DO ITP = 1,NT
         READ(ratfile_io,*) IT,TXTT(ITP),CONC(ITP),NQT(ITP),
     &                     (IQOFT(IQT,ITP),IQT=1,NQT(ITP))
C
C Different types on the same site share the same angular momentum cutoff:
         NLT(ITP) = NLQ(IQOFT(1,ITP))
      ENDDO
C Debug:
      WRITE(6,'("From ratfile: NLT(ITP=",i0,")=",i0)') (ITP,NLT(ITP),
     &   ITP=1,NT)
C
      CALL TXTTEXT( TXTT, LTXTT, NT, _NTMAX_ )
C
C Safeguard:
C
      IF ( IP .GE. 2 ) THEN
         CALL FINDEND( HEADER, LHEADER )
         CALL FINDEND( SYSTEM, LSYSTEM )
         WRITE(6,*) ' HEADER:     ', HEADER(1:LHEADER)
         WRITE(6,*) ' SYSTEM:     ', SYSTEM(1:LSYSTEM)
!          WRITE(6,*) '<RDHEAD>  passed '
      ENDIF 
C
      READ(ratfile_io,'(10x,I10)') NTXG
      READ(ratfile_io,'(10x,I10)') IDIPOL
      READ(ratfile_io,'(10x,I10)') NPOL
C
C Verify that rates are given for +, - and 0 x-ray polarizations:
C
      IF ( NPOL .NE. 3 ) THEN
         WRITE(6,'("NPOL=",i0,
     &      " but it has to be 3 in this version, aborting.")') NPOL
      ENDIF
C
C ======================================================================
c Reads HEADER for first spectrum in the *.rat file
C
C Loop over types of atoms:
C
      LOOP_READ_ITXG: DO ITXG = 1, NTXG
C
C Recognize which absorption edge is contained:
C
         kedge = .FALSE.
         l23edge = .FALSE.
         m45edge = .FALSE.
C
         READ(ratfile_io,'(a)') KEYWORD
C
         IF ( ITXG .EQ. 1 ) THEN
            IF ( INDEX(KEYWORD,' K') .GE. 3  .AND.
     &           INDEX(KEYWORD,' K') .LE. LEN(KEYWORD)-1 ) THEN
               kedge = .TRUE.
C
            ELSEIF ( INDEX(KEYWORD,' L2,3') .GE. 3  .AND.
     &               INDEX(KEYWORD,' L2,3') .LE. LEN(KEYWORD)-5 ) THEN
               l23edge = .TRUE.
C
            ELSEIF ( INDEX(KEYWORD,' M4,5') .GE. 3  .AND.
     &               INDEX(KEYWORD,' M4,5') .LE. LEN(KEYWORD)-5 ) THEN
               m45edge = .TRUE.
               WRITE(6,'("WARNING: support for M4,5 absorption edge",
     &            " is EXPERIMENTAL: press ENTER to continue...")')
               READ(*,*)
C
            ELSE
               WRITE(6,*) 'KEYWORD =''', KEYWORD, ''''
               STOP 'X-ray edge not identified/supported, aborting!'
            ENDIF
C
         ELSE
C
C Once recognized, the same absorption edge is expected for any subsequent ITXG:
C
            IF ( kedge ) THEN
               IF ( INDEX(KEYWORD,' K') .LT. 3  .OR.
     &              INDEX(KEYWORD,' K') .GT. LEN(KEYWORD)-1 ) THEN
                 STOP 'K-edge signature lost for further ITXG !'
               ENDIF
C
            ELSEIF ( l23edge ) THEN
               IF ( INDEX(KEYWORD,' L2,3') .LT. 3  .AND.
     &              INDEX(KEYWORD,' L2,3') .GT. LEN(KEYWORD)-5 ) THEN
                  STOP 'L2,3-edge signature lost for further ITXG !'
               ENDIF
C
            ELSEIF ( m45edge ) THEN
               IF ( INDEX(KEYWORD,' M4,5') .LT. 3  .AND.
     &              INDEX(KEYWORD,' M5,5') .GT. LEN(KEYWORD)-5 ) THEN
                  STOP 'M4,5-edge signature lost for further ITXG !'
               ENDIF
C
            ELSE
               STOP 'X-ray signature not known for further ITXG !'
            ENDIF
         ENDIF
C
         READ(ratfile_io,'(10x,I10)') ITXRSGRP(ITXG)
         READ(ratfile_io,'(10x,I10)') NCXRAY(ITXG)
         READ(ratfile_io,'(10x,I10)') LCXRAY(ITXG)

         IF ( ITXG .EQ. 1 ) THEN
            IF ( IP .GE. 2 ) THEN
               WRITE(6,102) IDIPOL, NCXRAY(ITXG), LCXRAY(ITXG), NPOL
            ENDIF
         ELSE
            ITXGP = ITXG-1
C Safeguard:
            IF ( NCXRAY(ITXG) .NE. NCXRAY(ITXGP) .OR.
     &         LCXRAY(ITXG) .NE. LCXRAY(ITXGP) ) THEN
               WRITE(6,*) ' NC:  ',NCXRAY(ITXG),NCXRAY(ITXGP)
               WRITE(6,*) ' LC:  ',LCXRAY(ITXG),LCXRAY(ITXGP)
               WRITE(6,*) ' ITXG:',       ITXG ,       ITXGP 
               STOP ' data inconsistent in <READ_RXAS>, aborting '
            ENDIF 
         ENDIF
C
         READ(ratfile_io,'(////,8X,I4,/,8X,20I4)',ERR=999) NCST,
     &      ( NKPCOR(ICST), ICST = 1, NCST )
         DO i = 1, 2
            READ(ratfile_io,*)
         ENDDO
         IF ( IP .GE. 1 ) THEN
            WRITE(6,980) ITXG
            IF ( IP .GE. 2 ) THEN
               WRITE(6,981) NCST, ( NKPCOR(ICST), ICST = 1, NCST )
               WRITE(6,982)
            ENDIF
         ENDIF
C
c ----------------------------------------------------------------------
c Reads core energies:
C
C Debug:
         WRITE(6,'("Parsing core energies for ITXG=",i0,
     &      " NCST=",i0)') ITXG,NCST
C
         LOOP_ICST: DO ICST = 1, NCST 
                         
            READ(ratfile_io,'(5I4,2X,I4,F12.6,F12.4,F12.3)') JCST,
     &          NCXRAY(ITXG), LCXRAY(ITXG), KAPCOR(ICST),
     &          MM05COR(ICST), IKMCOR(ICST,1), XNORM(1), 
     &          ECORE(ITXG,ICST)
C
C Safeguard: verify the correct core state line is parsed:
            IF ( JCST .NE. ICST ) STOP 'JCST .NE. ICST !'
C
            MM05COR(ICST) =  (MM05COR(ICST)-1)/2
C
            IF ( IP .GE. 1 ) THEN
               WRITE(6,984) ICST, NCXRAY(ITXG), LCXRAY(ITXG),
     &            KAPCOR(ICST),
     &            (2*MM05COR(ICST)+1), IKMCOR(ICST,1),XNORM(1), 
     &            ECORE(ITXG,ICST),
     &            ECORE(ITXG,ICST) * RYOVEV
            ENDIF
C
            IF ( NKPCOR(ICST) .EQ. 2 ) THEN 
               READ(ratfile_io,'(22X,I4,F12.6)') IKMCOR(ICST,2), 
     &            XNORM(2)
               IF ( IP .GE. 1 ) THEN
                  WRITE(6,986) IKMCOR(ICST,2), XNORM(2)
               ENDIF 
            ENDIF
C
C Put into consistent units:
C
            IF ( EUNIT .EQ. 'eV' ) THEN
               ECORE(ITXG,ICST) = ECORE(ITXG,ICST) * RYOVEV
            ENDIF
C
c Debug: report what has been parsed:
            WRITE(6,'(" ECORE(ITXG=",i0,",ICST=",i0,")=",F15.8," ",A)') 
     &         ITXG,ICST,ECORE(ITXG,ICST),EUNIT
         ENDDO LOOP_ICST
C
C ----------------------------------------------------------------------
C
         READ(ratfile_io,'(/,20X,I4,/,19X,2I5)') NE0, IFMT, KUNIT
C
         FMT710='(2X,2F 7.4,4X,F8.4, 12X,2E14.6 )'
         FMT740='(I2,I2,1X,2E11.4,1X,2E11.4,2X,2E11.4 )'         
         IF ( IFMT .GE. 1 ) THEN
            FMT710='(2X,2F11.8,4X,F8.4, 12X,2E14.6 )'
            FMT740='(I2,I2,1X,2E15.8,1X,2E15.8,2X,2E15.8 )'
         ENDIF
C
         IF ( IFMT .GE. 2 ) THEN
            READ(ratfile_io,'(10X,F15.8)',err=61) VUC
            GOTO 62
 61         VUC = RZERO
 62         CONTINUE
         ELSE
            VUC = RZERO
         ENDIF  
         IF ( IFMT .NE. 0 ) READ(ratfile_io,*)
C
c ----------------------------------------------------------------------
c Reads the spectral transition rate for each energy ERYD
C
         DO IE = 1, NE0
C
C The absorption rate calculation energy can be complex:
C
            READ(ratfile_io,FMT=FMT710,ERR=999,END=997) ERYD

            LOOP_READ_ICST: DO ICST = 1,NCST
               LOOP_READ_IPOL: DO IPOL = 1,NPOL
C
                 READ(ratfile_io,FMT=FMT740,ERR=999) 
     &               I1, I2, ABSRTA, ABSRTB, ABSRATE
C
C Safeguard: verify indexing for core states and polarizations:
                  IF ( I1 .NE. ICST ) STOP ' ( I1 .NE. ICST ) !'
                  IF ( I2 .NE. IPOL ) STOP ' ( I2 .NE. IPOL ) !'
C
C Store the real part of the (complex) absorption rate ABSRTA:
C
                  RD(IE,ITXG,ICST,IPOL) = DBLE( ABSRTA )
C
               ENDDO LOOP_READ_IPOL
            ENDDO LOOP_READ_ICST
C
            IF ( ITXG .EQ. 1 ) THEN
C
C For the first type, as the energy mesh is then anyway kept the same:
C store the real part of each (complex) sampling point:
C
               EMESH_FROM_RATFILE(IE) = DBLE( ERYD )
C
            ELSE
C
C Safeguard: verify that the parsed energy mesh remains consistent
C for subsequent ITXG:
              IF ( ABS(DBLE(ERYD) - EMESH_FROM_RATFILE(IE)) .GT.
     &           0.000001_DP ) THEN
                STOP 'Inconsistent energy ranges for further ITXG !'
              ENDIF
C
            ENDIF 
            READ(ratfile_io,*,END=997)
C
C Last available point, i.e. how many in total:
            NE_FROM_RATFILE = IE
C
         ENDDO
c ----------------------------------------------------------------------
  997    CONTINUE
C
C Safeguard:
         IF ( NE_FROM_RATFILE .NE. NE0 ) THEN
            WRITE(6,*) ' only NE_FROM_RATFILE=', NE_FROM_RATFILE,
     &          ' E-values tabulated, instead of NE0=',NE0
            IF ( NTXG .NE. 1 ) STOP 'NTXG > 1 may cause trouble !'
         ENDIF
C
      ENDDO LOOP_READ_ITXG
C
      IF ( IP .GE. 1 ) THEN
         WRITE(6,'(i5,a,f10.5,a,f10.5,A)') NE_FROM_RATFILE, 
     &      ' energy points read: from EMESH_FROM_RATFILE(1)=', 
     &      EMESH_FROM_RATFILE(1),
     &      ' to EMESH_FROM_RATFILE(NE_FROM_RATFILE)=', 
     &      EMESH_FROM_RATFILE(NE_FROM_RATFILE), ' [Ry]'
      ENDIF 
C
C Find KTYP for each core state:
C
      NKTYP   = 0
      KTYP(1) = 0
      DO ICST=1,NCST  
         DO K=1,NKTYP
            IF ( KAPCOR(ICST) .EQ. KTYP(K) ) THEN
              KTYPCOR(ICST) = K
              GOTO 210
            ENDIF
         ENDDO
C
         NKTYP         = NKTYP + 1
         IF ( NKTYP .GT. 2 ) THEN
            STOP ' NKTYP > 2  >>>>>> check core states'
         ENDIF
         KTYP(NKTYP)   = KAPCOR(ICST)
         KTYPCOR(ICST) = NKTYP
C
 210     CONTINUE
      ENDDO
C
c Calculates average core energies for each KAPPA:
C
      DO K = 1, NKTYP
         SUMK(K)   = RZERO
         ECOREAV(K)= RZERO
      ENDDO
C
      DO ITXG = 1, NTXG
C
C Weight reflecting the concentration:
C
         WGT = DBLE(NQT(ITXG)) * CONC(ITXG)
         DO ICST=1,NCST
            K          = KTYPCOR(ICST) 
            ECOREAV(K) = ECOREAV(K) + WGT * ECORE(ITXG,ICST) 
            SUMK(K)    = SUMK(K)    + WGT
         ENDDO
      ENDDO
C
      DO K=1,NKTYP
         ECOREAV(K)= ECOREAV(K) / SUMK(K) 
      ENDDO
C
C Few more checks that everything is all-right:
C
      IF ( IP .GE. 2 ) THEN
         WRITE(6,'(A,F10.4,:,6x,F10.4,1x,A)') ' original  E-range: ',
     &      EMESH_FROM_RATFILE(1), 
     &      EMESH_FROM_RATFILE(NE_FROM_RATFILE), EUNIT
      ENDIF 
C
      DO ITXG = 2,NTXG
         IF ( ITXRSGRP(ITXG) .LE. ITXRSGRP(ITXG-1) ) THEN
           WRITE(*,*) 'ITXRSGRP(I) <= ITXRSGRP(I-1) for I=',ITXG
            STOP
         ENDIF
      ENDDO
C
C Safeguard: monotonically increasing energies?
C
      DO IE = 2, NE_FROM_RATFILE
        IF ( EMESH_FROM_RATFILE(IE) - EMESH_FROM_RATFILE(IE-1) 
     &     .LE. RZERO ) THEN
          WRITE(6,*) 'IE-1 =', IE-1, '   EMESH_FROM_RATFILE(IE-1) =', 
     &       EMESH_FROM_RATFILE(IE-1)
          WRITE(6,*) 'IE   =', IE,   '   EMESH_FROM_RATFILE(IE)   =', 
     &       EMESH_FROM_RATFILE(IE)
          STOP '(1)  Energies are not sorted in ascending order !'
        ENDIF 
      ENDDO 
C
c ----------------------------------------------------------------------
c Post-processing the parsed data:
C
C Put energy points into consistent units:
C
      IF ( EUNIT .EQ. 'eV' ) THEN
         DO IE = 1, NE_FROM_RATFILE
            EMESH_FROM_RATFILE(IE) = EMESH_FROM_RATFILE(IE) * RYOVEV
         ENDDO
      ENDIF
C
C Put the Fermi energy too into consistent units:
C
      IF ( EUNIT .EQ. 'eV' ) THEN
         EF = EF * RYOVEV
      ENDIF
      WRITE(6,'("Parsed Fermi energy: EF=",F15.8," ",A,
     &   " relative to VMTZ")') EF,EUNIT
C
      SUM = RZERO
      DO ITXG = 1,NTXG
        IT = ITXRSGRP(ITXG)
        IF ( IT .LT. ITXG ) STOP 'IT < ITXG'
        IF ( IP .GE. 2 ) THEN
          WRITE(*,*) 'adjusting table for ITXG IT:', ITXG, IT
        ENDIF
        LTXTT(ITXG) = LTXTT(IT)
        TXTT(ITXG)  = TXTT(IT)
        NQT(ITXG)   = NQT(IT)
        CONC(ITXG)  = CONC(IT)
        SUM = SUM + CONC(IT)
      ENDDO
C
C Verify / enforce that the concentration adds up to 100%:
C
      IF ( ABS(SUM) .LT. 0.000001_DP ) THEN
        WRITE(6,'("WARNING: setting CONC(:) to 100%")')
        DO ITXG = 1,NTXG
          CONC(ITXG)  = 1.0_DP
        ENDDO
      ENDIF
C
c Removes spaces from atomic labels TXTT:
C
      DO ITXG = 1, NTXG
 41     CONTINUE 
        DO i = 1, LTXTT(ITXG)
C
          IF ( TXTT(ITXG)(i:i) .EQ. ' ' ) THEN
            DO ii = i, LTXTT(ITXG) - 1
              TXTT(ITXG)(ii:ii) = TXTT(ITXG)(ii+1:ii+1) 
            ENDDO
            LTXTT(ITXG) = LTXTT(ITXG) - 1
            GOTO 41
          ENDIF 
C
        ENDDO
      ENDDO
C
C ----------------------------------------------------------------------
C Debug:
      WRITE(6,'("Before referring to EF as new energy zero: ",
     &   "EMESH_FROM_RATFILE(IE=",i0,")=",
     &   F15.8)') (IE,EMESH_FROM_RATFILE(IE),IE=1,NE_FROM_RATFILE)
C
C Adopt the Fermi energy as reference zero:
C
      EREFER = EF
      EF = EF - EREFER
C
c Shifts the energy zero to Fermi energy:
C
C 1) for the core states' energies:
C
      DO ITXG=1,NTXG
         DO ICST=1,NCST 
            ECORE(ITXG,ICST) = ECORE(ITXG,ICST) - EREFER
         ENDDO
      ENDDO
C
      DO K = 1, NKTYP
         ECOREAV(K) = ECOREAV(K) - EREFER
      ENDDO
C
C 2) for the final states energies parsed:
C
      DO IE = 1, NE_FROM_RATFILE
         EMESH_FROM_RATFILE(IE) = EMESH_FROM_RATFILE(IE) - EREFER 
      ENDDO
C
C Safeguard: again verify monotonically increasing energies:
C
      DO IE = 2, NE_FROM_RATFILE
         IF ( EMESH_FROM_RATFILE(IE) - EMESH_FROM_RATFILE(IE-1) 
     &      .LE. RZERO ) THEN
            STOP '(2)  Energies are not sorted in ascending way !'
         ENDIF
      ENDDO
C
c Writes polarized absorption coefficient 
c   assigned to each edge:
C
C Safeguard: does the number of core states match expectation for this edge?
C
      IF ( l23edge.and.(NCST.NE.6) ) THEN
         WRITE(6,*) 'Illegal NCST =', NCST
         WRITE(6,*) ' for l23edge=',l23edge
         STOP 'PROGRAM ended BADLY ! l23edge implies NCST=6'
C
      ELSEIF ( m45edge.and.(NCST.NE.10) ) THEN
         WRITE(6,*) 'Illegal NCST =', NCST
         WRITE(6,*) ' for m45edge=',m45edge
         STOP 'PROGRAM ended BADLY ! m45edge implies NCST=10'
      ENDIF
C
      IF ( NTXG .NE. 1 ) THEN
         WRITE(6,*) 'Illegal NTXG =', NTXG
         WRITE(6,*) 'This version assumes NTXG=1 always'
         STOP 'PROGRAM ended BADLY !'
      ENDIF
C
C Ondrej Sipr:
c I am using the opposite +/- convention as Hubert Ebert / SPRKKR
C because then the XMCD can be defined as mu(+) - mu(-) and it has got 
c same signs at L3/L2 edges as in nearly all publications
c (including ours)
C
      DO IE = 1, NE_FROM_RATFILE
C
C Adding up cross sections 
C under assumption that indexing refers to same sampling energy
C regardless of any difference in core level splitting:
C
         XMU_PL_L3(IE) = RZERO
         XMU_MN_L3(IE) = RZERO
         XMU_Z0_L3(IE) = RZERO
C
         XMU_PL_L2(IE) = RZERO
         XMU_MN_L2(IE) = RZERO
         XMU_Z0_L2(IE) = RZERO
C
C WARNING: what about the k edge case, does it make sense to treat it
C          the same way as the l23edge (as the original version was doing)?
         if( l23edge.OR.kedge ) then
C
C Dipole-allowed transitions for the L_{2,3} case:
C
C The L_3 edge has as initial states the four 2P 3/2 levels ICST = 3,4,5,6:
C +---------------------------+
C |Example case of BCC Fe, L_3|
C +---------------------------+--------------------------------------------------
C  ICST  N   L  KAP  MUE  IKM     NORM        E(Ry)       E(eV)    <SIGMA_z>  I0
C
c    3   2   1  -2  -3/2   5    1.000000    -49.6459    -675.467     -0.9975  721
c 
c    4   2   1  -2  -1/2   6    0.998655    -49.6694    -675.786     -0.4006  721
c                          3    0.001345
c 
c    5   2   1  -2   1/2   7    0.998504    -49.6940    -676.122      0.2588  721
c                          4    0.001496
c 
c    6   2   1  -2   3/2   8    1.000000    -49.7201    -676.476      0.9975  721
c
c (notice the grouping into four shallower & contiguous initial energies values)
C ---------------------------------------------------------------------------------
            DO ICST = 3, 6
               XMU_PL_L3(IE) = XMU_PL_L3(IE) + RD(IE,1,ICST,2)
               XMU_MN_L3(IE) = XMU_MN_L3(IE) + RD(IE,1,ICST,1)
               XMU_Z0_L3(IE) = XMU_Z0_L3(IE) + RD(IE,1,ICST,3)
            ENDDO
C ===============================================================================
C +---------------------------+
C |Example case of BCC Fe, L_2|
C +---------------------------+--------------------------------------------------
C The L_2 edge has as initial states the two 2P 1/2 levels ICST = 1,2:
c  ICST  N   L  KAP  MUE  IKM     NORM        E(Ry)       E(eV)    <SIGMA_z>  I0
c    1   2   1   1  -1/2   3    0.998655    -50.6114    -688.603      0.3993  721
c                          6    0.001345
c
c    2   2   1   1   1/2   4    0.998504    -50.5873    -688.275     -0.2575  721
c                          7    0.001496
c (notice the grouping into two deeper & contiguous initial energies values)
C ---------------------------------------------------------------------------------
            DO ICST = 1, 2
               XMU_PL_L2(IE) = XMU_PL_L2(IE) + RD(IE,1,ICST,2)
               XMU_MN_L2(IE) = XMU_MN_L2(IE) + RD(IE,1,ICST,1)
               XMU_Z0_L2(IE) = XMU_Z0_L2(IE) + RD(IE,1,ICST,3)
            ENDDO
C ===============================================================================
C ===============================================================================
         ELSEIF( m45edge ) then
C
C +---------------------------+
C |Example case of HCP Gd, M_5|
C +---------------------------+---------------------------------------------------
C The M_5 edge has as initial states the six 3D 5/2 levels ICST = 5,6,7,8,9,10:
C (Table 4.4 page 71  in "Electronic Structure and magneto-optical properties of 
C solids" by V.Antonov, B.Harmon, A.Yaresko, Kluwer Academics (2004)
C --------------------------------------------------------------------------------
c  ICST  N   L  KAP  MUE  IKM     NORM        E(Ry)       E(eV)    <SIGMA_z>  I0
c    5   3   2  -3  -5/2  13    1.000000    -84.5272   -1150.051     -0.9941  691
c 
c    6   3   2  -3  -3/2  14    0.999365    -84.5562   -1150.446     -0.6358  691
c                          9    0.000635
c 
c    7   3   2  -3  -1/2  15    0.998999    -84.5860   -1150.851     -0.2601  691
c                         10    0.001001
c 
c    8   3   2  -3   1/2  16    0.998946    -84.6165   -1151.266      0.1351  691
c                         11    0.001054
c 
c    9   3   2  -3   3/2  17    0.999258    -84.6478   -1151.692      0.5522  691
c                         12    0.000742
c 
c   10   3   2  -3   5/2  18    1.000000    -84.6800   -1152.130      0.9941  691
c (notice the grouping into six shallower & contiguous initial energies values)
            DO ICST = 5, 10
               XMU_PL_L3(IE) = XMU_PL_L3(IE) + RD(IE,1,ICST,2)
               XMU_MN_L3(IE) = XMU_MN_L3(IE) + RD(IE,1,ICST,1)
               XMU_Z0_L3(IE) = XMU_Z0_L3(IE) + RD(IE,1,ICST,3)
            ENDDO
C ===============================================================================
C +---------------------------+
C |Example case of HCP Gd, M_4|
C +---------------------------+---------------------------------------------------
C The M_4 edge has as initial states the four 3D 3/2 levels ICST = 1,2,3,4:
C (Table 4.4 page 71  in "Electronic Structure and magneto-optical properties of 
C solids" by V.Antonov, B.Harmon, A.Yaresko, Kluwer Academics (2004)
C --------------------------------------------------------------------------------
c  ICST  N   L  KAP  MUE  IKM     NORM        E(Ry)       E(eV)    <SIGMA_z>  I0
c    1   3   2   2  -3/2   9    0.999366    -86.9336   -1182.792      0.6337  691
c                         14    0.000634
c 
c    2   3   2   2  -1/2  10    0.998999    -86.9051   -1182.404      0.2594  691
c                         15    0.001001
c 
c    3   3   2   2   1/2  11    0.998946    -86.8757   -1182.005     -0.1344  691
c                         16    0.001054
c 
c    4   3   2   2   3/2  12    0.999258    -86.8456   -1181.595     -0.5501  691
c                         17    0.000742 
c (notice the grouping into four deeper & contiguous initial energies values)
C ---------------------------------------------------------------------------------
            DO ICST = 1, 4
               XMU_PL_L2(IE) = XMU_PL_L2(IE) + RD(IE,1,ICST,2)
               XMU_MN_L2(IE) = XMU_MN_L2(IE) + RD(IE,1,ICST,1)
               XMU_Z0_L2(IE) = XMU_Z0_L2(IE) + RD(IE,1,ICST,3)
            ENDDO
C ===============================================================================
C ===============================================================================
         ELSE
            STOP 'Unknown absorption edge, cannot prepare XMU_PL/MN/Z0'
         ENDIF
      ENDDO
C
C Debug: raw output file for simple plotting:
C
      DO IE = 1, NE_FROM_RATFILE
         WRITE(73,'(f10.5,4(3x,3f10.5))')
     &      EMESH_FROM_RATFILE(IE), 
     &      RD(IE,1,3,1), RD(IE,1,3,2), RD(IE,1,3,3),
     &      RD(IE,1,4,1), RD(IE,1,4,2), RD(IE,1,4,3),
     &      RD(IE,1,5,1), RD(IE,1,5,2), RD(IE,1,5,3),
     &      RD(IE,1,6,1), RD(IE,1,6,2), RD(IE,1,6,3)
         WRITE(72,'(f10.5,4(3x,3f10.5))')
     &      EMESH_FROM_RATFILE(IE), 
     &      RD(IE,1,1,1), RD(IE,1,1,2), RD(IE,1,1,3),
     &      RD(IE,1,2,1), RD(IE,1,2,2), RD(IE,1,2,3)
      ENDDO
C
C First find out the energy sample closest to the value parsed for the Fermi level:
C
      DIST = 1000000000.0_DP
      DO IE = 1, NE_FROM_RATFILE
         IF ( ABS(EMESH_FROM_RATFILE(IE) - EF) .LT. DIST ) THEN 
            DIST = ABS(EMESH_FROM_RATFILE(IE) - EF)
            IE_closest_to_EF = IE
         ENDIF
      ENDDO
      WRITE(6,'(" EMESH_FROM_RATFILE(IE=",i0,")=",F15.8,
     &   " lies closest to EF=",F15.8," ",A)') 
     &   IE_closest_to_EF,EMESH_FROM_RATFILE(IE_closest_to_EF),EF,EUNIT
! Debug:
!       WRITE(6,'(" Previous and next samples at: ",
!      &   " EMESH_FROM_RATFILE(IE_closest_to_EF-1)=",F15.8,
!      &   " EMESH_FROM_RATFILE(IE_closest_to_EF+1)=",F15.8)')
!      &   EMESH_FROM_RATFILE(IE_closest_to_EF-1),
!      &   EMESH_FROM_RATFILE(IE_closest_to_EF+1)
C      
      IF ( .NOT.add_padded_samples_below_EF ) THEN
         WRITE(6,'("NOT adding any samples below EF")')
         GOTO 555
      ELSE
         WRITE(6,'("Adding up to NE_MAX_PADDED=",i0,
     &      " padding samples")') NE_MAX_PADDED
      ENDIF
C
c Extend E-range below E_F
C Debug:
      WRITE(6,'("Before extending: EMESH_FROM_RATFILE(IE=",i0,")=",
     &    F15.8)') (IE,EMESH_FROM_RATFILE(IE),IE=1,NE_FROM_RATFILE)
C
C Find out how many additional energy points lie below those
C which were actually sampled, despite lying below the Fermi level:
C
      NE_PADDED   = NE_MAX_PADDED - IE_closest_to_EF + 1
      IF ( NE_PADDED .LE. 0 ) THEN
         STOP ' ( NE_PADDED .LE. 0 )  !'
      ELSE
         WRITE(6,'("NE_PADDED=",i0," out of NE_MAX_PADDED=",i0)') 
     &      NE_PADDED,NE_MAX_PADDED
      ENDIF
C
C Arbitrarily chosen Delta_Energy?
C
      DELE = 40.0_DP*DBLE(NE_PADDED)/DBLE(NE_MAX_PADDED)
C
C Shift all the actually sampled values, to make room for the padding:
C
      DO JE = NE_FROM_RATFILE, 1, -1
         IE    = JE + NE_PADDED
         EMESH_FROM_RATFILE(IE) = EMESH_FROM_RATFILE(JE)
         WRITE(6,'("Overwriting IE=",i0," with JE=",i0)') IE,JE
C
         XMU_PL_L3(IE) = XMU_PL_L3(JE)
         XMU_MN_L3(IE) = XMU_MN_L3(JE)
         XMU_Z0_L3(IE) = XMU_Z0_L3(JE)
         XMU_PL_L2(IE) = XMU_PL_L2(JE)
         XMU_MN_L2(IE) = XMU_MN_L2(JE)
         XMU_Z0_L2(IE) = XMU_Z0_L2(JE)
      ENDDO
C
C Apply the padding:
C
      DO IE = 1, NE_PADDED
         EMESH_FROM_RATFILE(IE) = EMESH_FROM_RATFILE(NE_PADDED+1) - 
     &     DELE*(DBLE(NE_PADDED-IE + 1)/DBLE(NE_PADDED))**3
C
C In all cases, 'fake' samples are set to zero:
         XMU_PL_L3(IE) = RZERO
         XMU_MN_L3(IE) = RZERO
         XMU_Z0_L3(IE) = RZERO
         XMU_PL_L2(IE) = RZERO
         XMU_MN_L2(IE) = RZERO
         XMU_Z0_L2(IE) = RZERO
      ENDDO
C
C Account for the offset of padded values:
C
      NE_FROM_RATFILE  = NE_FROM_RATFILE  + NE_PADDED
      IE_closest_to_EF = IE_closest_to_EF + NE_PADDED 
C Debug:
      WRITE(6,'("After adding padded to zero samples: ",
     &   "EMESH_FROM_RATFILE(IE=",
     &   i0,")=",F15.8)') (IE,EMESH_FROM_RATFILE(IE),
     &   IE=1,NE_FROM_RATFILE)
C
555   CONTINUE
C
C Set to 0 all cross sections below  E_F:
C
      IF ( set_to_zero_below_EF ) THEN
         WRITE(6,'("Setting to ZERO all cross section samples",
     &     " below IE_closest_to_EF=",i0," (regardless of padding)")') 
     &     IE_closest_to_EF
         DO IE  = 1, IE_closest_to_EF-1
            XMU_PL_L3(IE) = RZERO
            XMU_MN_L3(IE) = RZERO
            XMU_Z0_L3(IE) = RZERO
            XMU_PL_L2(IE) = RZERO
            XMU_MN_L2(IE) = RZERO
            XMU_Z0_L2(IE) = RZERO
         ENDDO
      ELSE
         WRITE(6,'("NOT setting to zero any data below EF")')
      ENDIF
C
C Safeguard:
C
      DO IE = 2, NE_FROM_RATFILE
         IF ( EMESH_FROM_RATFILE(IE) - EMESH_FROM_RATFILE(IE-1)
     &      .LE. RZERO ) THEN
            STOP '(3)  Energies are not sorted in ascending way !'
         ENDIF 
      ENDDO
C
      IF ( IP .GE. 1 ) THEN
         WRITE(6,9142) EREFER, EF, EMESH_FROM_RATFILE(IE_closest_to_EF),
     &     IE_closest_to_EF
      ENDIF
C
      DO IE = 1, NE_FROM_RATFILE
         IF ( EMESH_FROM_RATFILE(IE) .LT. 2.0_DP*EBAND ) THEN
            fneg = .FALSE.
            IF ( XMU_PL_L3(IE) .LT. RZERO ) THEN
               WRITE(6,'(a,f10.5,a,F15.8)')
     &         'Negative  XAS:   for E =', EMESH_FROM_RATFILE(IE),
     &         '    XMU_PL_L3 =', XMU_PL_L3(IE)
             fneg = .TRUE.
            ENDIF
            IF ( XMU_MN_L3(IE) .LT. RZERO ) THEN
               IF ( .NOT.fneg ) THEN
                  WRITE(6,'(a,f10.5,a,F15.8)')
     &           'Negative  XAS:   for E =', EMESH_FROM_RATFILE(IE),
     &           '    XMU_MN_L3 =', XMU_MN_L3(IE)
               ELSE 
                  WRITE(6,'(24x,10x,a,F15.8)')
     &           '    XMU_MN_L3 =', XMU_MN_L3(IE)
               ENDIF
               fneg = .TRUE.
            ENDIF
            IF ( XMU_Z0_L3(IE) .LT. RZERO ) THEN
               IF ( .NOT.fneg ) THEN
                  WRITE(6,'(a,f10.5,a,F15.8)')
     &              'Negative  XAS:   for E =', EMESH_FROM_RATFILE(IE),
     &              '    XMU_Z0_L3 =', XMU_Z0_L3(IE)
               ELSE
                  WRITE(6,'(24x,10x,a,F15.8)')
     &               '    XMU_Z0_L3 =', XMU_Z0_L3(IE)
               ENDIF
               fneg = .TRUE.
            ENDIF
            IF ( XMU_PL_L2(IE) .LT. RZERO ) THEN
               IF ( .NOT.fneg ) THEN
                  WRITE(6,'(a,f10.5,a,F15.8)')
     &              'Negative  XAS:   for E =', EMESH_FROM_RATFILE(IE),
     &              '    XMU_PL_L2 =', XMU_PL_L2(IE)
               ELSE
                  WRITE(6,'(24x,10x,a,F15.8)')
     &              '    XMU_PL_L2 =', XMU_PL_L2(IE)
               ENDIF
               fneg = .TRUE.
            ENDIF
            IF ( XMU_MN_L2(IE) .LT. RZERO ) THEN
               IF ( .NOT.fneg ) THEN
                  WRITE(6,'(a,f10.5,a,F15.8)')
     &           'Negative  XAS:   for E =', EMESH_FROM_RATFILE(IE),
     &           '    XMU_MN_L2 =', XMU_MN_L2(IE)
               ELSE
                  WRITE(6,'(24x,10x,a,F15.8)')
     &           '    XMU_MN_L2 =', XMU_MN_L2(IE)
               ENDIF
               fneg = .TRUE.
            ENDIF
            IF ( XMU_Z0_L2(IE) .LT. RZERO ) THEN
               IF ( .NOT.fneg ) THEN
                  WRITE(6,'(a,f10.5,a,F15.8)')
     &           'Negative  XAS:   for E =', EMESH_FROM_RATFILE(IE),
     &           '    XMU_Z0_L2 =', XMU_Z0_L2(IE)
               ELSE 
                  WRITE(6,'(24x,10x,a,F15.8)')
     &           '    XMU_Z0_L2 =', XMU_Z0_L2(IE)
               ENDIF 
               fneg = .TRUE.
            ENDIF 
         ENDIF 
      ENDDO
C
      RETURN 
C
  999 STOP 'ERROR reading datafile in <READ_RXAS>'
C
  102 FORMAT(1H1,//,
     & 4X,' IDIPOL:     ',I5,'  0: rate=A+B',/,
     & 4X, '                    1: rate=A',//,
     & 4X,' CORE quantum-numbers  N=',I2,'  L=',I2,/,
     & 4X,' NPOL:       ',I5,/)
 980  FORMAT(//,' CORE STATES :    ITXG =',i3 )
  981 FORMAT(//, ' NCST:  ',I4,/,' NKPCOR:',20I4)
  982 FORMAT( /,' ICST  N   L  KAP  MUE  IKM     NORM   ',
     &          '     EMESH_FROM_RATFILE(Ry)   EMESH_FROM_RATFILE(eV)')
  984 FORMAT(5I4,'/2',I4,F12.6,F12.4,F12.3)
  986 FORMAT(22X,I4,F12.6)      
 9142 FORMAT(/,'  EREFER    ',F10.4,' eV = E_Fermi (rel. VMTZ)',/,
     &     '  EF        ',F10.4,/,
     &   '  EMESH_FROM_RATFILE(IE_closest_to_EF)    ',
     &   F10.4,/,'  IE_closest_to_EF       ',I10,/)
C
      END SUBROUTINE READ_RXAS
C
      END MODULE io_utilities
C
c **********************************************************************
      PROGRAM sum_rules_theory
c **********************************************************************
      USE math_utilities,ONLY: QINT2,trapez
      USE io_utilities,ONLY: READ_RXAS
      IMPLICIT NONE
C Portable double precision definition:
      INTEGER,PARAMETER :: DP = selected_REAL_kind(15, 307)
C
      REAL(DP),PARAMETER :: RZERO=0.0_DP
      INTEGER,PARAMETER :: ratfile_io=15, 
     &   raw_file_L2_IO=62, raw_file_L3_IO=63
C
C
      INTEGER NQT(_NTMAX_), LTXTT(_NTMAX_), NE_FROM_RATFILE, i, 
     &  IE, NKTYP, IE_first_above_zero, IE_last_below_EBAND,
     &  NE_used, NTXG, NCST
      REAL(DP) EMESH_FROM_RATFILE(_NEMAX_), EBAND,
     &   CONC(_NTMAX_),
     &   XMU_PL_L3(_NEMAX_), XMU_MN_L3(_NEMAX_),XMU_Z0_L3(_NEMAX_), 
     &   XMU_PL_L2(_NEMAX_), XMU_MN_L2(_NEMAX_),XMU_Z0_L2(_NEMAX_),
     &   DIF_L3(_NEMAX_), DIF_L2(_NEMAX_), 
     &   SUM_L3(_NEMAX_), SUM_L2(_NEMAX_), 
     &   DELA2, DELA3, 
     &   EMESH_used(_NEMAX_), 
     &   att, spin, orb,ratio, 
     &   ecoreav(9),
     &   ECORE(_NTMAX_,10) ! 6 for the 3 'p' initial states of L_{2,3} edge, each spin-polarized;
                           ! 10 for the 5 'd' initial states of M_{4,5} edge, each spin-polarized.
      CHARACTER(LEN=10) TXTT(_NTMAX_)
      CHARACTER(LEN=40) ratfile
      CHARACTER(LEN=40) argument
      LOGICAL FILE_EXISTS,set_to_zero_below_EF,
     &             add_padded_samples_below_EF,
     &             kedge,l23edge,m45edge
C
C Usage instructions:
C
!       IF ( iargc() .NE. 4 ) THEN
      IF ( command_argument_count() .NE. 4 ) THEN
         WRITE(6,*) 'Four command-line arguments needed !'
C
         WRITE(6,*) '--> rat_file (IN)'
         WRITE(6,*) '-->E_band-edge ',
     &     '(i.e. energy where IDOS_d_up + IDOS_d_down = 10)'
         WRITE(6,*) '-->set_to_zero_below_EF (T or F)'
         WRITE(6,*) '-->add_padded_samples_below_EF (T or F)'
         STOP 'Insufficient input'
      ENDIF         
C
C SPRKKR XAS/XMCD output file parsing:
C
!       CALL GETARG(1, ratfile )
      CALL get_command_argument(1, argument)
      ratfile = trim(argument)
      INQUIRE (FILE=ratfile,EXIST=FILE_EXISTS)
C
      IF ( .NOT.FILE_EXISTS ) THEN
         WRITE(6,'("Pre-existing ratfile=",a,
     &      " does not exist, aborting.")')
     &      ratfile
         STOP
      ENDIF
C
      OPEN( ratfile_io, file=ratfile, STATUS='OLD',ACTION='READ',ERR=20)
      GOTO 30
C
20    WRITE(6,'("ratfile=",a,
     &   " could not be opened for reading, aborting")')
     &   ratfile
      STOP
C
30    CONTINUE
C
C Cutoff energy i.e. E_band-edge parsing:
C
      CALL get_command_argument(2, argument)
!       CALL GETARG(2, argument )
C
!       OPEN (88,STATUS='SCRATCH')
!       WRITE(88,'(a)') argument
!       BACKSPACE(88)
!       READ(88,*) EBAND
! Simpler:
      READ(argument,*) EBAND
C
      CALL get_command_argument(3, argument)
      READ(argument,*) set_to_zero_below_EF
C
      CALL get_command_argument(4, argument)
      READ(argument,*) add_padded_samples_below_EF
! Debug:
      WRITE(6,'("Working with: set_to_zero_below_EF=",L1,
     &   " add_padded_samples_below_EF=",L1)') 
     & set_to_zero_below_EF,add_padded_samples_below_EF
C
C Read-in SPRKKR *.rat output:
C
C Sane defaults, in all cases:
      EMESH_FROM_RATFILE(1:_NEMAX_) = RZERO
      XMU_PL_L3(1:_NEMAX_) = RZERO
      XMU_PL_L2(1:_NEMAX_) = RZERO
      XMU_MN_L3(1:_NEMAX_) = RZERO
      XMU_MN_L2(1:_NEMAX_) = RZERO
      XMU_Z0_L3(1:_NEMAX_) = RZERO
      XMU_Z0_L2(1:_NEMAX_) = RZERO
C
      CALL READ_RXAS( ratfile_io,EMESH_FROM_RATFILE,
     &   XMU_PL_L3, XMU_MN_L3, XMU_Z0_L3,
     &   XMU_PL_L2, XMU_MN_L2, XMU_Z0_L2, 
     &   ECORE, NE_FROM_RATFILE,
     &   NTXG, NCST, NQT, CONC, TXTT, LTXTT,ecoreav, NKTYP, EBAND,
     &   set_to_zero_below_EF, add_padded_samples_below_EF,
     &   kedge,l23edge,m45edge )
C
      CLOSE(ratfile_io)
C
C Verify ordering of energy points, as parsed from the *.rat file:
C
      DO IE = 2, NE_FROM_RATFILE
        IF ( EMESH_FROM_RATFILE(IE)-EMESH_FROM_RATFILE(IE-1) 
     &     .LE. 0.0_DP ) THEN
          STOP '(4)  Energies are not sorted in ascending way !'
        ENDIF 
      ENDDO  
C
C Determine the integration interval between 1) initial and 2) final samples:
C
C 1) starting point: the first EMESH_FROM_RATFILE(IE) above those set to ZERO 
C within READ_RXAS:
C
      IF ( set_to_zero_below_EF ) THEN
         DO IE = 2, NE_FROM_RATFILE
           IF ( (EMESH_FROM_RATFILE(IE).GT.RZERO).AND.
     &          (EMESH_FROM_RATFILE(IE-1).LE.RZERO) ) THEN
             IE_first_above_zero = IE
             GOTO 210
           ENDIF 
         ENDDO
         STOP 'ERROR:  IE_first_above_zero not found !'
 210     CONTINUE 
      ELSE
         IE_first_above_zero = 1
         WRITE(6,'("NOT skipping any energy point below EF")')
      ENDIF
      WRITE(6,'("IE_first_above_zero=",i0,
     &   " out of NE_FROM_RATFILE=",i0)') 
     &   IE_first_above_zero,NE_FROM_RATFILE
C
C 2) ending point: 
C the last IE still below the EBAND cutoff of e_band-edge:
C
      DO IE = NE_FROM_RATFILE-1, 1, -1
        IF ( ( EMESH_FROM_RATFILE(IE).LT.EBAND ).AND.
     &       ( EMESH_FROM_RATFILE(IE + 1).GE.EBAND ) ) THEN
          IE_last_below_EBAND = IE
          GOTO 220
        ENDIF 
      ENDDO
      STOP 'ERROR:  IE_last_below_EBAND not found !'
 220  CONTINUE
C
      WRITE(6,'("IE_last_below_EBAND=",i0,
     &   " out of NE_FROM_RATFILE=",i0)') 
     &   IE_last_below_EBAND, NE_FROM_RATFILE
C
c Interpolates last spectral point: in case the energy mesh did not contain it,
C a value for each relevant quantities there is estimated anyway.
C
C L_3 edge:
C
      CALL QINT2( XMU_PL_L3, EMESH_FROM_RATFILE, NE_FROM_RATFILE, 
     &   XMU_PL_L3(IE_last_below_EBAND+1), EBAND, 1 )
C     
      CALL QINT2( XMU_MN_L3, EMESH_FROM_RATFILE, NE_FROM_RATFILE, 
     &   XMU_MN_L3(IE_last_below_EBAND+1), EBAND, 1 )
C
      CALL QINT2( XMU_Z0_L3, EMESH_FROM_RATFILE, NE_FROM_RATFILE, 
     &  XMU_Z0_L3(IE_last_below_EBAND+1), EBAND, 1 )
C
C L_2 edge:
C
      CALL QINT2( XMU_PL_L2, EMESH_FROM_RATFILE, NE_FROM_RATFILE, 
     &   XMU_PL_L2(IE_last_below_EBAND+1), EBAND, 1 )
C
      CALL QINT2( XMU_MN_L2, EMESH_FROM_RATFILE, NE_FROM_RATFILE, 
     &   XMU_MN_L2(IE_last_below_EBAND+1), EBAND, 1 )
C
      CALL QINT2( XMU_Z0_L2, EMESH_FROM_RATFILE, NE_FROM_RATFILE, 
     &   XMU_Z0_L2(IE_last_below_EBAND+1), EBAND, 1 )
C
C Set up the actual portion of the energy mesh which will be used in the integrations:
C
      NE_used = IE_last_below_EBAND - IE_first_above_zero + 2
C
      DO IE = 1, NE_used-1
         EMESH_used(IE) = EMESH_FROM_RATFILE(IE+IE_first_above_zero-1)
      ENDDO
C
      EMESH_used(NE_used) = EBAND
C
C Evaluates FUNCTIONs which will be integrated: sums and differences
C
      DO IE = IE_first_above_zero, IE_last_below_EBAND+1
         DIF_L3(IE-IE_first_above_zero+1) = XMU_PL_L3(IE) 
     &                                    - XMU_MN_L3(IE)
         DIF_L2(IE-IE_first_above_zero+1) = XMU_PL_L2(IE)
     &                                    - XMU_MN_L2(IE)
C
         SUM_L3(IE-IE_first_above_zero+1) = XMU_PL_L3(IE) 
     &                                    + XMU_MN_L3(IE)
     &                                    + XMU_Z0_L3(IE)
         SUM_L2(IE-IE_first_above_zero+1) = XMU_PL_L2(IE) 
     &                                    + XMU_MN_L2(IE)
     &                                    + XMU_Z0_L2(IE)
      ENDDO
C
C Simple shift of the values, consistent with the above, just for debugging:
C
      DO IE = IE_first_above_zero, IE_last_below_EBAND+1
         XMU_PL_L3(IE-IE_first_above_zero+1) = XMU_PL_L3(IE)
         XMU_MN_L3(IE-IE_first_above_zero+1) = XMU_MN_L3(IE)
         XMU_Z0_L3(IE-IE_first_above_zero+1) = XMU_Z0_L3(IE)
C
         XMU_PL_L2(IE-IE_first_above_zero+1) = XMU_PL_L2(IE)
         XMU_MN_L2(IE-IE_first_above_zero+1) = XMU_MN_L2(IE)
         XMU_Z0_L2(IE-IE_first_above_zero+1) = XMU_Z0_L2(IE)
      END DO
C
c Writes out all the FUNCTIONs which will be integrated - just for debugging:
C
      if( l23edge ) then
C L3 edge:
         OPEN( raw_file_L3_IO, file='raw_file_L3.dat')
C Header:
         WRITE(raw_file_L3_IO,'("# Energy Re(z) from EF, ",
     &   " XMU_PL_L3, XMU_MN_L3, XMU_Z0_L3, DIF_L3, SUM_L3")')
      elseif (m45edge ) then
C M5 edge:
         OPEN( raw_file_L3_IO, file='raw_file_M5.dat')
C Header:
         WRITE(raw_file_L3_IO,'("# Energy Re(z) from EF, ",
     &   " XMU_PL_M5, XMU_MN_M5, XMU_Z0_M5, DIF_M5, SUM_M5")')
      else
         OPEN( raw_file_L3_IO, file='raw_file_K_a.dat')
      endif
C
      WRITE(raw_file_L3_IO,'(f9.5,3x,3f12.6,3x,2f12.6)') (RZERO,i=1,6)
      DO IE = 1, NE_used
         WRITE(raw_file_L3_IO,'(f9.5,3x,3f12.6,3x,2f12.6)') 
     &      EMESH_used(IE), 
     &      XMU_PL_L3(IE), XMU_MN_L3(IE), XMU_Z0_L3(IE),
     &      DIF_L3(IE), SUM_L3(IE)
      ENDDO
      CLOSE(raw_file_L3_IO)
C
C
      if( l23edge ) then
C L2 edge:
         OPEN( raw_file_L2_IO, file='raw_file_L2.dat')
C Header:
         WRITE(raw_file_L2_IO,'("# Energy Re(z) from EF, ",
     &   " XMU_PL_L2, XMU_MN_L2, XMU_Z0_L2, DIF_L2, SUM_L2")')
      elseif (m45edge ) then
C M4 edge:
         OPEN( raw_file_L2_IO, file='raw_file_M4.dat')
C Header:
         WRITE(raw_file_L2_IO,'("# Energy Re(z) from EF, ",
     &   " XMU_PL_M4, XMU_MN_M4, XMU_Z0_M4, DIF_M4, SUM_M4")')
      else
         OPEN( raw_file_L2_IO, file='raw_file_K_b.dat')
      endif
C
C
      WRITE(raw_file_L2_IO,'(f9.5,3x,3f12.6,3x,2f12.6)') (RZERO,i=1,6)
      DO IE = 1, NE_used
         WRITE(raw_file_L2_IO,'(f9.5,3x,3f12.6,3x,2f12.6)') 
     &      EMESH_used(IE), 
     &      XMU_PL_L2(IE), XMU_MN_L2(IE), XMU_Z0_L2(IE),
     &      DIF_L2(IE), SUM_L2(IE)
      ENDDO
      CLOSE(raw_file_L2_IO)
C
c Integrates the spectra from zero to E_band
C
      if( l23edge ) then
         OPEN( raw_file_L3_IO, file='integrand_DIF_L3.dat')
      elseif (m45edge ) then
         OPEN( raw_file_L3_IO, file='integrand_DIF_M5.dat')
      else
         OPEN( raw_file_L3_IO, file='integrand_DIF_K_a.dat')
      endif
!       WRITE(666,'("# Integrand DIF_L3 over NE_used=",i0)') 
!     &   NE_used
      WRITE(raw_file_L3_IO,'(F15.8," ",F15.8)') (EMESH_used(IE),
     &   DIF_L3(IE),IE=1,NE_used)
      CLOSE(raw_file_L3_IO)
C
      if( l23edge ) then
         OPEN( raw_file_L2_IO, file='integrand_DIF_L2.dat')
      elseif (m45edge ) then
         OPEN( raw_file_L2_IO, file='integrand_DIF_M4.dat')
      else
         OPEN( raw_file_L2_IO, file='integrand_DIF_K_b.dat')
      endif
!       WRITE(raw_file_L2_IO,'("# Integrand DIF_L2 over NE_used=",i0)') 
!      &   NE_used
      WRITE(raw_file_L2_IO,'(F15.8," ",F15.8)') (EMESH_used(IE),
     &   DIF_L2(IE),IE=1,NE_used)
      CLOSE(raw_file_L2_IO)
C
      DELA3 = trapez( RZERO, RZERO, NE_used, EMESH_used, DIF_L3 )
      DELA2 = trapez( RZERO, RZERO, NE_used, EMESH_used, DIF_L2 )
C
      att = trapez( RZERO, RZERO, NE_used, EMESH_used, SUM_L3 )
     &      + 
     &      trapez( RZERO, RZERO, NE_used, EMESH_used, SUM_L2 )
C
C Compute the ratios
C
      spin = 3.0_DP * ( DELA3 - 2.0_DP * DELA2 ) / att
      orb = 2.0_DP * ( DELA3 + DELA2 ) / att
      ratio = 2.0_DP/3.0_DP * ( DELA3 + DELA2 ) / 
     &   ( DELA3 - 2.0_DP*DELA2 )
C
c Inverts the signs - I want the moments to be positive
C
      spin = -spin
      orb = -orb
C
c Writes out the ratios
C
      WRITE(6,'(a,a,f9.4,/,a,a,I1,a,f4.2,/,5x,a,/,5x,3f13.4)')
     &     'XMCD sum rules from theoretical spectra:',
     &     '    E_band =', EBAND,
     &     TXTT(1)(1:LTXTT(1)), '   Nsites=', NQT(1),
     &     '     Conc=', CONC(1),
     &     '     (S+7*T)/n       O/n       O/(S+7*T)',
     &     spin, orb, ratio
C
      END PROGRAM sum_rules_theory
C
