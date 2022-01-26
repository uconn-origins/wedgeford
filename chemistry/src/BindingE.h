C..............................................................................
C
C..............................................................................
      INTEGER nBEspec
      PARAMETER (nBEspec=900)
      REAL*8 bindenerg
      CHARACTER*13 bindspec
      DIMENSION bindenerg(nBEspec),bindspec(nBEspec) 
      COMMON /BNDE/ bindenerg,bindspec