C Common blocks for 'dvodpk.f' ODEs solver:
C
C Type declarations for labeled COMMON block DVOD01 --------------------
C
      DOUBLE PRECISION ACNRM, CCMXJ, CONP, CRATE, DRC, EL,
     1     ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2     RC, RL1, TAU, TQ, TN, UROUND
      INTEGER ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     1        L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     2        LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     3        N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     4        NSLP, NYH
C
C Type declarations for labeled COMMON block DVOD02 --------------------
C
      DOUBLE PRECISION HU
      INTEGER NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
C
C Type declarations for labeled COMMON block DVPK01 --------------------
C
      DOUBLE PRECISION DELT, SQRTN, RSQRTN
      INTEGER JPRE, JACFLG, LOCWP, LOCIWP, LVSAV, KMP, MAXL, MNEWT,
     1      NLI, NPS, NCFL
C
      COMMON /DVOD01/ ACNRM, CCMXJ, CONP, CRATE, DRC, EL(13),
     1                ETA, ETAMAX, H, HMIN, HMXI, HNEW, HSCAL, PRL1,
     2                RC, RL1, TAU(13), TQ(5), TN, UROUND,
     3                ICF, INIT, IPUP, JCUR, JSTART, JSV, KFLAG, KUTH,
     4                L, LMAX, LYH, LEWT, LACOR, LSAVF, LWM, LIWM,
     5                LOCJS, MAXORD, METH, MITER, MSBJ, MXHNIL, MXSTEP,
     6                N, NEWH, NEWQ, NHNIL, NQ, NQNYH, NQWAIT, NSLJ,
     7                NSLP, NYH
      COMMON /DVOD02/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST
      COMMON /DVPK01/ DELT, SQRTN, RSQRTN, JPRE, JACFLG, LOCIWP,
     1                LOCWP, LVSAV, KMP, MAXL, MNEWT, NLI, NPS, NCFL

C  /BLK2/:
      INTEGER iwork
      REAL*8 rwork
      
      COMMON /BLK2/ rwork(nspec*nss+721+32*nspec),iwork(30+5*nspec)
