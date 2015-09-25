      integer inNY, inNZ
      parameter (inNY=3*(348-1), inNZ=450)
      integer nfile0, nfile1
      parameter (nfile0=0, nfile1=18564)
      real dtinf
      parameter (dtinf=0.0016)
      integer nintpt
      parameter (nintpt=4)
      real rintpw
      parameter (rintpw=1.0)

      real xi(inny,innz), yi(inny,innz), zi(inny,innz)
      common /infcrd/ xi, yi, zi
      real vxi(inny,innz), vyi(inny,innz), vzi(inny,innz)
      common /infvel/ vxi, vyi, vzi

      integer maxlist, maxloc
      parameter (maxlist=20000, maxloc=2000)
      integer iloc(2,maxloc)
      integer intj(nintpt,maxlist), intk(nintpt,maxlist)
      real wgtint(nintpt,maxlist)
      common /interpol/ iloc, intj, intk, wgtint
