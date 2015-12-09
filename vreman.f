c-----------------------------------------------------------------------
c
c     Based on the turbChannel example
c
c-----------------------------------------------------------------------
       subroutine uservp (ix,iy,iz,eg)
       include 'SIZE'
       include 'TOTAL'
       include 'NEKUSE'

       common /cdsmag/ ediff(lx1,ly1,lz1,lelv)
       integer ie,f,eg
       real aux
       ie = gllel(eg)

       udiff = ediff(ix,iy,iz,ie)
c      udiff = param(2)
       utrans=1.

       aux = ediff(ix,iy,iz,ie) - param(2)
       if (istep.eq.19) then
         write(6,*) aux
       endif
       return
       end
c-----------------------------------------------------------------------
       subroutine userf  (ix,iy,iz,eg)
       include 'SIZE'
       include 'TOTAL'
       include 'NEKUSE'

       integer e,f,eg
c     e = gllel(eg)

       ffx = 0.0
       ffy = 0.0
       ffz = 0.0

       return
       end
c-----------------------------------------------------------------------
       subroutine userq  (ix,iy,iz,ieg)
       include 'SIZE'
       include 'TOTAL'
       include 'NEKUSE'

       qvol = 0.

       return
       end
c-----------------------------------------------------------------------
       subroutine userchk
       include 'SIZE'
       include 'TOTAL'
       include 'RESTART'
       include 'ZPER'  ! for nelx,nely,nelz

       parameter(lt=lx1*ly1*lz1*lelv)
       real vort(lt,3),vort_pol(lt,3),w1(lt),w2(lt),rayon

       real x0(3)
       save x0

       integer icalld
       save    icalld
       data    icalld /0/

       character*1   snam1(80)
       character*1   f1nam1(80),f2nam1(80)
       character*80  f1name
       equivalence  (f1nam1,f1name)
       character*80  f2name
       equivalence  (f2nam1,f2name)

       real atime,timel
       save atime,timel

       integer icount
       save    icount

       common /cdsmag/ ediff(lx1,ly1,lz1,lelv)

       common /avg/    uravg(lx1,ly1,lz1,lelv)
      &             ,  utavg(lx1,ly1,lz1,lelv)
      &             ,  uzavg(lx1,ly1,lz1,lelv)

       common /lesleo/ sij (lx1*ly1*lz1,ldim,ldim)
      $              , dx2 (lx1,ly1,lz1,lelv,lp)
      $              , dy2 (lx1,ly1,lz1,lelv,lp)
      $              , dz2 (lx1,ly1,lz1,lelv,lp)
      $              , cdiff
       real sij,dx2,dy2,dz2,cdiff

       integer e
       logical ifverbose

       n=nx1*ny1*nz1*nelv

       pi    = 4.*atan(1.0)
       rho   = 1.
       dnu   = param(2)

       cdiff = 0.07

       ! add SGS term to RHS (explicit treatment in plan4)
       ! Explicit viscosity.
       if(ifsplit .and. param(30).gt.0) ifexplvis = .true.

c     Grid spacing.
c      if ((istep.lt.1).and.(nid.eq.0)) then
c      if (istep.lt.1) then
       if (nid.eq.(istep-1)) then
         call set_grid_spacing(dx2,dy2,dz2)
       endif

       ! Compute eddy viscosity using Vreman model
c     Waiting for several timesteps in order to be sure that gradients
c     are not zero
       if (istep.gt.17) then
           ifuservp = .true.
       else
           ifuservp = .false.
       endif
       if(istep.gt.16) then
         if(nid.eq.0) write(6,*) 'Calculating eddy visosity'
         do e=1,nelv
            call eddy_visc(ediff,e)
         enddo
c        call copy(t,ediff,n)
       endif

c
c     Below is just for postprocessing ...
c
c     Compute perturbation energy
       e2 = glsc3(vy,bm1,vy,n)+glsc3(vz,bm1,vz,n)
       e2 = e2/volvm1

       ifxyo = .true.                        ! Turn on xyz output
       if (istep.gt.iostep) ifxyo = .false.  ! Turn off xyz output after 
first dump

c     Average calculations.
       if(icalld.eq.0) then
         call rzero(uravg,n)
         call rzero(utavg,n)
         call rzero(uzavg,n)
         atime = 0.
         timel = time
         icalld = 1
       endif

       dtime = time - timel
       atime = atime + dtime

       if (atime.ne.0. .and. dtime.ne.0.) then
         beta      = dtime/atime
         alpha     = 1.-beta
         ifverbose = .false.

         call avg1(uravg,vx   ,alpha,beta,n,'uravg',ifverbose)
         call avg1(utavg,vy   ,alpha,beta,n,'utavg',ifverbose)
         call avg1(uzavg,vz   ,alpha,beta,n,'uzavg',ifverbose)
       endif

       timel = time



       ntot = nx1*ny1*nz1*nelv

c     Vorticity calculations.

       if (mod(istep,iostep).eq.0) then
           call comp_vort3(vort,w1,w2,vx,vy,vz)
c         Calculations done in polar coordinates.
           do i=1,ntot
               rayon = SQRT(xm1(i,1,1,1)**2 + ym1(i,1,1,1)**2)
               if (rayon>0.0) then
                   vort_pol(i,1)=(xm1(i,1,1,1)*vort(i,1) +
      $                ym1(i,1,1,1)*vort(i,2))/rayon
                   vort_pol(i,2)=(-ym1(i,1,1,1)*vort(i,1) +
      $                xm1(i,1,1,1)*vort(i,2))/rayon
                   vort_pol(i,3)=vort(i,3)
               else
                   vort_pol(i,1)=0.0
                   vort_pol(i,2)=0.0
                   vort_pol(i,3)=vort(i,3)
               endif
           enddo
c          call outpost(vort(1,1),vort(1,2),vort(1,3),
c     $         pr,t(1,1,1,1,1),'vrt')
           call outpost(vort_pol(1,1),vort_pol(1,2),vort_pol(1,3),
      $         pr,t(1,1,1,1,1),'pol')
       endif

       if (istep.eq.1) then

       endif

c     We export the grid spacing and the average values of velocities.
       if(istep.eq.18) then
           ifxyo = .true.
           call outpost(dx2(1,1,1,1,1),dy2(1,1,1,1,1),dz2(1,1,1,1,1),
      $         ediff(1,1,1,1),t(1,1,1,1,1),'les')
           call outpost(uravg(1,1,1,1),utavg(1,1,1,1),uzavg(1,1,1,1),
      $         pr,t(1,1,1,1,1),'avg')
           if (nid.eq.0) then
           open(223,file='diffusion.dat',status='old',action='write',
      $      form='formatted',position="append")
           do i=1,ntot
c              write(223,*) xm1(i,1,1,1),ym1(i,1,1,1),zm1(i,1,1,1),nid,
c     $ ediff(i,1,1,1),dx2(i,1,1,1,nid+1),dy2(i,1,1,1,nid+1),
c     $            dz2(i,1,1,1,nid+1)
               write(223,*) nid
           enddo
           close(223)
           endif
           if (nid.eq.1) then
open(225,file='diffusion1.dat',status='old',action='write',
      $      form='formatted',position="append")
           do i=1,ntot
c              write(223,*) xm1(i,1,1,1),ym1(i,1,1,1),zm1(i,1,1,1),nid,
c     $ ediff(i,1,1,1),dx2(i,1,1,1,nid+1),dy2(i,1,1,1,nid+1),
c     $            dz2(i,1,1,1,nid+1)
               write(225,*) nid
           enddo
           close(225)
           endif

       endif
       return
       end
c-----------------------------------------------------------------------
       subroutine userbc (ix,iy,iz,iside,ieg)
c     NOTE ::: This subroutine MAY NOT be called by every process
       include 'SIZE'
       include 'TOTAL'
       include 'NEKUSE'
c
       ux=0.0
       uy=0.0
       uz=0.0
       temp=0
c
       z0 = 0.05
       rayon = SQRT(x**2 + y**2) ! Radius in polar coordinates
       R = 1.0                   ! Radius of the jet
       theta0 = 0.05*R           ! Thickness of the mixing layer
       axialamp = 0.05        ! Turbulence amplitude (axial excitation)
       helicamp = 0.05        ! Turbulence amplitude (helical direction)
       faxial = 0.55             ! Axial frequency = St*U0/D et D=1.0
       fhelic = faxial/2.0       ! Helical frequency = faxial/2.0
       phihelic = 0.0            ! Phase lag between axial and helical
c
       if (z.le.z0) then
           ux=0.0
           uy=0.0
           if (rayon.gt.(0.0)) then
               uz=0.5*(1-tanh(0.25*(rayon/R - R/rayon)*R/theta0))
           else
               uz=1.0
           endif
       endif
c
       return
       end
c-----------------------------------------------------------------------
       subroutine useric (ix,iy,iz,ieg)
       include 'SIZE'
       include 'TOTAL'
       include 'NEKUSE'

       rayon = SQRT(x**2 + y**2)  ! Radius in polar coordinates
       R = 0.5                    ! Radius of the jet
       theta0 = 0.05*R            ! Thickness of the mixing layer

       ux=0.0
       uy=0.0
       if (rayon.gt.(0.0)) then
           uz=0.5*(1-tanh(0.25*(rayon/R - R/rayon)*R/theta0))
       else
           uz=1.0
       endif
       temp=0
       return
       end
c-----------------------------------------------------------------------
       subroutine usrdat
       include 'SIZE'
       include 'TOTAL'
       common /cdsmag/ ediff(lx1,ly1,lz1,lelv)

       n = nx1*ny1*nz1*nelt
       call cfill(ediff,param(2),n)  ! initialize viscosity

       ! enable stress formulation if we use PN-PN-2
       if(lx1.ne.lx2 .and. param(30).gt.0) IFSTRS = .true.

       return
       end
c-----------------------------------------------------------------------
       subroutine usrdat2
       include 'SIZE'
       include 'TOTAL'


       return
       end
c-----------------------------------------------------------------------
       subroutine usrdat3
       include 'SIZE'
       include 'TOTAL'
c
       return
       end
c-----------------------------------------------------------------------
       subroutine set_obj  ! define objects for surface integrals
c
       include 'SIZE'
       include 'TOTAL'

       integer e,f,eg

       nobj = 1
       iobj = 0
       do ii=nhis+1,nhis+nobj
          iobj = iobj+1
          hcode(10,ii) = 'I'
          hcode( 1,ii) = 'F'
          hcode( 2,ii) = 'F'
          hcode( 3,ii) = 'F'
          lochis(1,ii) = iobj
       enddo
       nhis = nhis + nobj

       if (maxobj.lt.nobj) call exitti('increase maxobj in SIZE$',nobj)
       nxyz  = nx1*ny1*nz1
       nface = 2*ndim

       do e=1,nelv
       do f=1,nface
          if (cbc(f,e,1).eq.'W  ') then
             iobj  = 1
             if (iobj.gt.0) then
                nmember(iobj) = nmember(iobj) + 1
                mem = nmember(iobj)
                eg  = lglel(e)
                object(iobj,mem,1) = eg
                object(iobj,mem,2) = f
c              write(6,1) iobj,mem,f,eg,e,nid,' OBJ'
c   1          format(6i9,a4)

             endif
          endif
       enddo
       enddo
c     write(6,*) 'number',(nmember(k),k=1,4)
c
       return
       end
c-----------------------------------------------------------------------
       subroutine eddy_visc(ediff,e)
c
c     Compute eddy viscosity using Vreman's model.
c
       include 'SIZE'
       include 'TOTAL'
       include 'ZPER'

       real ediff(nx1*ny1*nz1,nelv)
       real betaij(lx1*ly1*lz1,ldim,ldim)
       real Bbeta(lx1*ly1*lz1)
       real Sommeij2(lx1*ly1*lz1)
       integer e,i,ldim1,ldim2

       common /lesleo/ sij (lx1*ly1*lz1,ldim,ldim)
      $              , dx2 (lx1*ly1*lz1,lelv,lp)
      $              , dy2 (lx1*ly1*lz1,lelv,lp)
      $              , dz2 (lx1*ly1*lz1,lelv,lp)
      $              , cdiff
       real sij,dx2,dy2,dz2,cdiff
       integer ntot

       ntot = nx1*ny1*nz1

c     Velocity gradient.
       call comp_gije(sij,vx(1,1,1,e),vy(1,1,1,e),vz(1,1,1,e),e)

c     Calculation of Beta.
       do ldim2=1,ldim
         do ldim1=1,ldim
           do i=1,ntot
             betaij(i,ldim1,ldim2) =
      $        dx2(i,e,nid+1)*sij(i,1,ldim1)*sij(i,1,ldim2)
      $        + dy2(i,e,nid+1)*sij(i,2,ldim1)*sij(i,2,ldim2)
      $        + dz2(i,e,nid+1)*sij(i,3,ldim1)*sij(i,3,ldim2)
           enddo
         enddo
       enddo


c     Calculation of Bbeta.
       do i=1,ntot
         Bbeta(i) = betaij(i,1,1)*betaij(i,2,2) - (betaij(i,1,2)**2)
      $           + betaij(i,2,2)*betaij(i,3,3) - (betaij(i,2,3)**2)
      $           + betaij(i,3,3)*betaij(i,1,1) - (betaij(i,3,1)**2)
       enddo

c     Calculons Sommeij2 = Sum_{i,j} (alpha_{ij}^2) pour le modèle de
c     Vreman. On met d'abord Sommeij2 à zero.
       call rzero(Sommeij2,ntot)
       do ldim2=1,ldim
         do ldim1=1,ldim
           do i=1,ntot
             Sommeij2(i) = Sommeij2(i) + (sij(i,ldim1,ldim2)**2)
           enddo
         enddo
       enddo

c     ediff calculation.
       do i=1,ntot
         if (Sommeij2(i).gt.0.0) then
           ediff(i,e) = param(2) + cdiff*sqrt(Bbeta(i)/Sommeij2(i))
c          ediff(i,e) = param(2)
         else
           ediff(i,e) = param(2)
         endif
       enddo

       return
       end
c-----------------------------------------------------------------------
       subroutine set_grid_spacing(dx2,dy2,dz2)
c
c     Compute D^2, the grid spacing used in the DS sgs model.
c
       include 'SIZE'
       include 'TOTAL'


       real dx2(nx1,ny1,nz1,nelv,np)
       real dy2(nx1,ny1,nz1,nelv,np)
       real dz2(nx1,ny1,nz1,nelv,np)
       real longueur,largeur,dhlong,dhlarg,deltai,deltaj

       integer e,i,j,k,im,jm,km,ip,jp,kp

c      return               ! Comment this line for a non-trivial Delta defn

       write(6,*) 'LadHyX'

       open (224,file='dx2dy2dz2.dat',status='old',action='write',
      $      form='formatted',position="append")
       do e=1,nelv
          do k=1,nz1
            km = max(1  ,k-1)
            kp = min(nz1,k+1)
            do j=1,ny1
              jm = max(1  ,j-1)
              jp = min(ny1,j+1)
              do i=1,nx1
                im = max(1  ,i-1)
                ip = min(nx1,i+1)
                longueur = 0.5*(abs(xm1(ip,jp,k,e) - xm1(im,jp,k,e))
      $           + abs(xm1(ip,jm,k,e) - xm1(im,jm,k,e)))/(ip - im)
                largeur = 0.5*(abs(ym1(ip,jp,k,e) - ym1(ip,jm,k,e))
      $           + abs(ym1(im,jp,k,e) - ym1(im,jm,k,e)))/(jp - jm)
                dhlong = 0.5*(abs(ym1(ip,jp,k,e) - ym1(im,jp,k,e))
      $           + abs(ym1(ip,jm,k,e) - ym1(im,jm,k,e)))/(ip - im)
                dhlarg = 0.5*(abs(xm1(ip,jp,k,e) - xm1(ip,jm,k,e))
      $           + abs(xm1(im,jp,k,e) - xm1(im,jm,k,e)))/(jp - jm)
                deltai = sqrt(longueur**2 + dhlong**2)
                deltaj = sqrt(largeur**2 + dhlarg**2)

                dx2(i,j,k,e,nid+1) = ((longueur**2)/deltai + (dhlarg**2)/
      $                              deltaj)**2
                dy2(i,j,k,e,nid+1) = ((dhlong**2)/deltai + (largeur**2)/
      $                              deltaj)**2
                dz2(i,j,k,e,nid+1) = ((zm1(i,j,kp,e)-zm1(i,j,km,e))/
      $                              (kp-km))**2
                write(224,*) i,j,k,e,nid,dx2(i,j,k,e,nid+1),
      $                      dy2(i,j,k,e,nid+1),dz2(i,j,k,e,nid+1)
              enddo
            enddo
          enddo
       enddo
       close(224)

       return
       end
c-----------------------------------------------------------------------
c
c automatically added by makenek
       subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
       integer*8 glo_num(1)
       return
       end
