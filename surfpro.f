c========================================================================
c-------------- THE ROUTINE FOR PROJECTING ON A SURFACE ---------------
c========================================================================

      subroutine surfprojection(intpt,alpha,beta)  ! This routine to modify mesh coordinates
      include 'SIZE'
      include 'TOTAL'
      integer intpt
      real    alpha,
     &        beta

      integer i_closest, 
     &        iter,
     &        i,j,k,
     &        e,ie,f,iside, ! Element,local element pre face and lex face 
     &        ix,iy,iz,  ! indices
     &        kx1,kx2,ky1,ky2,kz1,kz2,
     &        iinterp(3) ! Interpolation indices 
      real rinterp(3), ! Interpolation distances
     &     midpoint(3), wrk(3),sn(3),
     &     lambda
     &     err_surf,err_new

      common /surfstat/ wrksurf(nsurf,3),iwrkelem(nbdry,3)
      real wrkbdry(nwork,3)

c ---- Function 1 starts here ---- create_working_surface---
       call rzero(wrksurf(1,1),3*nsurf)
       call rzero(wrkbdry(1,1),3*nwork)
       call izero(iwrkelem(1,1),3*nbdry)

          !! Reading the list of surface elements
          open(unit=13,file='surf.i',status='OLD')
          do i = 1,nsurf
              read(13,*) wrksurf(i,1),wrksurf(i,2),wrksurf(i,3)
          enddo
          close(unit=13)

          !! Reading the list of boundary elements
          open (unit=12,file='bdry.i',status='OLD')
          do i = 1,nbdry
            read(12,*) iwrkelem(i,1),iwrkelem(i,2)
          enddo
          close(unit=12)

          ! Iterating over the Boundary elements
          do i = 1,nbdry
          e = iwrkelem(i,1) ! Global element 
          ie = gllel(e)     ! Local element

          !Only perform this procedure if the element
          !Corresponds to the processor

          if(lglel(ie).eq.e) then 

          f = iwrkelem(i,2) ! preprocessor face
          iside = eface(f) ! symmetric notation face

          call rzero(wrkbdry,3*nwork)
          call rzero(midpoint,3)
          call rzero(wrk,3)

          !Getting the iteration limits
          call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)

          ! Calculating the radius of interest
          call nekasgn(kx1,ky1,kz1,ie)
          wrk(1) = x
          wrk(2) = y
          wrk(3) = z
          call nekasgn(kx2,ky2,kz2,ie)

          ! This calculation could be bad if the face is very 
          ! distorted... 

          radius = 1.0/2*sqrt((wrk(1)-x)**2+(wrk(2)-y)**2+(wrk(3)-z)**2)

          ! Calculating the midpoint on the face
          midpoint(1) = (wrk(1)-x)/2 + x
          midpoint(2) = (wrk(2)-y)/2 + y
          midpoint(3) = (wrk(3)-z)/2 + z

          ! Finding the points of interest
          iter = 1
          do j = 1,nsurf
              if (sqrt((wrksurf(j,1)-midpoint(1))**2+
     $        (wrksurf(j,2)-midpoint(2))**2+
     $        (wrksurf(j,3)-midpoint(3))**2).le.beta*radius) then 
                  wrkbdry(iter,1) =wrksurf(j,1)
                  wrkbdry(iter,2) =wrksurf(j,2)
                  wrkbdry(iter,3) =wrksurf(j,3)
                  iter = iter+1              
              end if
              if(iter.eq.nwork) 
     $        write(*,*) 'ERR:Increase wrkbdry in SIZE or reduce radius'
          enddo
c--------- First function ends here !! 

          ! Iterating over the face
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
c---------- This corresponds to the initialization of the interp array
          ! Get coordinates
          call nekasgn(ix,iy,iz,ie) 

          ! Get surface normal
          call getsurfnorm(sn,ix,iy,iz,f,ie) 

          ! Iterating through wrkbdry
          i_closest = 1             ! integer holding the optimal surf point
          err_surf = 999999.9       ! min error on surface
          call cfill(rinterp,999999.9,3)
          call izero(iinterp,3)

c---------- End initialization of the interp array
          ! Finding the interpolation points
              do k = 1,iter-1  
                wrk(1) = x-wrkbdry(k,1) ! vector from gllpoint to surf 
                wrk(2) = y-wrkbdry(k,2)
                wrk(3) = z-wrkbdry(k,3)
                lambda = sqrt(wrk(1)**2+wrk(2)**2+wrk(3)**2)
                if(lambda.eq.0.0) then
                    i_closest = k
                    err_surf = 0.0
                    continue 
                end if
                call cmult(wrk,1.0/lambda,3) ! Normalizing
                
                ! flip vector if necessary 
                if((wrk(1)*sn(1)+wrk(2)*sn(2)+wrk(3)*sn(3)).lt.0.0) 
     &            call chsign(wrk,3) 

                ! Calculating error 
                call calcerror(err_new,lambda,sn,wrk,radius,alpha)

                ! Updating the interpolation points
                call interp_up(iinterp,rinterp,intpt,err_new,k)
              enddo ! k
                !set_new_pt(wrkbdry,intp,rinterp,n,ix,iy,iz,ie)
    
          enddo ! ix 
          enddo ! iy 
          enddo ! iz 

          ! Redistributing the internal nodes
          call fix_gll(ie,f) 

          end if! if global node corresponds to this processor 
          enddo ! i
          ! Fixing the internal GLL-points
          param(59) = 1.0
          call fix_geom()
c============== END ROUTINE FOR PROJECTING ON A SURFACE ===============
          
          return
          end
C=======================================================================
      subroutine fix_gll(e,f) ! Need to redistribute the nodes above 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
c
      ! Redistribute the gll-points between the given face and the
      ! opposite to make sure that all points lie within the element.
c
      integer   f_opp(6) !element and face
      save      f_opp
      !data      f_opp / 4,3,2,1,6,5 /  ! opposite face
      data      f_opp / 3,4,1,2,6,5 /   ! opposite face
      integer iter,inneriter,
     &        ix,iy,iz,nnx,nny,nnz,     ! iteration indices
     &        e,f,                      ! Element and pre face
     &        kx1,kx2,ky1,ky2,kz1,kz2,  ! it limits face1
     &        jx1,jx2,jy1,jy2,jz1,jz2   ! it limiits face2
      real zg(nx1),wg(nx1)              ! GLL points and weights
      real wrk1(nx1*nx1,3),             ! working arrays
     &     wrk2(nx1*nx1,3) 

      call zwgll(zg,wg,nx1)  ! Get the GLL-points
      call cadd(zg,1.0,nx1)  ! truncating
      call cmult(zg,0.5,nx1) !normalizing
      !do iter = 1,nx1        ! Just an equidistant distribution
        !zg(iter) = (iter-1)*1.0/(nx1-1)
      !enddo
      ! Assigning the face indices
      call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
      call facind(jx1,jx2,jy1,jy2,jz1,jz2,nx1,ny1,nz1,f_opp(f))
      ! Get the first face (Filling wrk1)
      call getface(kx1,kx2,ky1,ky2,kz1,kz2,wrk1,nx1,e)
      ! Get the second face (Filling wrk2)
      call getface(jx1,jx2,jy1,jy2,jz1,jz2,wrk2,nx1,e)
      ! Get the vector wrk2 = wrk2 - wrk1
      call sub2(wrk2,wrk1,nx1*nx1*3)
      ! Redistributing 
      iter = 1
      do iz=kz1,kz2 ! iterating over the face
      do iy=ky1,ky2
      do ix=kx1,kx2
        nnx = ix ! Indices for the line between the faces
        nny = iy
        nnz = iz
        if(kx1.ne.jx1) nnx = jx1  
        if(ky1.ne.jy1) nny = jy1  
        if(kz1.ne.jz1) nnz = jz1  
          inneriter = 1
          do jz=iz,nnz ! iterating between the faces
          do jy=iy,nny
          do jx=ix,nnx
      xm1(jx,jy,jz,e) =xm1(ix,iy,iz,e)+zg(inneriter)*wrk2(iter,1)
      ym1(jx,jy,jz,e) =ym1(ix,iy,iz,e)+zg(inneriter)*wrk2(iter,2)
      zm1(jx,jy,jz,e) =zm1(ix,iy,iz,e)+zg(inneriter)*wrk2(iter,3)
          inneriter = inneriter+1
          enddo
          enddo
          enddo
      iter = iter+1
      enddo
      enddo
      enddo

       return 
       end
c---------------------------------------------------------------
      subroutine getface(kx1,kx2,ky1,ky2,kz1,kz2,wrk,n,e)
c     OUTPUT:
c     assigning the values of the face in element e corresponding 
c     to the index limits kx1,kx2... to the array wrk.
c     INPUT:
c     kx1,kx2... are index limits for the face of interest
c     n is the number of GLL-points in each direction
c     e is the current global element
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer n,e,iter
     &        ix,iy,iz,
     &        kx1,kx2,ky1,ky2,kz1,kz2 ! it limits face
      real wrk(n*n,3) ! working array

      iter = 1
      do iz=kz1,kz2 ! iterating over the face
      do iy=ky1,ky2
      do ix=kx1,kx2
          !call nekasgn(ix,iy,iz,e)
          wrk(iter,1) = xm1(ix,iy,iz,e)
          wrk(iter,2) = ym1(ix,iy,iz,e)
          wrk(iter,3) = zm1(ix,iy,iz,e)
          iter = iter+1
      enddo
      enddo
      enddo

      return
      end
c---------------------------------------------------------------
      subroutine getsurfnorm(sn,ix,iy,iz,f,ie)
c     Providing the surface normal sn at point ix,iy,iz of element ie
c     and face f 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
          real sn(3)
          integer ix,iy,iz
          integer ie,f,iside

          iside = eface(f) ! Providing the lexicographic notation
          ! Get the normal vector
          if (f.eq.1.or.f.eq.2) then ! "r face"
              sn(1) = unx(iy,iz,iside,ie)
              sn(2) = uny(iy,iz,iside,ie)
              sn(3) = unz(iy,iz,iside,ie)
          elseif (f.eq.3.or.f.eq.4) then ! "s face"
              sn(1) = unx(ix,iz,iside,ie)
              sn(2) = uny(ix,iz,iside,ie)
              sn(3) = unz(ix,iz,iside,ie)
          elseif (f.eq.5.or.f.eq.6) then ! "t face"
              sn(1) = unx(ix,iy,iside,ie)
              sn(2) = uny(ix,iy,iside,ie)
              sn(3) = unz(ix,iy,iside,ie)
          end if
          return 
          end
c===============================================================
         subroutine calcerror(error,lambda,sn,wrk,radius,alpha)
           real error,  ! the error to be calculated
     &          lambda, ! the absolute distance from x_gll to x_surf
     &          sn(3),  ! sn is the element face normal
     &          wrk(3), ! wrk is the normalized vector between x_gll and x_surf
     &          radius, ! lengthscale of the current element
     &          alpha   ! Weight coefficient 0<alpha<1
             
             ! Many options (Should be some experimentation)
             ! No good effect of these two options ... 
             error = (1.0-alpha)*sqrt((wrk(1)-sn(1))**2+ ! angular diff
     &                    (wrk(2)-sn(2))**2+
     &                    (wrk(3)-sn(3))**2)

             error = alpha*lambda/radius+error ! an inbetween soln

             !error = lambda ! std L2 norm
         return
         end
c===============================================================
        subroutine interp_up(iinterp,rinterp,n,error,k)
        real error,     ! new error
     &       rinterp(n) ! current errors
        integer k,      ! current surf point
     &          n,      ! Iteration pts
     &          i,      ! iteration index
     &          ipos,   ! temp
     &          iinterp(n) !current interp. points
        ipos = 0
        do i = 1,n
        if(error.lt.rinterp(i)) ipos = i 
        enddo
        if ( ipos .gt. 0 ) then !
          do i=1,ipos-1,1
           iinterp(i) = iinterp(i+1)
           rinterp(i)= rinterp(i+1)
          enddo 
          iinterp(ipos) = k
          rinterp(ipos) = error
      end if

      return
      end
c===============================================================
      subroutine set_new_pt(wrkbdry,iinterp,rinterp,n,ix,iy,iz,e)
c     OUTPUT:
c     defining the position of the new gll-point on the surface
c     INPUT:
c     iinter contains the position of the three closest points 
c     rinter contains the distance to the three closest points 
c     ix,iy,iz are the indexes of the gll point in element e
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer i,n,ix,iy,iz,e,
     &        iinterp(n)
      real rinterp(n),d,xx,yy,zz
      !common /surfproj/ wrkbdry(nwork,3)
      real wrkbdry(nwork,3)
      integer itemp
      d = 0.0
      xx = 0
      yy = 0 
      zz = 0
      do i = 1,n
        itemp = iinterp(i) 
        rtemp = rinterp(i)
        d  = d  +              1.0/rtemp**2 
        xx = xx + wrkbdry(itemp,1)/rtemp**2
        yy = yy + wrkbdry(itemp,2)/rtemp**2
        zz = zz + wrkbdry(itemp,3)/rtemp**2
        if (rtemp.eq.0.0) then 
            xx = wrkbdry(itemp,1)
            yy = wrkbdry(itemp,2)
            zz = wrkbdry(itemp,3)
            d = 1.0
            continue
       end if
      enddo
      ! Redistributing the points
      xm1(ix,iy,iz,e) = xx/d
      ym1(ix,iy,iz,e) = yy/d
      zm1(ix,iy,iz,e) = zz/d
      !if(nio.eq.0) write(*,*) sqrt((x-xx/d)**2+(y-yy/d)**2+(z-zz/d)**2)
      return 
      end

c-----------------------------------------------------------------------
