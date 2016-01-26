c========================================================================
c-------------- THE ROUTINE FOR PROJECTING ON A SURFACE ---------------
c========================================================================

      subroutine surfprojection(intpt,alpha,beta)  ! This routine to modify mesh coordinates
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
c     OUTPUT:
c     Distributes the GLL-nodes onto the surface prescribed in
c     surf.i.
c     INPUT:
c     intpt is the number of interpolation points to use (Max 3.)
c     alpha is the amount of L2/
c     e is the current global element

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
     &     lambda,mcx(2),
     &     err_surf,err_new

      integer iwrkelem(nbdry,3),iscorner
      real wrksurf(nsurf,3),wrkbdry(nwork,3)

      !nx1 = lx1
      !ny1 = ly1
      !nz1 = lz1
       
       call rzero(mcx,2) ! Max deviation
       call rzero(wrksurf(1,1),3*nsurf)
       call rzero(wrkbdry(1,1),3*nwork)
       call izero(iwrkelem(1,1),3*nbdry)

       if(nid.eq.0) then 
          !! Reading the list of surface elements
          open(unit=13,file='surf.i',status='OLD')
          if (nid.eq.0) write(*,*) 'Opening: surf.i' 
          do i = 1,nsurf
              read(13,*) wrksurf(i,1),wrksurf(i,2),wrksurf(i,3)
          end do
          close(unit=13)

          !! Reading the list of boundary elements
          open (unit=12,file='bdry.i',status='OLD')
          if (nid.eq.0) write(*,*) 'Opening: bdry.i' 
          do i = 1,nbdry
            read(12,*) iwrkelem(i,1),iwrkelem(i,2)
          end do
          close(unit=12)
        endif
        call bcast(i,ISIZE)
        call bcast(wrksurf,WDSIZE*3*nsurf)
        call bcast(iwrkelem,ISIZE*2*nbdry)
        !write(*,*) nid,iwrkelem(787,1),iwrkelem(7,1)

          !do i = 1,nsurf
              !write(*,*) wrksurf(i,1),wrksurf(i,2),wrksurf(i,3)
          !end do

          !do i = 1,nbdry
            !write(*,*) iwrkelem(i,1),iwrkelem(i,2)
          !end do

          ! Iterating over the Boundary elements
          do i = 1,nbdry
          e = iwrkelem(i,1) ! Global element 
          ie = gllel(e)     ! Local element

          !write(*,*) 'elem',e,ie
          !Only perform this procedure if the element
          !Corresponds to the processor

          if(lglel(ie).eq.e) then 

          f = iwrkelem(i,2) ! preprocessor face
          iside = eface(f) ! symmetric notation face
          !write(*,*) 'face and iside',f,iside

          call rzero(wrkbdry,3*nwork)
          call rzero(midpoint,3)
          call rzero(wrk,3)

          !Getting the iteration limits
          call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)

          !call getFace

          ! Calculating the radius of interest
          !call nekasgn(kx1,ky1,kz1,ie)
          ! One corner of the face
          wrk(1) = xm1(kx1,ky1,kz1,ie)
          wrk(2) = ym1(kx1,ky1,kz1,ie)
          wrk(3) = zm1(kx1,ky1,kz1,ie)

          ! opposite corner of the face
          midpoint(1) = xm1(kx2,ky2,kz2,ie)
          midpoint(2) = ym1(kx2,ky2,kz2,ie)
          midpoint(3) = zm1(kx2,ky2,kz2,ie)

          ! This calculation could be bad if the face is very 
          ! distorted... 

          radius = 1.0/2*sqrt(
     $              (wrk(1)-midpoint(1))**2+
     $              (wrk(2)-midpoint(2))**2+
     $              (wrk(3)-midpoint(3))**2)

          !write(*,*) 'radius :',radius
          if(radius.le.0.0) then 
            write(*,*) 'ERROR: radius too small:',radius
            write(*,*) 'LIM:',kx1,ky1,kz1,kx2,ky2,kz2
            write(*,*) 'pt1:', wrk(1), wrk(2), wrk(3)
            write(*,*) 'pt2:', midpoint(1),midpoint(2),midpoint(3)
            write(*,*) 'elem, prepro face:',ie,e, f
            call exitt        
          endif

          ! Calculating the midpoint on the face
          midpoint(1) = (wrk(1)-midpoint(1))/2 + midpoint(1)
          midpoint(2) = (wrk(2)-midpoint(2))/2 + midpoint(2)
          midpoint(3) = (wrk(3)-midpoint(3))/2 + midpoint(3)

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
              if(iter.eq.(nwork-4)) then
              write(*,*) 'ERR:Increase wrkbdry in SIZE or reduce radius'
              call exitt
              end if
          enddo
          ! Add corner nodes 
          !call addcorners(wrkbdry(iter,1),ie,f)

          ! Iterating over the face
          do iz=kz1,kz2
          do iy=ky1,ky2
          do ix=kx1,kx2
          iscorner = 0
          call checkcorner(iscorner,ix,iy,iz,kz1,kz2,ky1,ky2,kx1,kx2)
          if(iscorner.eq.-1) then 
              write(*,*) 'corner:',ix,iy,iz
          else 


          ! Get coordinates
          call nekasgn(ix,iy,iz,ie) 
          ! Get surface normal
          call getsurfnorm(sn,ix,iy,iz,f,ie) 

          ! Iterating through wrkbdry
          i_closest = 1             ! integer holding the optimal surf point
          err_surf = 999999.9       ! min error on surface
          call cfill(rinterp,999999.9,3)
          call izero(iinterp,3)

          ! Finding the interpolation points
              do k = 1,iter-1
                wrk(1) = x-wrkbdry(k,1) ! vector from gllpoint to surf 
                wrk(2) = y-wrkbdry(k,2)
                wrk(3) = z-wrkbdry(k,3)
                lambda = sqrt(wrk(1)**2+wrk(2)**2+wrk(3)**2)
                !if(lambda.eq.0.0) then
                    !i_closest = k
                    !err_surf = 0.0
                    !continue 
                !end if
                call cmult(wrk,1.0/lambda,3) ! Normalizing
                
                ! flip vector if necessary 
                if((wrk(1)*sn(1)+wrk(2)*sn(2)+wrk(3)*sn(3)).lt.0.0) 
     &            call chsign(wrk,3) 

                ! Calculating error 
                call calcerror(err_new,lambda,sn,wrk,radius,alpha)

                ! Updating the interpolation points
                call interp_up(iinterp,rinterp,intpt,err_new,k)
              enddo ! k
              !write(*,*)'interpolation',lglel(ie),
             !& rinterp(1),rinterp(2),rinterp(3)
          ! set the new value of node ix,iy,iz
          !call set_new_pt(wrkbdry,iinterp,rinterp,intpt,ix,iy,iz,ie)
          call set_new_pt2(wrkbdry,iinterp,rinterp,intpt,
     &ix,iy,iz,sn,ie,mcx,radius)
    
          endif ! iscorner
          enddo ! ix 
          enddo ! iy 
          enddo ! iz 

          ! Redistributing the internal nodes
          call fix_gll(ie,f) 

          end if! if global node corresponds to this processor 
          enddo ! i
          ! Writing useful info to logfile 
          wrk(1) = glmax(mcx(1),1)
          wrk(2) = glmax(mcx(2),1)
          if(nid.eq.0) then
              write(*,*) 'Largest relative change:', wrk(1) 
              write(*,*) 'Largest absolute change:', wrk(2) 
            if (wrk(1).ge.1.0) then 
             write(*,*) 'WARNING: Significant change in geom \n'
            endif
          endif
c    &in geometry could lead to problems\n','
c    &try refining either surf.i or original mesh'
             
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

             error = lambda ! std L2 norm
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
      subroutine set_new_pt2(wrkbdry,iint,rint,n,ix,iy,iz,sn,ie,mcx,rad)
c     OUTPUT:
c     defining the position of the new gll-point on the surface
c     INPUT:
c     iint contains the position of the three closest points 
c     rint contains the distance to the three closest points 
c     ix,iy,iz are the indexes of the gll point in element ie
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      integer i,n,ix,iy,iz,e,f,
     &        iint(n)
      real rint(n),d,k,sn(3),tol,mcx(2),rad
      !common /surfproj/ wrkbdry(nwork,3)
      real wrkbdry(nwork,3),s,m
      integer itemp
          call nekasgn(ix,iy,iz,ie)
          !write(*,*)'interpoints', iint(1),iint(2),iint(3) 
      s = 0.15
      k = 0.0
      d = 0.0
      tol = 1e-4
      m = 2.0
      do i = 1,n
      itemp = iint(i)         
      rtemp = abs(rint(i))+tol
      d = d +                        1.0/rtemp**m 
      k = k + sn(1)*(x-wrkbdry(itemp,1))/rtemp**m
     &      + sn(2)*(y-wrkbdry(itemp,2))/rtemp**m
     &      + sn(3)*(z-wrkbdry(itemp,3))/rtemp**m
      end do

      if(k/d/(2*rad)*lx1.le.0.5) then
      xm1(ix,iy,iz,ie) = x-k*sn(1)/d*s
      ym1(ix,iy,iz,ie) = y-k*sn(2)/d*s
      zm1(ix,iy,iz,ie) = z-k*sn(3)/d*s
      mcx(1) = max(mcx(1),k/d/(2*rad)*lx1)
      mcx(2) = max(mcx(2),k/d)
      endif
      !write(*,*)'change:' , k
      return
      end
c===============================================================
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
c-----------------------------------------------------------------------
      subroutine checkcorner(is,ix,iy,iz,kz1,kz2,ky1,ky2,kx1,kx2)
          integer is,ix,iy,iz,kz1,kz2,ky1,ky2,kx1,kx2

          if(ix.eq.kx1.or.ix.eq.kx2) then
              if(iy.eq.ky1.or.iy.eq.ky2) then
                if(iz.eq.kz1.or.iz.eq.kz2) then
                    is = 1
                else 
                    is = 0
                endif
            endif
        endif

        return
        end
c-----------------------------------------------------------------------
       subroutine addcorners(wrk,ie,f)
        include 'SIZE'
        include 'GEOM'
       real wrk(4,3) 
       integer ie,ix,iy,iz,kz1,kz2,ky1,ky2,kx1,kx2,iter,f
           
       call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
       
       iter = 1
       do ix = kx1,kx2,nx1-1
       do iy = ky1,ky2,ny1-1
       do iz = kz1,kz2,nz1-1
       wrk(iter,1) = xm1(ix,iy,iz,ie)
       wrk(iter,2) = ym1(ix,iy,iz,ie)
       wrk(iter,3) = zm1(ix,iy,iz,ie)
       iter = iter + 1
       enddo
       enddo
       enddo

       return
       end






