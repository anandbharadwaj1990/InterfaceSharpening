program main

implicit none

integer :: i,j,imax, jmax
integer :: nIB, n, sub_it
real :: delx,dely, delt
real :: xc,yc,rad, tol 
real, dimension(:,:),allocatable :: xg,yg,phi,x,y, sphi, source, psi, psi_m1
real, dimension(:,:),allocatable :: iflux, jflux, divw
real, dimension(:,:,:),allocatable :: w, ncap
real :: phi_x, phi_y, eps, phi_int, div, psi_xx, psi_yy, phi_xy, kappa
real :: phix_p, phix_m, phiy_p, phiy_m, xnorm
real, dimension(2) :: w_int, n_f 
real :: coef, phi_temp, xtheta
real, dimension(:,:),allocatable :: rIB, noIB
real :: eps_i,dist, distmin, cdot, psi_f
integer :: indmin

imax = 101
jmax = 101
delx = 1.0/(1.0*(imax-1))
dely = 1.0/(1.0*(jmax-1))
delt = 1.0e-3 
eps = delx
eps_i = 8.0*delx
xtheta = 0.1

!Read in the interface and its normals 
open(10,file='LS.dat')
read(10,*)nIB
allocate(rIB(nIB,1:3))
allocate(noIB(nIB,1:3))
do n=1,nIB 
   read(10,*)rIB(n,1:3),noIB(n,1:3)
end do
close(10)

allocate(xg(imax,jmax))
allocate(yg(imax,jmax))
allocate(phi(1:imax+1,1:jmax+1))
allocate(psi(1:imax+1,1:jmax+1))
allocate(divw(1:imax+1,1:jmax+1))
allocate(psi_m1(1:imax+1,1:jmax+1))
allocate(sphi(1:imax+1,1:jmax+1))
allocate(source(1:imax+1,1:jmax+1))
allocate(x(1:imax+1,1:jmax+1))
allocate(y(1:imax+1,1:jmax+1))
allocate(w(1:imax+1,1:jmax+1,1:2))
allocate(ncap(1:imax+1,1:jmax+1,1:2))
allocate(iflux(1:imax,2:jmax))
allocate(jflux(2:imax,1:jmax))

!Generate grid
do j=1,jmax
   do i=1,imax
      xg(i,j) = (i-1)*delx
      yg(i,j) = (j-1)*dely 
   end do
end do

!Generate cell centres and get signed distance function
do j=2,jmax
   do i=2,imax
      x(i,j) = 0.25*(xg(i,j) + xg(i-1,j) + xg(i-1,j-1) + xg(i,j-1))
      y(i,j) = 0.25*(yg(i,j) + yg(i-1,j) + yg(i-1,j-1) + yg(i,j-1))
      distmin = 1.0e6
      do n=1,nIB
	 dist = sqrt((x(i,j) - rIB(n,1))**2 + (y(i,j) - rIB(n,2))**2)
	 if(dist .lt. distmin) then 
            distmin = dist 
	    indmin = n 
          end if 
      end do
      cdot = (x(i,j)-rIB(indmin,1))*noIB(indmin,1) + (y(i,j) - rIB(indmin,2))*noIB(indmin,2)
      phi(i,j) = distmin*sign(1.0,cdot) 
      if(phi(i,j) .lt. -eps_i) then 
	 psi(i,j) = 1.0 
      else if(phi(i,j) .gt. eps_i) then 
	 psi(i,j) = 0.0 
      else 
         psi(i,j) = 0.5*(tanh(-phi(i,j)/(2.0*eps_i)) + 1.0) 
      end if 
   end do
end do

!Set boundary conditions for psi 
psi(1,:) = psi(2,:) 
psi(:,1) = psi(:,2) 
psi(imax+1,:) = psi(imax,:) 
psi(:,jmax+1) = psi(:,jmax) 


!Calculate interface normals from psi
do j=2,jmax 
   do i=2,imax 
      ncap(i,j,1) = (psi(i+1,j)-psi(i-1,j))/(2.0*delx) 
      ncap(i,j,2) = (psi(i,j+1)-psi(i,j-1))/(2.0*dely) 
      if(norm2(ncap(i,j,1:2)) .gt. 1.0e-5) ncap(i,j,1:2) = ncap(i,j,1:2)/norm2(ncap(i,j,1:2))
   end do
end do

!Writing out the interface normals
open(11,file='2d_phi.dat')
do j=2,jmax 
   do i=2,imax 
      write(11,'(5e30.15)')x(i,j),y(i,j),psi(i,j),ncap(i,j,1:2)
   end do 
end do
close(11)
tol = 1.0
sub_it = 0
!Writing the initial interface
open(10,file='Initial.csv')
write(10,*)'x,y,f_lf'
do j=2,jmax
   do i=2,imax
      write(10,'(1e30.15,a,1e30.15,a,1e30.15)')x(i,j),',',y(i,j),',',psi(i,j)
   end do 
end do

close(10)

do while(tol .gt. 1.0e-8 .and. sub_it .lt. 50000)
   sub_it = sub_it + 1
   psi_m1 = psi

   !Setting zero-Neumann boundary condition
   psi(1,:) = psi(2,:) 
   psi(:,1) = psi(:,2) 
   psi(imax+1,:) = psi(imax,:) 
   psi(:,jmax+1) = psi(:,jmax) 
 
   !Calculate the fluxes for the i-face 
   do j=2,jmax 
      do i=1,imax 
         psi_f = 0.5*(psi(i,j) + psi(i+1,j))
         n_f = 0.5*(ncap(i,j,:) + ncap(i+1,j,:))        	
         if(norm2(n_f) .gt. 1.0e-5) n_f = n_f/norm2(n_f) 
         iflux(i,j) = psi_f*(1.0-psi_f)*n_f(1)*dely
      end do     
   end do 

   !Calculate the fluxes for the j-face 
   do j=1,jmax 
      do i=2,imax 
         psi_f = 0.5*(psi(i,j) + psi(i,j+1))
         n_f = 0.5*(ncap(i,j,:) + ncap(i,j+1,:))        
         if(norm2(n_f) .gt. 1.0e-5) n_f = n_f/norm2(n_f) 
         jflux(i,j) = psi_f*(1.0-psi_f)*n_f(2)*delx
      end do     
   end do

   !Calculate the source terms
   do j=2,jmax
      do i=2,imax 
         psi_xx = (psi(i+1,j) - 2.0*psi(i,j) + psi(i-1,j))/delx**2
         psi_yy = (psi(i,j+1) - 2.0*psi(i,j) + psi(i,j-1))/dely**2
         source(i,j) = eps*(psi_xx + psi_yy)
      end do 
   end do 
 
   !March the equation using pseudo time stepping 
   do j=2,jmax 
      do i=2,imax 
   	psi(i,j) = psi(i,j) - delt/(delx*dely)*(iflux(i,j) - iflux(i-1,j) + jflux(i,j) - jflux(i,j-1)) + & 
			      delt*source(i,j) 
      end do 
   end do 

   tol = maxval(abs(psi-psi_m1))
   print*,sub_it,tol
end do

!Writing out the final interface after sharpening
open(10,file='Final.csv')
write(10,*)'x,y,f_lf'
do j=2,jmax
   do i=2,imax
      write(10,'(1e30.15,a,1e30.15,a,1e30.15)')x(i,j),',',y(i,j),',',psi(i,j)
   end do 
end do

close(10)
end program main
