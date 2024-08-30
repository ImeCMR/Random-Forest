 program radialdis
 implicit none
!
! Defining variables
!
 integer:: i,j,ii,ibin,nbin,nconfig,iconfig,natom
 real,dimension(500)::x,y,z,bin
 real::distance,rr,xbox,ybox,zbox,xi,yi,zi,dx,dy,dz,binSize,volume,rho,pi,dv,gr
 character(len=20)::icdt,ordf
!
! Open input data file
!
 open(unit=10,file='rdf.dat',status='old')
!
! Reading the file
!
 read(10,*) natom
 read(10,*) nbin
 read(10,*) nconfig
 read(10,*) distance
 read(10,'(a20)') icdt          ! Input coordinate file (armd.cdt)
 read(10,'(a20)') ordf          ! Output file name
!
! Calculating the bin size
!
 binSize = distance/nbin
!
! Intialise all the bins
!
 do i = 1,nbin+5
   bin(i) = 0.0
 end do
!
! Read the armd.cdt file
!
 open(unit=11,file='armd.cdt',status='old')
 do iconfig = 1,nconfig             !  going through all the configurations
   read(11,*)
   read(11,'(7x,3f8.3)') xbox,ybox,zbox
   do i = 1,natom
     read(11,'(30x,3f8.3)') x(i),y(i),z(i)
   end do
   read(11,*)
!------------------------- iconfig_th configuration is been read
!
!------------RDF calculation starts here
!
   do i = 1,natom-1
     xi = x(i)
     yi = y(i)
     zi = z(i)
     ii = i+1
     do j = ii,natom
       dx = x(j) - xi
       dy = y(j) - yi
       dz = z(j) - zi
!
!----------Minimum Image Convension
!
       dx = dx - anint(dx/xbox)*xbox
       dy = dy - anint(dy/ybox)*ybox
       dz = dz - anint(dz/zbox)*zbox
!
!----------Calculating the distances (rr)
!
       rr = sqrt(dx*dx + dy*dy + dz*dz)
       if (rr <= distance) then
         ibin = int(rr/binSize) + 1
         bin(ibin) = bin(ibin) + 2             ! This works if all the atoms are identical
       end if  
     end do
   end do
 end do
!
!
 volume = xbox*ybox*zbox
 rho = natom/volume
 pi = 22.0/7.0
!
! Average number of Ar atoms in the bin
!
 do i = 1,nbin
   bin(i) = bin(i)/natom/nconfig
   rr = (i-1)*binSize + 0.5*binSize
   dv = 4.0*pi*rr*rr*binSize
   gr = bin(i)/dv/rho
   write(20,'(f8.4,f12.6)') rr,gr
 end do
!
!
 end program radialdis
