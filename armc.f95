  program armc
  implicit none
!
  integer::natom,i,j,k,nprint,nsteps,ncont,ilts,isteps,istart,id,nacc=0
  real,dimension(500)::x,y,z
  real::dseed,temp,ep,sg,pe,ke,te,xbox,ybox,zbox,xx,yy,zz,dx,dy,dz,dd,mar
  real:: kb,rcut,ep4,sg2,rcut2,tmpi,expe,kerq,cal,scal,avpe,avke,ape2,sdpe,sdke,sdtmp
  real:: sdte,accr,pote,dlt,dlte,peo,bf
  character(len=20)::gout,eout,xout,ipdb,fpdb

! averge => avpe,avte so on
! sd     => standard dev

!
  open(unit=10,file='armc.dat',status='old')
!
  read(10,*) istart               ! CDT generate (0) OR from last MD run (1)
  read(10,*) ilts                 ! rigid (0) OR random (1)
  read(10,*) natom                ! No of atoms in the system
  read(10,*) nsteps                ! No of MD steps
  read(10,*) nprint               ! Printing intervels
  read(10,*) ep                   ! epsilon
  read(10,*) sg                   ! sigma
  read(10,*) temp                 ! initial temperature
  read(10,*) expe                 ! Expected total energy (kJ/mol)
  read(10,*) mar                  ! mass of the Ar (g/mol)
  read(10,*) kb                   ! Boltzman constant (kJ/mol/K)
  read(10,*) rcut                 ! Cutoff distance
  read(10,*) dlt                   ! MC  step size (nm)
  read(10,*) xbox                 !
  read(10,*) ybox                 ! Simulation box dimension
  read(10,*) zbox                 !
  read(10,'(a)') gout             ! general out file 
  read(10,'(a)') eout             ! energy out file
  read(10,'(a)') xout             ! coordinate out file 
  read(10,'(a)') ipdb             ! Initial Coordinate out file in PDB format
  read(10,'(a)') fpdb             ! final CDT file in PDB format
!
  open(unit=20,file=gout,status='unknown')
  open(unit=21,file=eout,status='unknown')
  open(unit=22,file=xout,status='unknown')
  open(unit=23,file=ipdb,status='unknown')
  open(unit=24,file=fpdb,status='unknown')
!
  write(20,500) istart,ilts,natom,nsteps,nprint,ep,sg,temp,expe,mar,kb,rcut,&
                xbox,ybox,zbox,gout,eout,xout,ipdb,fpdb
 500 format(/,2x,' Start OR re-start op    = ',i5,/,&
              2x,' Rigid or solid          = ',i5,/,&
              2x,' Number of atoms         = ',i5,/,&
              2x,' No of MC steps          = ',i7,/,&
              2x,' Printing intervels      = ',i5,/,&
              2x,' Ar epsilon              = ',f8.4,/,&
              2x,' Ar sigma                = ',f8.4,/,&
              2x,' Initial tempereture     = ',f8.4,/,&
              2x,' Expected Tot Energy     = ',f9.2,/,&
              2x,' Mass of Ar              = ',f8.4,/,&
              2x,' Boltzman constant       = ',f8.6,/,&
              2x,' Cutoff distance         = ',f8.4,/,&
              2x,' xbox                    = ',f8.4,/,& 
              2x,' ybox                    = ',f8.4,/,& 
              2x,' zbox                    = ',f8.4,/,&
              2x,' Genaral out file        = ',a20,/,&
              2x,' Energy out file         = ',a20,/,&
              2x,' Coordinate out file     = ',a20,/,&
              2x,' Coordinate out file     = ',a20,/,&
              2x,' Final coordinate        = ',a20,/)

 if ( istart ==0) then
!
!---------Genarating Initial Coordinates
 
  if( ilts == 0) then
  ncont = 0
  do i = 1,6
    xx = (i-1)*0.39
    do j = 1,6
     yy = (j-1)*0.39
      do k = 1,6
        zz = (k-1)*0.39
        ncont =ncont + 1
        x(ncont) = xx
        y(ncont) = yy
        z(ncont) = zz
       end do  
    end do
  end do
!
  else
!
   ncont = 0
   99 do i = 1,5000
        xx = (rand() - 0.5)*xbox     ! pick a point randomly
        yy = (rand() - 0.5)*ybox
        zz = (rand() - 0.5)*zbox
        if( ncont == 0) then
          ncont = ncont + 1
          x(ncont) = xx
          y(ncont) = yy
          z(ncont) = zz
        else
        do j = 1,ncont
         dx = x(j) - xx     
         dy = y(j) - yy     
         dz = z(j) - zz
         dd = sqrt(dx*dx + dy*dy + dz*dz )     
         if( dd < 0.39 ) go to 99
        end do
         ncont = ncont + 1
         x(ncont) = xx
         y(ncont) = yy
         z(ncont) = zz
!      if(mod(ncont,25) == 0) write(6,'(i5)') ncont
         if( ncont ==natom) go to 88
     end if
   end do
   88 continue
  end if
!
  else          ! for istart
          
   open(unit=11,file='finl.pdb',status='old')
   read(11,*)
   read(11,*)
   do i = 1,natom
      read(11,'(30x,3f8.3)') xx,yy,zz
      x(i) = xx*0.10    
      y(i) = yy*0.10    
      z(i) = zz*0.10
   end do
  end if
    
!          
!-Initial Coordinate file------
!
  do i = 1,natom
    xx = x(i)*10
    yy = y(i)*10
    zz = z(i)*10
    write(23,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0
  end do  
 501 format('ATOM',2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2)

! 
!
! ----- Initial potential energy of the system
!
  ep4   = 4.0*ep
  sg2   = sg*sg
  rcut2 = rcut*rcut 
  pote  = 0.0
!
  do i = 1,natom
     call eng(i,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
     pote = pote + pe
  end do 

  pote = 0.5*pote             ! because of do loop

!  write(20,'(/,a,f12.5,a)') ' Initial potential energy = ',pe,' kJ/mol'
  write(20,*) ' initial potential energy = ',pote,' kJ/mol'
  write(20,*) ' '

!------------------------------------------------------------------------------
!  
!

  id   = 0
  avpe  = 0.0
  ape2  = 0.0
!  
!------Monte carlo run starts here
!
  do isteps = 1,nsteps
!
!---------- Pick an atom randomly


    id = int(rand()*natom) + 1

!-------------Potential energy before the mom=vement
    
      call eng(id,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
 
      peo = pe 

!------Move the atom "ID' randomly

    xx = x(id)  ! if movement rejected store the coordinates
    yy = z(id) 
    zz = z(id)

!

    x(id) = x(id) + dlt*(rand()-0.5)      ! new coordinate
    y(id) = y(id) + dlt*(rand()-0.5)
    z(id) = z(id) + dlt*(rand()-0.5)


!----periodic boundry condition
!
      x(id) = x(id) - anint(x(i)/xbox)*xbox
      y(id) = y(id) - anint(y(i)/ybox)*ybox
      z(id) = z(id) - anint(z(i)/zbox)*zbox
!
!---call eng
!
      call eng(id,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)  
      dlte = pe - peo
      
!
      bf   = exp(-dlte/(kb*temp))
      if ( bf > rand() ) then
        pote = pote + dlte
        nacc = nacc + 1
      else
        x(id) = xx
        y(id) = yy
        z(id) = zz
      end if  


!
!
      avpe  = avpe + pote
      ape2  = ape2 + pote*pote
!-------------------------------------------------------------------------------
    
      if(mod(isteps,nprint) == 0) then
        accr  = real(nacc)/isteps
!        write(6,'(a,i5)') ' No of steps = ',isteps     
        id = id + 1     
        write(21,'(i7,2f13.4)') isteps,pote,accr
        write(22,'(a,i8)') 'MODEL',id
        ! box length and angles
        write(22,'(a,3f8.3,3f6.1)') 'REMARKS',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0 
        do i = 1,natom                 ! to write coordinates
          xx = x(i)*10.0             ! convert in to angstrom (nm-->angstrom)
          yy = y(i)*10.0             ! convert in to angstrom (nm-->angstrom)
          zz = z(i)*10.0             ! convert in to angstrom (nm-->angstrom)
          write(22,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0
        end do
        write(22,'(a)') 'ENDMDL'
        write(20,503) isteps,pote,accr
     end if        
     if(mod(isteps,2000)==0) write(6,'(a,i8)') ' No of steps = ',isteps
  end do  
!
! ------------------------MC run is over here
!
  503 format(/,2x,' No of steps              = ',i7,/,&
               2x,' Potential Energy         = ',f12.4,' kJ/mol',/,&
               2x,' Acceptance Ratio         = ',f12.4,' K',/)
!
!--------final configuration
!
   write(24,'(a)') ' Final Configuration '
   write(24,'(a,3f8.3,3f6.1)') 'REMARKS ',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0
   do i = 1,natom
     xx = x(i)*10 
     yy = y(i)*10 
     zz = z(i)*10 
     write(24,501) i,'Ar  ','Ar  ',i,xx,yy,zz,0.0,0.0
   end do
   write(24,'(a)') 'ENDMOL'
   close(unit=24)  
!
!-----------Averages and stander deviations
!
    sdpe   = sqrt((nsteps*ape2 - avpe*avpe)/(nsteps*nsteps))
    avpe  = avpe/nsteps
    write(20,505) avpe,sdpe,accr

 505 format(//,2x,' Average properties-------------------- ',//,&
               2x,' Potential energy    = ',f12.4,' +/-',f8.4,' kJ/mol',/,& 
               2x,' Acceptence ration   = ',f12.4,/) 
  end program
!***************************************************************************


  subroutine eng(id,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
  implicit none
!

  integer::i,j,k,natom,ii,id
  real,dimension(500)::x,y,z
  real::ep4,sg2,rcut2,xbox,ybox,zbox,pe,cc,c2,c6,c12,r2,xi,yi,zi,dx,dy,dz,mar
!
!------Initialize
    pe = 0.0
    xi   = x(id)
    yi   = y(id)
    zi   = z(id)
    do j = 1,natom
      if( j /= id ) then
       dx = x(j) - xi
       dy = y(j) - yi
       dz = z(j) - zi
!---Minimum image convension
       dx = dx - xbox*anint(dx/xbox)
       dy = dy - ybox*anint(dy/ybox)
       dz = dz - zbox*anint(dz/zbox)
!
       r2 = dx*dx + dy*dy + dz*dz
       if( r2 <= rcut2 ) then
         c2  = sg2/r2
         c6  = c2*c2*c2
         c12 = c6*c6
         pe  = pe + ep4*(c12-c6)
       end if 
      end if 
    end do

!
  end subroutine
