program armc
implicit none

  integer::natom,nsteps,nprint,i,j,k,nacc=0,ilts,istep,id,istart,ncont
  real,dimension(500)::x,y,z
  real::dseed,temp,ep,sg,pe,te,xbox,ybox,zbox,mar,kb,rcut,dx,dy,dz,dd,xx
  real:: yy,zz,ep4,sg2,rcut2,tmpi,expe,kerq,scal
  real:: avpe,ape2,accr,pote,dlt,sdpe,dlte,peo,bf
  character(len=20)::gout,eout,xout,ipdb,fpdb

  !ep for epsilon & sg for sigma,pe-potential energy,te-total, av-average , ape2-average of the sqaured pe terms

  open(unit=10,file='armc.dat',status='old')
  
  read(10,*) istart ! CDT generate (0) Or from last MD run(1)
  read(10,*) ilts   ! rigid (0) OR random (1)
  read(10,*) natom  ! No.of atoms in the system
  read(10,*) nsteps ! No.of MC steps
  read(10,*) nprint ! printing interval
  read(10,*) ep     ! epsilon
  read(10,*) sg     ! sigma
  read(10,*) temp   ! intial temp
  read(10,*) expe   ! Expected total energy
  read(10,*) mar    ! Mass of Argon(g/mol)
  read(10,*) kb     ! Boltzmann constant (Kj/mol/K)
  read(10,*) rcut   ! cutoff distance
  read(10,*) dlt    ! MC step size (nm)
  read(10,*) xbox
  read(10,*) ybox   ! simulation box dimensions
  read(10,*) zbox
  read(10,'(a)') gout  ! name of the general output
  read(10,'(a)') eout  ! energy out file
  read(10,'(a)') xout  ! coordinate out file(trajectory file)
  read(10,'(a)') ipdb  ! initial CDT file in PDB format
  read(10,'(a)') fpdb  ! final CDT file in PDB format


  open(unit=20,file= gout,status='unknown')
  open(unit=21,file= eout,status='unknown')
  open(unit=22,file= xout,status='unknown')
  open(unit=23,file= ipdb,status='unknown')
  open(unit=24,file= fpdb,status='unknown')

  write(20,500)istart,ilts,natom,nsteps,nprint,ep,sg,temp,expe,mar,kb,rcut,xbox,ybox,zbox,gout,eout,xout,ipdb,fpdb
  500 format(/,2x,'start OR re-start op     =',i5,/,&
               2x,'Initial CDT option       =',i5,/,&
               2x,'number of atoms          =',i5,/,&
               2x,'No. of MC steps          =',i5,/,&
               2x,'printing interval        =',i5,/,&
               2x,'Ar epsilon               =',f8.4,/,&
               2x,'Ar sigma                 =',f8.4,/,&
               2x,'initial temperature      =',f8.4,/,&
               2x,'Expected total energy    =',f9.2,/,&
               2x,'Mass of Ar               =',f8.4,/,&
               2x,'Boltzmann Constant       =',f8.4,/,&
               2x,'cutoff distance          =',f8.4,/,&
               2x,'xbox                     =',f8.4,/,&
               2x,'ybox                     =',f8.4,/,&
               2x,'zbox                     =',f8.4,/,&
               2x,'general output file      =',a20,/,&
               2x,'energy output file       =',a20,/,&
               2x,'coordinate output file   =',a20,/,&
               2x,'initial coordinates      =',a20,/,&
               2x,'final coordinates        =',a20,/)
      
if( istart == 0 ) then   
!------------ Generating intial Coordinates

 if(ilts == 0) then
 ncont = 0
 do i =1,6
   xx = (i-1)*0.39
   do j =1,6
      yy = (j-1)*0.39 !----0.39 is the rmin value in LJ potential
      do k =1,6
         zz = (k-1)*0.39
        ncont = ncont+1
        x(ncont) = xx
        y(ncont) = yy
        z(ncont) = zz
       end do
     end do
  end do

    else
            
    ncont = 0
    99 do i = 1,5000
       xx = (rand() - 0.5)*xbox
       yy = (rand() - 0.5)*ybox
       zz = (rand() - 0.5)*zbox
       if(ncont==0) then
         ncont = ncont+1
         x(ncont) = xx
         y(ncont) = yy
         z(ncont) = zz
       else      
       do j = 1,ncont
         dx = x(j) - xx   !x(j) for previously generated coordinates   
         dy = y(j) - yy      
         dz = z(j) - zz
         dd = sqrt(dx*dx + dy*dy + dz*dz)
         if(dd < 0.39 ) then
           go to 99
        end if    
       end do ! for j loop

                
          ncont = ncont+1       
          x(ncont) = xx
          y(ncont) = yy
          z(ncont) = zz
          if (ncont == natom) go to 88
        end if
          end do !for i loop
         88 continue 
         end if

   else !istart else

     open(unit=11,file='finl.pdb',status='old')
     read(11,*) !reading the 1st line(statement)
     read(11,*) !reading 2nd line(data)
     do i= 1,natom
       read(11,'(30x,3f8.3)') xx,yy,zz
       x(i) = xx*0.10
       y(i) = yy*0.10
       z(i) = zz*0.10
     end do

    end if !for istart 
       

! ----Initial CDT file

  do i = 1,natom
    xx = x(i)*10
    yy = y(i)*10
    zz = z(i)*10
    write(23,501) i,'Ar   ','Ar   ',i,xx,yy,zz,0.0,0.0
  end do
 501 format('ATOM',2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2)  


!-------Initial potential energy and forces of the system
  ep4   = 4.0*ep 
  sg2   = sg*sg
  rcut2 = rcut*rcut
  pote  = 0.0

  do i = 1,natom
     call eng(i,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
     pote = pote + pe
  end do
  pote = 0.5*pote

  write(20,*) 'Initial potential energy = ',pote,'KJ/mol'
  write(20,*) ' ' 

  id   = 0.0
  avpe = 0.0
  ape2 = 0.0

!---------------------MC simulation

 do istep = 1,nsteps


!--------pick an atom randomly

  id = int(rand()*natom) + 1

!------------Potential Energy before the movement

 
     call eng(id,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
     peo = pe        !peo-pe old

!--------Move the atom "id" randomly


  xx = x(id)
  yy = y(id) !---storing as the new coordinates could be rejected
  zz = z(id)


  x(id) = x(id) + dlt*(rand()-0.5)
  y(id) = y(id) + dlt*(rand()-0.5)
  z(id) = z(id) + dlt*(rand()-0.5)

!--------------periodic Boundary Conditions

     x(id) = x(id) - anint(x(i)/xbox)*xbox
     y(id) = y(id) - anint(y(i)/ybox)*ybox
     z(id) = z(id) - anint(z(i)/zbox)*zbox

!---------------Call energy

  
  call eng(id,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
  dlte = pe - peo

  bf = exp(-dlte/(kb*temp))
  if (bf > rand()) then
      pote = pote + dlte
      nacc = nacc + 1
  else
      x(id) = xx
      y(id) = yy
      z(id) = zz
  end if    


   avpe = avpe + pote
   ape2 = ape2 + pote*pote
  
   if(mod(istep,nprint) == 0)then
       accr = real(nacc)/istep
      ! write(6,'(a,i5)') 'No.of steos = ',istep      
       id = id+1 
       write(21,'(i7,2f13.4)') istep,pote,accr
       write(22,'(a,i8)') 'MODEL',id
       write(22,'(a,3f8.3,3f6.1)') 'REMARKS', xbox*10,ybox*10,zbox*10,90.0,90.0,90.0
      do i = 1,natom
        xx = x(i)*10.0
        yy = y(i)*10.0
        zz = z(i)*10.0

      write(22,501) i,'Ar   ','Ar   ',i,xx,yy,zz,0.0,0.0  
      end do
      write(22,'(a)') 'ENDMOL'
      write(20,503) istep,pote,accr
      end if   
      if(mod(istep,2000)==0) write(6,'(a,i8)') 'No of steps =',istep
 end do ! ------MC run is over

 503 format(/,2x,'No. of steps     =',i7,/,&
              2x,'Potential energy =',f12.4,'KJ/mol',/,&
              2x,'Acc. Ratio       =',f12.4,'K',/)

!------Final configuration

    write(24,'(a)') 'Final configuration'
    write(24,'(a,3f8.3,3f6.1)') 'REMARKS',xbox*10,ybox*10,zbox*10,90.0,90.0,90.0
    do i = 1,natom
      xx = x(i)*10
      yy = y(i)*10
      zz = z(i)*10
      write(24,501) i,'Ar   ','Ar   ',i,xx,yy,zz,0.0,0.0
    end do
    write(24,'(a)') 'ENDMOL'
    close(unit=24)

!---------Averages and standard deviations
   sdpe = sqrt((nsteps*ape2 - avpe*avpe)/nsteps**2)
   avpe = avpe/nsteps
   write(20,505) avpe,sdpe,accr
505 format(//,2x,'Average properties--------',//,&
     2x,'Potential energy  =',f12.4,' +/- ',f8.4, 'KJ/mol',/,&
     2x,'Accr              =',f12.4,/)






  end program armc               
! ====================================================================
  subroutine eng(id,natom,mar,ep4,sg2,rcut2,x,y,z,xbox,ybox,zbox,pe)
  implicit none

  integer::i,j,k,natom,ii,id
  real,dimension(500)::x,y,z
  real::ep4,sg2,rcut2,xbox,ybox,zbox,pe,cc,c2,c6,c12,r2,xi,yi,zi,dx,dy,dz,mar 

!--------initialize--------
 
  pe = 0.0

   xi = x(id)
   yi = y(id)
   zi = z(id)
   do j = 1,natom
     if(j/=id) then
     dx = x(j) - xi  
     dy = y(j) - yi  
     dz = z(j) - zi

!--------MIC---------
    dx = dx - xbox*anint(dx/xbox)
    dy = dy - ybox*anint(dy/ybox)
    dz = dz - zbox*anint(dz/zbox)

    r2 = dx*dx + dy*dy + dz*dz
    if (r2<= rcut2) then
    c2  = sg2/r2
    c6  = c2*c2*c2
    c12 = c6*c6
    pe  = pe + ep4*(c12 - c6)
   end if
  end if 
 end do


  end subroutine
