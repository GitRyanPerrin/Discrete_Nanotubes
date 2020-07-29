!==================================================================================
! CiDi = Circular Geometry with Discrete Model
! v03 : includes magnetic field
!==================================================================================
  MODULE Mod_Precision
    implicit none
    integer, parameter :: dp=8   ! Double precission parameter
  END MODULE Mod_Precision
!==================================================================================
  MODULE Mod_Data      
    use Mod_Precision
    implicit none
 ! general data
    real(kind=dp),    parameter :: pi  = 3.14159265358979324_dp
    complex(kind=dp), parameter :: ci=(0.0_dp,1.0_dp)   ! complex i
 ! input sample data
!   real(kind=dp), parameter :: Rext=1.0_dp, Rint=0.6_dp ! exterior and interior radii
    real(kind=dp), parameter :: Rext=1.0_dp, Rint=0.8_dp ! exterior and interior radii
!   real(kind=dp), parameter :: Rext=1.0_dp, Rint=0.95_dp ! exterior and interior radii
    integer,       parameter :: Nr=10, Nphi=50   ! number of rings and angles
!   integer,       parameter :: Nr=1, Nphi=500 ! number of rings and angles
    integer,       parameter :: Nsites=Nphi*Nr    ! number of sites (or states)
!   integer,       parameter :: Nsites=Nphi*(Nr-1)+1+(Nphi-1)*Nint(Rint/(Rint+0.0001))  ! includes Rint=0
    integer, parameter :: Nstat=Nsites*2      ! number of single particle states with spin
    real(kind=dp), parameter :: ts=1.0_dp            ! energy unit = hbar**2/(2*m*Rext**2)
!   real(kind=dp), parameter :: gamma=-0.01474_dp  ! GaAs geff=-0.44, meff=0.067, gamma=geff*meff/2
!   real(kind=dp), parameter :: gamma=-0.17_dp  ! InAs geff=-14.9, meff=0.023, gamma=geff*meff/2
    real(kind=dp), parameter :: gamma=-0.2_dp  ! InAs geff=-14.9, meff=0.023, gamma=geff*meff/2
!   real(kind=dp), parameter :: gamma=0.0_dp  ! InAs geff=-14.9, meff=0.023, gamma=geff*meff/2
!   real(kind=dp), parameter :: emag=2.50_dp       ! Cyclotron energy in units ts
!   real(kind=dp), parameter :: alphaR=0.0_dp, betaD=0.0_dp ! (SOI param)/Rext in units ts
!   real(kind=dp), parameter :: alphaR=0.0_dp, betaD=0.0_dp ! (SOI param)/Rext in units ts
!   real(kind=dp), parameter :: alphaR=0.2_dp, betaD=0.1_dp ! (SOI param)/Rext in units ts
!   real(kind=dp), parameter :: alphaR=0.70_dp, betaD=0.70_dp ! (SOI param)/Rext in units ts
    real(kind=dp), parameter :: alphaR=0.7_dp, betaD=0.3_dp ! (SOI param)/Rext in units ts
!   real(kind=dp), parameter :: uc0=2.2000_dp     ! Coulomb interaction strength
    real(kind=dp), parameter :: uc0=0.0000_dp     ! Coulomb interaction strength
    real(kind=dp), parameter :: eta=0.001_dp      ! "screening" parameter
    real(kind=dp) emag       ! Cyclotron energy in units ts
  END MODULE Mod_Data      
!==================================================================================
  MODULE Mod_Single_Particle_States
    use Mod_Data
    implicit none
    real(kind=dp)  d_r, d_phi                 ! radial and angular interval
    real(kind=dp), dimension(Nsites,4) :: pos ! position of sites in polar and cartezian coord.
! pos(k,1)=r, pos(k,2)=phi, pos(k,3)=x, pos(k,4)=y
    integer, dimension(Nstat,2) :: qn         ! quantum numbers
! qn(a,1)=k=site index,   qn(a,2)=spin
    real(kind=dp), dimension(Nstat)          :: Ea        ! single particle energies
    complex(kind=dp), dimension(Nstat,Nstat) :: Hmat, Psi ! Hamiltonian single and eigenvector matrix
    integer, dimension(Nstat)                :: Pari      ! Parity of single particle states
  END MODULE Mod_Single_Particle_States
!==================================================================================
  MODULE Mod_Many_Particle_States
    use Mod_Single_Particle_States
    integer, parameter :: Nel=1                  ! number of electrons
!   integer, parameter :: Nel=5                  ! number of electrons
    integer, parameter :: Nsing=12               ! number of single-particle states used
    integer  Nmany, Nmanyg                 ! number of many-body states and degenerated groups
    integer, allocatable, dimension(:,:) :: ioc  ! occupation numbers of single particle states 
    integer, allocatable, dimension(:,:) :: cmat ! matrix c_a |nu>
    integer, allocatable, dimension(:,:,:) :: ccmat ! matrix c_a*c_b |nu>
    complex(kind=dp), dimension(Nsing,Nsing,Nsing,Nsing) :: VCmat ! Coulomb integrals
    real(kind=dp), allocatable, dimension(:)      :: EEa         ! many particle energies
    complex(kind=dp), allocatable, dimension(:,:) :: HHmat, PPsi ! Hamiltonian many and eigenvector matrix
    integer, allocatable, dimension(:)            :: PPari      ! Parity of many particle states
    integer, allocatable, dimension(:) :: Deg_group,Deg_order ! Degeneracy info 
    integer, parameter :: Ncsd=6  ! number of many-body states to calc. charge and spin density
    real(kind=dp), dimension(Ncsd,Nsites,0:3)     :: CSd,CSdg ! charge-spin density, inlcuding degeneracy
    real(kind=dp), dimension(Ncsd,3)     :: CDWamp ! amplitude of the CDW
  END MODULE Mod_Many_Particle_States
!==================================================================================
! end modules
!==================================================================================
  PROGRAM main
    use Mod_Single_Particle_States
    use Mod_Many_Particle_States
    implicit none
    integer k,mu,nplus,nminus
    character*2 filenumber
    real(kind=dp) Epar(2*Ncsd,2)   ! energies separated by parities

    complex(kind=dp) Sigmat

    call Open_Output_Files

    write(16,*)'# emag, Sz' 
    write(20,*)'# emag, EEa(mu),mu=1,2*Ncsd '

    if(Nel.eq.1.and.2*Ncsd.gt.Nsing) then 
       write(*,*) 'Nel=1 and Nsing < 2*Nscd.  I need > or ='
       stop
    endif 

  do emag=0,10,0.025_dp
! do emag=0.2,0.2,0.01_dp

    call Lattice_Sites
    do k=1,Nsites
!     write(*,*) k, pos(k,3),pos(k,4)
      write(10,*) pos(k,3),pos(k,4)
    enddo
 
    call Quantum_Numbers

    call Hamiltonian_Single
 
    call Single_Particle_Eigenstates

    call Occupation_Numbers

    call CC_Matrix

    call Coulomb_integrals

    call Hamiltonian_Many

    call Many_Particle_Eigenstates

    call Info_Many_Particle_Eigenstates

    call Charge_Spin_Density 

    write(20,50) emag,(EEa(mu),mu=1,2*Ncsd)

    nplus=0; nminus=0
    if(emag.gt.0.) then
      do mu=1,2*Ncsd
!       if(Pari(mu).gt.0) then 
        if(PPari(mu).gt.0) then 
           nplus=nplus+1; Epar(1,nplus)=EEa(mu)
        else
           nminus=nminus+1; Epar(2,nminus)=EEa(mu)
        endif   
      enddo
      write(71,50) emag,(Epar(1,mu),mu=1,nplus)
      write(72,50) emag,(Epar(2,mu),mu=1,nminus)
    endif

    write(16,50) emag,(sum(CSd(mu,1:Nsites,3)), mu=1,Ncsd) ! Sz

    do mu=1,Ncsd
      write(300+mu,100) emag, CDWamp(mu,1),CDWamp(mu,2),CDWamp(mu,3)
    enddo

    call Deallocate_Open_Matrices

    call Close_Output_Files
     
 enddo

   stop 

 50 format(20(e14.7,2x))
100 format(4(e14.7,2x))      

  END PROGRAM main
!==================================================================================
  SUBROUTINE Lattice_Sites     ! coordinates of the lattice sites
  use Mod_Single_Particle_States
  implicit none
  integer k,iphi,ir
  real(kind=dp) r,phi,x,y
  k=0    ! site index
! pos(k,1)=r, pos(k,2)=phi, pos(k,3)=x, pos(k,4)=y
  r=Rext
  d_r=0
  if(Nr.gt.1) d_r=(Rext-Rint)/(Nr-1)   ! radial interval
  d_phi=2*pi/Nphi                      ! angular interval
  do ir=1,Nr
!   if(Nr.gt.1) r=Rext-(Rext-Rint)*(ir-1)/(Nr-1)    ! radial pos.
    if(Nr.gt.1) r=Rext-d_r*(ir-1) ! radial pos.
    do iphi=1,Nphi
       phi=d_phi*(iphi-1)    ! angular pos.
       k=k+1
       pos(k,1)=r ;  pos(k,2)=phi
       x=r*cos(phi); y=r*sin(phi)
       pos(k,3)=x ;  pos(k,4)=y
       if(r.eq.0) then 
         return
       endif  
    enddo
  enddo
  END SUBROUTINE Lattice_Sites
!==================================================================================
  SUBROUTINE Quantum_Numbers   ! define quantum numbers qn(a,1)=k  qn(a,2)=spin
  use Mod_Single_Particle_States
  implicit none
  integer a,k,spin
  a=0
! do spin=1,-1,-2 ! separate spin blocks, more convenient to see opposite spin states
                  ! but possibly the diagonalization subroutine may not find degeneracies
    do k=1,Nsites
    do spin=1,-1,-2 ! opposite spin states can be seen only if spin degeneracy is lifted
      a=a+1
      qn(a,1)=k
      qn(a,2)=spin
    enddo
  enddo
  END SUBROUTINE Quantum_Numbers   
!==================================================================================
  SUBROUTINE Hamiltonian_Single
  use Mod_Single_Particle_States
  implicit none
  integer a1,a2,k1,k2,spin1,spin2,kont
  real (kind=dp) r1,phi1,x1,y1
  real (kind=dp) r2,phi2,x2,y2
  real (kind=dp) rr,pp,tphi,tr,phimed
  complex(kind=dp) SigMat,term
  integer j1,j2,jj,semn
  integer n1,n2,nn

  Hmat=(0.0_dp,0.0_dp)
! Hmat(1,1)=1.0_dp   ! impurity
! Hmat(2,2)=1.0_dp   ! impurity


  do a1=1,Nstat
    k1=qn(a1,1); spin1=qn(a1,2)                      ! site, spin
    r1=pos(k1,1); phi1=pos(k1,2)                     ! polar coord
    n1=Nint((Rext-r1)/d_r+1); j1=Nint(phi1/d_phi)+1  ! polar coord discretized
    if(d_r.eq.0) n1=1
    tphi=ts/d_phi**2*(Rext/r1)**2                    ! the angular hopping energy 

! basic Hamiltonian, diagonal terms: orbital magnetic r**2 and Zeeman
    Hmat(a1,a1)=Hmat(a1,a1)+2*tphi+(emag*r1/4)**2+gamma*emag*spin1/2  

! periodic angular potential
!  Hmat(a1,a1)=Hmat(a1,a1) - 0.1*sin(2*phi1)    ! V_per=V*sin(n*phi)
!  if(spin1.eq.-1) Hmat(a1,a1)=Hmat(a1,a1) + 0.1*sin(1*phi1)    ! V_per=V*sin(n*phi)

! linear potential
!  Hmat(a1,a1)=Hmat(a1,a1) + 0.5*pos(k1,3)    ! V_lin=V*x1  x1=pos(k1,3)
!  Hmat(a1,a1)=Hmat(a1,a1) + 0.1*pos(k1,3)    ! V_lin=V*x1  x1=pos(k1,3)

    do a2=a1,Nstat        ! use Hermiticity
!   do a2=1,Nstat         ! do not use Hermiticity
      k2=qn(a2,1); spin2=qn(a2,2)                      ! site, spin
      r2=pos(k2,1); phi2=pos(k2,2);                    ! polar coord
      n2=Nint((Rext-r2)/d_r+1); j2=Nint(phi2/d_phi)+1  ! polar coord discretized
      if(d_r.eq.0) n2=1
      rr=abs(r1-r2); pp=abs(phi1-phi2);                ! relative distances
      nn=abs(n1-n2); jj=abs(j1-j2)                     ! relative distances discretized

! basic Hamiltonian, kinetic hopping terms / begin
      if(spin1.eq.spin2) then  
! angular hopping
        if(n1.eq.n2) then 
          if(jj.eq.1.or.jj.eq.(Nphi-1)) Hmat(a1,a2)=Hmat(a1,a2)-tphi 
        endif
! radial hopping
        if(Nr.gt.1) then 
          tr=ts*(Rext/d_r)**2 
          if(jj.eq.0.and.nn.eq.0) Hmat(a1,a2)=Hmat(a1,a2)+2*tr 
          if(jj.eq.0.and.nn.eq.1) then 
            Hmat(a1,a2)=Hmat(a1,a2)-tr 
          endif
        endif
      endif
! basic Hamiltonian, kinetic hopping terms / end

! Rashba SOI / begin
! angular R SOI
      if( n1.eq.n2.and.( jj.eq.1.or.jj.eq.(Nphi-1) ) ) then 
        semn=1
        if((jj.eq.1.and.j1.gt.j2).or.(j1.eq.1.and.j2.eq.Nphi)) semn=-1
!       if(j1.gt.j2.or.(j1.eq.1.and.j2.eq.Nphi)) semn=-1

!       write(*,*)'a1,a2',a1,a2
!       write(*,*)'semn', j1,j2,semn,jj
!       write(*,*)

        term=(SigMat('r',spin1,spin2,phi1)+SigMat('r',spin1,spin2,phi2))/2
!       phimed=(phi1+phi2)/2                         ! version with average phi
!       if (j1.eq.1.and.j2.eq.Nphi) phimed=-d_phi/2  ! version with average phi
!       term=SigMat('r',spin1,spin2,phimed)          ! version with average phi
        term=ci*alphaR/(2*r1*d_phi)*term
        Hmat(a1,a2)=Hmat(a1,a2)-term*semn   
        
!       if(spin1.eq.spin2) Hmat(a1,a2)=Hmat(a1,a2)+0.25*ci*emag/d_phi*semn ! orbital magnetic d/(d phi) semn gresit
        if(spin1.eq.spin2) Hmat(a1,a2)=Hmat(a1,a2)-0.25*ci*emag/d_phi*semn ! orbital magnetic d/(d phi)

      endif
! radial R SOI
      if(j1.eq.j2.and.Nr.gt.1.and.nn.eq.1) then 
        semn=1
        if(n1.gt.n2) semn=-1
        term=SigMat('p',spin1,spin2,phi1)
        term=ci*alphaR/(2*d_r)*term
        Hmat(a1,a2)=Hmat(a1,a2)+term*semn     
      endif
! magnetic R SOI
       if(n1.eq.n2.and.j1.eq.j2) then 
         term=0.25*alphaR*emag*r1*SigMat('r',spin1,spin2,phi1)
         Hmat(a1,a2)=Hmat(a1,a2)+term
       endif

! Rashba SOI / end

! Dresselhaus SOI / begin
! angular D SOI
      if( n1.eq.n2.and.( jj.eq.1.or.jj.eq.(Nphi-1) ) ) then 
        semn=1
        if((jj.eq.1.and.j1.gt.j2).or.(j1.eq.1.and.j2.eq.Nphi)) semn=-1
!       if(j1.gt.j2.or.(j1.eq.1.and.j2.eq.Nphi)) semn=-1
        term=(SigMat('p',spin1,spin2,phi1)+SigMat('p',spin1,spin2,phi2))/2
        term=ci*betaD/(2*r1*d_phi)*Conjg(term)
        Hmat(a1,a2)=Hmat(a1,a2)-term*semn   
      endif
! radial D SOI
      if(j1.eq.j2.and.Nr.gt.1.and.nn.eq.1) then 
        semn=1
        if(n1.gt.n2) semn=-1
        term=-SigMat('r',spin1,spin2,phi1)
        term=ci*betaD/(2*d_r)*Conjg(term)
        Hmat(a1,a2)=Hmat(a1,a2)+term*semn     
      endif
! magnetic D SOI
       if(n1.eq.n2.and.j1.eq.j2) then 
         term=0.25*betaD*emag*r1*Conjg(SigMat('p',spin1,spin2,phi1))
         Hmat(a1,a2)=Hmat(a1,a2)+term
       endif

! Dresselhaus SOI / end


     Hmat(a2,a1)=Conjg(Hmat(a1,a2)) ! Hermitian conjugation

    enddo
  enddo

! Test Hermiticity
! write(*,*)'Test Herimiticity H single'
!  do a1=1,Nstat; do a2=a1,Nstat
!  if(Hmat(a1,a2).ne.Conjg(Hmat(a2,a1)))  then 
!      write(*,*)'states', a1,a2
!      write(*,*)'k, s, ',qn(a1,1), qn(a1,2), qn(a2,1),qn(a2,2)
!      write(*,*) Hmat(a1,a2)
!      write(*,*) Hmat(a2,a1)
!  endif
!  enddo;enddo
!  stop

  END SUBROUTINE Hamiltonian_Single
!==================================================================================
  SUBROUTINE Single_Particle_Eigenstates
! eigenstates of the single particle Hamiltonian
  use Mod_Many_Particle_States
!  use mkl95_lapack      
  use lapack95      
!  use libmkl_lapack95_lp64
  implicit none
  real(kind=dp) norm,sz
  integer a1,a2,ierr,spin,k
  real(kind=dp), dimension (Nsing,Nsites) :: SPdens  ! single particle density

  Psi=Hmat
! call heev(a=Psi, w=Ea, jobz='V', uplo='U', info=ierr)   ! wrong norm of Psi in degenerate cases
  call heevd(a=Psi, w=Ea, jobz='V', uplo='U', info=ierr)  ! divide-and-conquer algorithm
 
! write(*,*) 'ierr=',ierr

! test eigenstates
  SPdens=0
  do a1=1,Nstat
    write(*,*) a1,Ea(a1)
    norm=sqrt(sum(abs(Psi(1:Nstat,a1))**2))
    if(abs(norm-1.0_dp).gt.0.01) then 
       write(*,*)'WARNING norm of single-particle states not 1 for state', a1,norm
       stop
    endif
    if(a1.le.Nsing) then 
       do a2=1,Nstat
         k=qn(a2,1); spin=qn(a2,2)
         SPdens(a1,k)=SPdens(a1,k)+abs(Psi(a2,a1))**2
       enddo
    endif
!   do a2=1,Nstat
!   Pari(a1)=0
!   sz=Real(Psi(a2,a1))*Real(Psi(a2+Nphi,a1))
!   if(a2.lt.Nphi) Pari(a1)=Nint(sign(1.,sz))
!   write(12,60) a1,a2,Psi(a2,a1),pos(qn(a2,1),1),pos(qn(a2,1),2)/pi,qn(a2,2),Pari(a1)*qn(a2,2)
!   enddo
!   write(12,*)
! Parity
   sz=Real(Psi(1,a1))*Real(Psi(1+Nphi,a1)); Pari(a1)=Nint(sign(1.,sz))
! Spin z
    sz=0
    do a2=1,Nstat
    sz=sz+abs(Psi(a2,a1))**2*qn(a2,2)
    enddo
    write(11,*) a1,Ea(a1),sz,Pari(a1)

  enddo

! do k=1,Nsites
!   write(12,50) k,(SPdens(a1,k),a1=1,Nsing) 
! enddo
  50 format(i4,50(2x,e14.7))
  60 format(i5,2x,i5,4x,e14.7,1x,e14.7,4x,f8.4,2x,f8.4,2x,i2,2x,i2)

  END SUBROUTINE Single_Particle_Eigenstates
!==================================================================================
  SUBROUTINE Parity_Single_Particle_Eigenstates
  use Mod_Many_Particle_States
  END SUBROUTINE Parity_Single_Particle_Eigenstates
!==================================================================================
  SUBROUTINE Occupation_Numbers
! calculate occupation numbers for each many particle state
  use Mod_Many_Particle_States
  implicit none
  integer m,ia(Nsing),npart,a,alpha,ierr
  real(kind=dp) x
  write(*,*)
  write(*,*) 'Nel, Nsing ',Nel, Nsing

! calculate Nmany
  x=1.0_dp
  do a=0,Nel-1
    x=x*(Nsing-a)/(Nel-a)
  enddo
  Nmany=Nint(x)
  write(*,*) 'Nmany ',Nmany

! calculate occupation numbers for each many particle state
  allocate(ioc(Nsing,Nmany), stat=ierr)
  alpha=0
  do m=0,2**Nsing    ! index for many-body states with unspecified number of particles
    call CONV10(1,m,2,Nsing,ia)    ! convert m from decimal to binary
    npart=sum(ia(1:Nsing))
    if(npart.eq.Nel) then 
      alpha=alpha+1  ! index for basis many-body states with Nel particles
      do a=1,Nsing
        ioc(a,alpha)=ia(a)   ! occupation numbers for many-body state alpha
      enddo
    endif   
  enddo
  if(alpha.ne.Nmany) then 
    write(*,*) 'ERROR wrong counting of many-body states'
    stop
  endif
  write(*,*) 
!  do alpha=1,Nmany
!    write(*,10) (ioc(a,alpha),a=1,Nsing)
!10 format(10(i1))
!   enddo
  END SUBROUTINE Occupation_Numbers
!==================================================================================
  SUBROUTINE CC_Matrix
! calculate matrix ccmat(a,b,alpha)=c_a*c_b |alpha>
! calculate matrix cmat(a,alpha)=c_a |alpha>  used in Charge_Spin_Density 
  use Mod_Many_Particle_States
    implicit none
    integer,dimension(Nsing) :: jcc
    integer alpha,s,a,b,ksgn,n,ierr

!  matrix ccmat / begin

    allocate(ccmat(Nsing,Nsing,Nmany), stat=ierr)
    do 10 alpha=1,Nmany
      do a=1,Nsing
        jcc(a)=ioc(a,alpha)    ! store occupation numbers
      enddo
      do 20 a=1,Nsing
      do 20 b=1,a
        ccmat(a,b,alpha)=0
        if(a.eq.b) go to 20
        if(jcc(a).eq.0.or.jcc(b).eq.0) go to 20
        jcc(a)=0 ; jcc(b)=0    ! annnihilation on both states
        call CONV10(-1,n,2,Nsing,jcc)
        ccmat(a,b,alpha)=n+1  ! converts to decimal the new bit string
                              ! n+1 to avoid n=0
        jcc(a)=1 ; jcc(b)=1   ! restore initial value
        ! calc sign prefactor / begin
         ksgn=0
           if(a.gt.1) then
              do s=1,a-1
                ksgn=ksgn+jcc(s) 
              enddo
           endif
           if(b.gt.1) then
              do s=1,b-1
                ksgn=ksgn+jcc(s) 
              enddo
           endif
        ! calc sign prefactor / end
      ccmat(a,b,alpha)=(-1)**(ksgn-1)*ccmat(a,b,alpha)
   20 ccmat(b,a,alpha)=-ccmat(a,b,alpha)  ! antisymmetry
 10 continue
! test ccmat
  write(*,*)
  write(*,*)'test ccmat'
   do alpha=1,Nmany
     do a=1,Nsing
       do b=1,Nsing
       if(ccmat(a,b,alpha).ne.0) write(*,90)a,b,(ioc(s,alpha),s=1,Nsing),ccmat(a,b,alpha)
    90 format(i3,1x,i3,2x,10(i1),2x,i10)
       enddo
     enddo
  enddo

!  matrix ccmat / end

!  matrix cmat / begin

    allocate(cmat(Nsing,Nmany), stat=ierr)

    do 100 alpha=1,Nmany
       do a=1,Nsing
        jcc(a)=ioc(a,alpha)    ! store occupation numbers  jc(is)=ioc(is,im)
       enddo
       do 200 a=1,Nsing
         cmat(a,alpha)=0
         if(jcc(a).eq.0) go to 200
         jcc(a)=0
         call CONV10(-1,n,2,Nsing,jcc)
         cmat(a,alpha)=n+1
         jcc(a)=1    ! restore initial value for next case
         ! calc sign prefactor / begin
         ksgn=0
           if(a.gt.1) then
              do s=1,a-1
                 ksgn=ksgn+jcc(s)
               enddo
            endif
         ! calc sign prefactor / end
       cmat(a,alpha)=(-1)**ksgn*cmat(a,alpha)
    200 continue
 100 continue

! test cmat
  write(*,*)
  write(*,*)'test cmat'
   do alpha=1,Nmany
     do a=1,Nsing
       if(cmat(a,alpha).ne.0) write(*,900)a,(ioc(s,alpha),s=1,Nsing),cmat(a,alpha)
    900 format(i3,1x,2x,10(i1),2x,i10)
     enddo
  enddo


!  matrix cmat / end


  END SUBROUTINE CC_Matrix
!=====================================================================
  SUBROUTINE Coulomb_Integrals
! calculates the matrix of Coulomb integrals VCmat(a,b,c,d)
  use Mod_Many_Particle_States
  implicit none
  integer q1,q2,a,b,c,d,k1,k2,kk
  real(kind=dp) x1,x2,y1,y2,uc12
  integer kont(Nsing,Nsing,Nsing,Nsing)

  VCmat=0 ; if(uc0.eq.0.0_dp.or.Nel.eq.1) return

  kont=0

  do a=1,Nsing; do b=1,Nsing; do c=1,Nsing; do d=1,Nsing 

   write(*,*) a,b
   if(kont(a,b,c,d).eq.1) go to 50   ! use symmetries
      VCmat(a,b,c,d)=(0.0_dp,0.0_dp)

      do q1=1,Nstat; do q2=1,Nstat 
         k1=qn(q1,1); k2=qn(q2,1)
         kk=1; if(k1.eq.k2) kk=0
         x1=pos(k1,3); y1=pos(k1,4)
         x2=pos(k2,3); y2=pos(k2,4)
         uc12=uc0*kk/(sqrt((x1-x2)**2+(y1-y2)**2)+eta)
!        uc12=0; if(k1.ne.k2) uc12=uc0/(sqrt((x1-x2)**2+(y1-y2)**2))  ! version with eta=0
         VCmat(a,b,c,d)=VCmat(a,b,c,d)+uc12*Conjg(Psi(q1,a))*Conjg(Psi(q2,b))*Psi(q1,c)*Psi(q2,d)
      enddo; enddo

      kont(a,b,c,d)=1
      VCmat(b,a,d,c)=VCmat(a,b,c,d); kont(b,a,d,c)=1           ! use symmetries
      VCmat(c,d,a,b)=Conjg(VCmat(a,b,c,d)); kont(c,d,a,b)=1    ! use symmetries
      VCmat(d,c,b,a)=Conjg(VCmat(a,b,c,d)); kont(d,c,b,a)=1    ! use symmetries

  50 continue    
  enddo; enddo; enddo; enddo

!  test symmetry
!10 write(*,*)'a b c d ='
!   read(*,*) a,b,c,d
!   write(*,*) VCmat(a,b,c,d)
!   write(*,*) VCmat(b,a,d,c)
!   write(*,*) Conjg(VCmat(c,d,a,b))
!   write(*,*) Conjg(VCmat(d,c,b,a))
!   go to 10

  END SUBROUTINE Coulomb_Integrals
!=====================================================================
  SUBROUTINE Hamiltonian_Many
! Matrix elements of the many-body Hamiltonian
  use Mod_Many_Particle_States
  implicit none
  integer alpha,beta,a,b,c,d,ierr,ba,dc,s

  allocate(HHmat(Nmany,Nmany), stat=ierr)

  HHmat=(0.0_dp,0.0_dp)

! diagonal term 
  do alpha=1,Nmany
    do a=1,Nsing
      HHmat(alpha,alpha)=HHmat(alpha,alpha)+Ea(a)*ioc(a,alpha)
    enddo
  enddo

! Coulomb term
! do alpha=1,Nmany; do beta=1,Nmany      ! do not use Hermiticity
  do alpha=1,Nmany; do beta=alpha,Nmany  ! use Hermiticity
      do a=1,Nsing; do b=1,Nsing; do c=1,Nsing; do d=1,Nsing
        ba=ccmat(b,a,alpha); dc=ccmat(d,c,beta)
        if(abs(ba).eq.abs(dc).and.ba.ne.0) then 
          s=sign(1,ba)*sign(1,dc)
          HHmat(alpha,beta)=HHmat(alpha,beta)+0.5*VCmat(a,b,c,d)*s 
        endif
      enddo; enddo; enddo; enddo
  HHmat(beta,alpha)=Conjg(HHmat(alpha,beta)) ! use Hermiticity
  enddo; enddo

!  do alpha=1,Nmany; do beta=1,Nmany
!   write(*,*) alpha,beta,HHmat(alpha,beta)
!  enddo; enddo
      
  END SUBROUTINE Hamiltonian_Many
!=====================================================================
  SUBROUTINE Many_Particle_Eigenstates
! eigenstates of the many-body Hamiltonian
  use Mod_Many_Particle_States
!  use mkl95_lapack      
  use lapack95    
!  use libmkl_lapack95_lp64
  implicit none
  real(kind=dp) norm,coef
  integer ierr,mu,alpha,a

  allocate(EEa(Nmany), stat=ierr) 
  allocate(PPsi(Nmany,Nmany), stat=ierr)
  allocate(PPari(Nmany), stat=ierr) 

  PPsi=HHmat
! call heev(a=PPsi, w=EEa, jobz='V', uplo='U', info=ierr)  ! suspicious algorithm?
  call heevd(a=PPsi, w=EEa, jobz='V', uplo='U', info=ierr)  ! divide-and-conquer algorithm
 
! write(*,*) 'ierr=',ierr

! write many-body energyes
  do mu=1,Nmany
    write(*,*) mu,EEa(mu)
    write(13,*) mu,EEa(mu) 
   enddo
  END SUBROUTINE Many_Particle_Eigenstates
!=====================================================================
  SUBROUTINE Charge_Spin_Density 
! calc charge and spin densities for Nmu states
  use Mod_Many_Particle_States
  implicit none
  character*1 type
  integer a,b,q1,q2,k,s1,s2,k1,k2
  integer mu,alpha,beta,va,vb,s,ix,g
  real(kind=dp) r,phi,am1,am2,cm1,cm2,sm1,sm2
  complex(kind=dp), dimension(Nsing,Nsing,Nsites,0:3) :: Tab
  complex(kind=dp) SigMat,term(0:3),trx,fcx,U1mu
  character*2 filenumber

! calculation of Tab /begin
  Tab=0
! do a=1,Nsing; do b=1,Nsing
  do a=1,Nsing; do b=1,a      ! use Hermiticity
    do q1=1,Nstat; do q2=1,Nstat
      if(qn(q1,1).eq.qn(q2,1)) then 
        k=qn(q1,1); s1=qn(q1,2); s2=qn(q2,2)  ! site and spin labels
        Tab(a,b,k,0)=Tab(a,b,k,0)+Conjg(Psi(q1,a))*SigMat('1',s1,s2,0)*Psi(q2,b) ! charge
        Tab(a,b,k,1)=Tab(a,b,k,1)+Conjg(Psi(q1,a))*SigMat('x',s1,s2,0)*Psi(q2,b) ! spin x
        Tab(a,b,k,2)=Tab(a,b,k,2)+Conjg(Psi(q1,a))*SigMat('y',s1,s2,0)*Psi(q2,b) ! spin y
        Tab(a,b,k,3)=Tab(a,b,k,3)+Conjg(Psi(q1,a))*SigMat('z',s1,s2,0)*Psi(q2,b) ! spin z
      endif
      Tab(b,a,k,0)=Conjg(Tab(a,b,k,0))     ! use Hermiticity
      Tab(b,a,k,1)=Conjg(Tab(a,b,k,1))     ! use Hermiticity
      Tab(b,a,k,2)=Conjg(Tab(a,b,k,2))     ! use Hermiticity
      Tab(b,a,k,3)=Conjg(Tab(a,b,k,3))     ! use Hermiticity
    enddo; enddo
  enddo; enddo
! calculation of Tab /end

! test Tab / begin
!70  write(*,*)'site='
!    read(*,*)k
!    write(*,*)'a,b'
!    read(*,*)a,b
!      write(*,*)'0',Tab(a,b,k,0)
!      write(*,*)'0',Tab(b,a,k,0)
!      write(*,*)
!      write(*,*)'x',Tab(a,b,k,1)
!      write(*,*)'x',Tab(b,a,k,1)
!      write(*,*)
!      write(*,*)'y',Tab(a,b,k,2)
!      write(*,*)'y',Tab(b,a,k,2) 
!      write(*,*)
!      write(*,*)'z',Tab(a,b,k,3)
!      write(*,*)'z',Tab(b,a,k,3)
!      write(*,*)
!  go to 70
! test Tab / end

! calculation charge-spin density CSd / begin

  CSd=0

  do mu=1,Ncsd
    do k=1,Nsites
     
      do alpha=1,Nmany; do beta=1,Nmany   ! sum_{alpha,beta} / begin 

        fcx=Conjg(PPsi(alpha,mu))*PPsi(beta,mu)
        term=0

   if(abs(fcx).gt.(1.0d-16)) then 

        do a=1,Nsing; do b=1,Nsing  ! sum_{a,b} / begin
          va=cmat(a,alpha); vb=cmat(b,beta)
          if(abs(va).eq.abs(vb).and.va.ne.0) then 
            s=sign(1,va)*sign(1,vb)
            do ix=0,3; term(ix)=term(ix)+Tab(a,b,k,ix)*s; enddo
          endif
        enddo; enddo             ! sum_{a,b} / end
        
        do ix=0,3   ! add to sum_{alpha,beta} 
!         trx=Conjg(PPsi(alpha,mu))*PPsi(beta,mu)*term(ix)
          trx=fcx*term(ix)
          CSd(mu,k,ix)=CSd(mu,k,ix)+Real(trx)
        enddo

   endif

      enddo; enddo                       ! sum_{alpha,beta} / end

        do ix=1,3   
          CSd(mu,k,ix)=CSd(mu,k,ix)/2  ! S = 1/2 Sigma => spin dens. in units of hbar
        enddo

    enddo  ! k=1,Nsites
   write(*,*) 'calculated CSd mu=',mu
  enddo  ! mu=1,Ncsd

! calculation charge-spin density CSd / end


! calculation charge-spin density considering degeneracies CSdg / begin

  CSdg=0

  do mu=1,Ncsd
    g=Deg_group(mu)
    do k=1,Nsites; do ix=0,3
      CSdg(g,k,ix)=CSdg(g,k,ix)+CSd(mu,k,ix)/Deg_order(g)
    enddo; enddo
  enddo

! calculation charge-spin density considering degeneracies CSdg / end

! write charge-spin density 
  write(*,*) 'Max. group', g,'possibly incomplete'
  do mu=1,Ncsd
       ! open files / begin
         if(mu.lt.10) write(filenumber,10) mu; 10 format(i1)
         if(mu.ge.10) write(filenumber,20) mu; 20 format(i2)
         open(100+mu,file='CSd.otp.'//filenumber,status='unknown')
         open(200+mu,file='CSdg.otp.'//filenumber,status='unknown')
       ! open files / end
    g=Deg_group(mu)
    write(200+mu,*) '# deg. group, order',mu,g,Deg_order(g)
    do k=1,Nsites
    !pos(k,3)=x, pos(k,4)=y
    write(100+mu,300) k,pos(k,3),pos(k,4),(CSd(mu,k,ix),ix=0,3)
    write(200+mu,300) k,pos(k,3),pos(k,4),(CSdg(g,k,ix),ix=0,3)
    if (mod(k,Nphi).eq.0)  then 
        k1=k-Nphi+1
        write(100+mu,300) k1,pos(k1,3),pos(k1,4),(CSd(mu,k1,ix),ix=0,3)
        write(200+mu,300) k1,pos(k1,3),pos(k1,4),(CSdg(g,k1,ix),ix=0,3)
        write(100+mu,*)
        write(200+mu,*)
    endif
    enddo  
!   if(mu.eq.1) write(16,*)'# MES, E, Nel, Sx, Sy, Sz' 
!   write(16,200) mu,EEa(mu),EEa(mu),sum(CSd(mu,1:Nsites,0)),  & ! MES, Energy, Nel
!   sum(CSd(mu,1:Nsites,1)), sum(CSd(mu,1:Nsites,2)),sum(CSd(mu,1:Nsites,3)) ! Sx, Sy, Sz
    close(100+mu); close(200+mu)
  enddo  

!  do mu=1,g
!    do k=1,Nsites
!    !pos(k,3)=x, pos(k,4)=y
!    enddo  
!  enddo  

 300 format(i5,2(2x,f10.6),4(2x,e14.7))
 200 format(i5,2(2x,e14.7),2x,f4.1,3(2x,e14.7))

! calculation amplitude of the CDW / begin

  do mu=1,Ncsd

    am1=0; am2=0; sm1=0; sm2=0;

    do k=1,Nsites
    r=pos(k,1); phi=pos(k,2)  ! polar coordinates
    k1=Nint((Rext-r)/d_r+1); k2=Nint(phi/d_phi)+1  ! polar coord discretized
    if(Nr.eq.1) k1=1

!   write(*,*) 'xxx',k1,r
!   stop

     if(k1.eq.(Nr/2+1)) then 
         am1=am1+CSd(mu,k,0)/Nphi        ! mean value of power 1 charge
         am2=am2+(CSd(mu,k,0))**2/Nphi   ! mean value of power 2 charge
         sm1=sm1+CSd(mu,k,3)/Nphi        ! mean value of power 1 spinz
         sm2=sm2+(CSd(mu,k,3))**2/Nphi   ! mean value of power 2 spinz
         if(k2.eq.(Nphi/8)) cm1=CSd(mu,k,0)  ! value at phi=pi/4
     endif    
    enddo
    
    CDWamp(mu,1)=sqrt(am2-am1**2)  ! std evaluation charge
    CDWamp(mu,2)=(cm1-am1)/am1     ! Nowak Szafran
    CDWamp(mu,3)=sqrt(sm2-sm1**2)  ! std evaluation spinz

  enddo

! calculation amplitude of the CDW / end



  END SUBROUTINE Charge_Spin_Density
!=====================================================================
  SUBROUTINE Info_Many_Particle_Eigenstates
! write information about WF's and many-body states
  use Mod_Many_Particle_States
  implicit none
  integer a,mu,alpha,ierr
  real(kind=dp) norm,coef,tol

  allocate(Deg_group(Nmany), stat=ierr) ! Label for groups of degenerated states 
  allocate(Deg_order(Nmany), stat=ierr) ! Order of degeneracy for each group

! find degeneracies / begin

  tol=0.001_dp 
! tol=0.01_dp 
  Nmanyg=1; Deg_group(1)=1; Deg_order(1)=1 

  do mu=2,Nmany
    if(abs(EEa(mu)-EEa(mu-1)).lt.tol) then
      Deg_group(mu)=Deg_group(mu-1)
      Deg_order(Deg_group(mu))=Deg_order(Deg_group(mu))+1
    else
      Deg_group(mu)=Deg_group(mu-1)+1
      Deg_order(Deg_group(mu))=1
    endif
  enddo
  Nmanyg=Deg_group(Nmany)

! find degeneracies / end

  do mu=1,Nmany
    norm=sqrt(sum(abs(PPsi(1:Nmany,mu))**2))
    write(14,*) ' emag=',emag
    write(14,*) 'MES ',mu,'E ',EEa(mu),'norm ',norm
    if(abs(norm-1.0_dp).gt.0.01) then 
       write(*,*)'WARNING norm of many-particle states not 1 for state', mu
       stop
    endif
    write(14,*)'## degeneracy group, order',Deg_group(mu),Deg_order(Deg_group(mu))
    write(14,*)'## composition'
    norm=0
       do alpha=1,Nmany
         coef=abs(PPsi(alpha,mu))
         if(coef.gt.1e-1) then 
            PPari(mu)=1
            do a=1,Nsing; if(ioc(a,alpha).eq.1) PPari(mu)=PPari(mu)*Pari(a); enddo  ! calc parity 
            write(14,100) alpha,PPsi(alpha,mu),PPari(mu),(ioc(a,alpha),a=1,Nsing)
            norm=norm+coef**2   ! calculated only with displayed coefficients
         endif
       enddo
    write(14,*)'norm=',sqrt(norm)
    write(14,*)
    write(14,*)'====================================================='
  enddo

100 format(2x,i5,2x,'(',f10.7,',',f10.7,')',2x,i2,2x,20(i1))
  END SUBROUTINE Info_Many_Particle_Eigenstates
!=====================================================================
  FUNCTION SigMat(type,s1,s2,phi) ! Pauli matrices
    use Mod_Data
    character*1 type  ! x, y, z, r=radial, p=angular, 1=unity matrix
    integer s1,s2     ! s1,s2 = +/-1
    real(kind=dp) phi
    complex(kind=dp) SigMat

! linear x
    if(type.eq.'x') then
      if(s1.eq.s2) then
        SigMat=(0.0_dp,0.0_dp)       ! s1=s2=1 or s1=s2=-1
      else
        SigMat=(1.0_dp,0.0_dp)       ! s1=-s2=1 or -s1=s2=1
      endif
    return
    endif

! linear y
    if(type.eq.'y') then
      if(s1.eq.s2) then
        SigMat=(0.0_dp,0.0_dp)       ! s1=s2=1 or s1=s2=-1
      else
        SigMat=s2*(0.0_dp,1.0_dp)    ! s1=-s2=1 or -s1=s2=1
      endif
    return
    endif

! linear z
    if(type.eq.'z') then
      if(s1.eq.s2) then
        SigMat=s2*(1.0_dp,0.0_dp)    ! s1=s2=1 or s1=s2=-1
      else
        SigMat=s2*(0.0_dp,0.0_dp)    ! s1=-s2=1 or -s1=s2=1
      endif
    return
    endif

! radial r
    if(type.eq.'r') then
      if(s1.eq.s2) then
        SigMat=(0.0_dp,0.0_dp)       ! s1=s2=1 or s1=s2=-1
      else
        SigMat=exp(s2*phi*ci)        ! s1=-s2=1 or -s1=s2=1
      endif
    return
    endif

! angular p (phi)
    if(type.eq.'p') then
      if(s1.eq.s2) then
        SigMat=(0.0_dp,0.0_dp)       ! s1=s2=1 or s1=s2=-1
      else
        SigMat=s2*ci*exp(s2*phi*ci)  ! s1=-s2=1 or -s1=s2=1
      endif
    return
   endif

! unity matrix
    if(type.eq.'1') then
      if(s1.eq.s2) then
        SigMat=(1.0_dp,0.0_dp)   ! s1=s2 => 1
      else
        SigMat=(0.0_dp,0.0_dp)    
      endif
    return
    endif

    END FUNCTION SigMat
!=====================================================================
  SUBROUTINE CONV10(iopt,n,ib,k,ia)
! converts integer n from basis 10 to basis ib for iopt=+1
! converts integer n from basis ib to basis 10 for iopt=-1
    dimension ia(1)
    if(iopt.eq.1) then
      nx=n
      do 10 i=1,k
      if(nx.ge.ib**(k-i)) go to 20
   10 ia(i)=0
      return
   20 i1=i
      do 30 i=k,i1,-1
      ia(i)=mod(nx,ib)
   30 nx=(nx-ia(i))/ib
      return
    endif
    if(iopt.eq.-1) then
      n=0
      do 100 i=1,k
  100 n=n+ia(i)*ib**(k-i)
      return
    endif
  END SUBROUTINE CONV10
!=====================================================================
  SUBROUTINE Open_Output_Files
  use Mod_Many_Particle_States
    character*2 filenumber

    open(10,file="Sample.otp",status='unknown')
    open(11,file="Energy_single.otp",status='unknown')
    open(12,file="WF_single.otp",status='unknown')
    open(13,file="Energy_many.otp",status='unknown')
    open(14,file="WF_many.otp",status='unknown')
    open(15,file="U1mu.otp",status='unknown')
    open(16,file="CStotal.otp",status='unknown')

    do mu=1,Ncsd
      if(mu.lt.10) write(filenumber,10) mu; 10 format(i1)
      if(mu.ge.10) write(filenumber,20) mu; 20 format(i2)
      j=300+mu
      open(j,file='CDWamp.otp.'//filenumber,status='unknown')
    enddo

    open (20,file="Energy_B.otp",status='unknown')

  END SUBROUTINE Open_Output_files
!==================================================================================
  SUBROUTINE Close_Output_Files
    use Mod_Single_Particle_States
    close(10); close(11); close(12); close(13); close(14); close(15); !close(16)
    do mu=1,Ncsd; close(300+mu); enddo
    
  END SUBROUTINE Close_Output_Files
!==================================================================================
  SUBROUTINE Deallocate_open_matrices
    use Mod_Many_Particle_States
    implicit none
    integer ierr
    deallocate(ioc, stat=ierr)
    deallocate(ccmat, stat=ierr)
    deallocate(cmat, stat=ierr)
    deallocate(HHmat, stat=ierr)
    deallocate(EEa, stat=ierr)
    deallocate(PPsi, stat=ierr)
    deallocate(PPari, stat=ierr)
    deallocate(Deg_group, stat=ierr) 
    deallocate(Deg_order, stat=ierr) 
  END SUBROUTINE Deallocate_open_matrices
!==================================================================================
