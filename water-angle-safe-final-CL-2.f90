PROGRAM ELPTICINE

!This code calculates HB Parameters of 1st and 2nd nearest solvent water molecule of 
! N1 and N2 atoms of ellipticine molecule. Further it calculates other HB parmeters 
! of nearest labelled molecule.

!------------------------------------------------
!  A__C    = Covalent bond
!  A----B  = Heavey atom distances
!  B----C  = HB distances
! // HB ANGLES //
!   alpha1 = /_C
! hbbetaa1 = /_B
! hbgamma2 = /_A
!------------------------------------------------




!------------------------------------------------
!            S_b
!           A___C           
!           |  /          d_(OH)   = S_b
!       S_c | :   S_a     d_(O--H) = S_a
!           |/            d_(O--O) = S_c
!           B     
! S_a-S_b = Proton transfer coordinate
!------------------------------------------------

IMPLICIT NONE

INTEGER          :: q,sol,i,j,k,index,N1,N2,idoi,a,natoms,index1,index2, nsolute
INTEGER          :: ll, jj , NSTEP
INTEGER          :: rcount, ActualSize
INTEGER          :: sindex1, sindex2
INTEGER          :: icount, jcount, jarray, kcount, karray, lcount, larray
REAL*8           :: isum, ksum, sum
REAL*8           :: S,cmx,cmy,cmz,b
REAL*8           :: start,finish
REAL*8           :: alpha1,theta2,hbbetaa1,hbgamma2
REAL*8           :: dist_ref_O, LSI, cutoffr
REAL*8           :: delx, dely, delz
DIMENSION        :: jarray(100), karray(100), larray (100)
CHARACTER(len=30):: junk1 
REAL*8,DIMENSION(:),ALLOCATABLE::rdx,rdy,rdz,rx,ry,rz
CHARACTER(len=10),DIMENSION(:),ALLOCATABLE::ele, elesolv
REAL*8,DIMENSION(:),ALLOCATABLE:: InputData, rref
PARAMETER(nsolute=33, NSTEP=130841)
!PARAMETER(nsolute=33, NSTEP=50000)

!-----------------------------------------------------------------------------
!READ statement
CALL cpu_time(start)

q=1
IF (q==1) THEN 
    sol = 261                        !number of water molecules
      b = 14.861d0                     !box length of W-elp
      k = 3                            !number of atoms in a single molecule
ELSE IF (q==2) THEN
    sol = 405                        !number of methanol molecules
      b = 16.274d0                     !box length of W-elp
      k = 6                            !number of atoms in a single molecule
ELSE 
    sol = 733                        !number of ethylene glycol molecules
      b = 18.964d0                     !box length of W-elp
      k = 10                           !number of atoms in a single molecule
END IF

!----------------------------------------------------------------------------
!OPEN FILEs

OPEN(2,FILE='POSITION.xyz' ,STATUS='unknown',form='formatted')
WRITE(2,*)sol+nsolute
WRITE(2,*)"14.861 14.861 14.861"
CLOSE(2)

!OPEN(4,FILE='trial.xyz' ,STATUS='unknown',form='formatted')
!WRITE(4,*)"nsolute"
!WRITE(4,*)"water-elp"
!CLOSE(4)

!----------------------------------------------------------------------------


!----------------------------------------------------------------------------
!READ FILEs

IF (q==1) THEN 
!OPEN (unit=1, FILE='../mod-output.pos.xyz', STATUS='OLD', action='READ')
OPEN (unit=1, &
     !& FILE='/home/subbu/Downloads/Unmesh-Files/UPDATED_CODES/water/H2O-300K-pos-1.xyz', &
     & FILE='/home/subbu/Downloads/Unmesh-Files/UPDATED_CODES/water/Updated_rev/mod-output.pos.xyz', &
     & STATUS='OLD', action='READ')

ENDIF

PRINT*, "FORT.401", "----", " HB Solvent at outer sphere"
PRINT*, "FORT.301", "----", "ARRAY OF INDEX OF WATER MOLECULES AWAY FROM NITROGEN N2"
PRINT*, "FORT.302", "----", "ARRAY OF INDEX OF WATER MOLECULES AWAY FROM NITROGEN N1"
PRINT*, "FORT.303", "----", "ARRAY OF INDEX OF WATER MOLECULES AWAY FROM BOTH NITROGEN ATOMS"
PRINT*, "---------------------------"
PRINT*, "FORT.101", "----", "1st Nearest water molecule  to Pyrrole Nitrogen atoms"
PRINT*, "FORT.102", "----", "2nd Nearest water molecule  to Pyrrole Nitrogen atoms"
PRINT*, "FORT.103", "----", "2nd Nearest water molecule  to Pyrrole Nitrogen atoms"
PRINT*, "---------------------------"
PRINT*, "FORT.111", "----", "1st Nearest water molecule  to Pyridine Nitrogen atoms"
PRINT*, "FORT.112", "----", "2nd Nearest water molecule  to Pyridine Nitrogen atoms"
PRINT*, "FORT.113", "----", "HB water molecule  Ow1 (Pyrrole) Nitrogen atoms"
PRINT*, "FORT.114", "----", "HB water molecule  Ow1 Pyridine  Nitrogen atoms"
PRINT*, "FORT.115", "----", "HB water molecule  Hw1 (Pyrrole) Nitrogen atoms"
PRINT*, "FORT.116", "----", "HB water molecule  Hw2 (Pyrrole) Nitrogen atoms"
PRINT*, "FORT.117", "----", "HB water molecule  Hw1 Pyridine  Nitrogen atoms"
PRINT*, "FORT.118", "----", "HB water molecule  Hw2 Pyridine  Nitrogen atoms"


OPEN(101,FILE='PTC_NU2.dat' ,STATUS='unknown',form='formatted')
WRITE(101,*)"#PTC at nu2"
CLOSE(101)
OPEN(102,FILE='2PTC_NU2.dat' ,STATUS='unknown',form='formatted')
WRITE(102,*)"#PTC at nu2 second solvent"
CLOSE(102)
OPEN(111,FILE='PTC_NU1.dat' ,STATUS='unknown',form='formatted')
WRITE(111,*)"#PTC at nu1"
CLOSE(111)
OPEN(112,FILE='2PTC_NU1.dat' ,STATUS='unknown',form='formatted')
WRITE(112,*)"#PTC at nu1 second solvent"
CLOSE(112)
!OPEN(113,FILE='2PTC_NU1.dat' ,STATUS='unknown',form='formatted')
!WRITE(113,*)"#PTC at nu1 second solvent"
!CLOSE(113)
OPEN(113,FILE='OOPTC_NU2.dat' ,STATUS='unknown',form='formatted')
WRITE(113,*)"#OOPTC at nu2"
CLOSE(113)
OPEN(114,FILE='OOPTC_NU1.dat' ,STATUS='unknown',form='formatted')
WRITE(114,*)"#OOPTC at nu1"
CLOSE(114)
OPEN(115,FILE='H1OPTC_NU2.dat' ,STATUS='unknown',form='formatted')
WRITE(115,*)"#HOPTC at nu2"
CLOSE(115)
OPEN(116,FILE='H2OPTC_NU2.dat' ,STATUS='unknown',form='formatted')
WRITE(116,*)"#HOPTC at nu2"
CLOSE(116)
OPEN(117,FILE='H1OPTC_NU1.dat' ,STATUS='unknown',form='formatted')
WRITE(117,*)"#HOPTC at nu1"
CLOSE(117)
OPEN(118,FILE='H2OPTC_NU1.dat' ,STATUS='unknown',form='formatted')
WRITE(118,*)"#HOPTC at nu1"
CLOSE(118)
OPEN(401,FILE='Bulk_PTC_NU-OOD.dat' ,STATUS='unknown',form='formatted')
WRITE(401,*)"#Bulk_PTC at nu2"
CLOSE(401)
OPEN(402,FILE='Bulk_PTC_NU-OOA.dat' ,STATUS='unknown',form='formatted')
WRITE(402,*)"#Bulk_PTC at nu2 second solvent"
CLOSE(402)
OPEN(637,FILE='LSI-Info-BulkO.dat' ,STATUS='unknown',form='formatted')
!WRITE(637,*)"#LSI-Info-BulkO"
CLOSE(637)


DO idoi=1,NSTEP
!DO idoi=1,1

!WRITE(*,*)"                               "
!WRITE(*,*)"                               "
!WRITE(*,*)idoi

ALLOCATE(rdx(nsolute))
ALLOCATE(rdy(nsolute))
ALLOCATE(rdz(nsolute))
ALLOCATE(rx(sol))
ALLOCATE(ry(sol))
ALLOCATE(rz(sol))
ALLOCATE(ele(nsolute))
ALLOCATE(elesolv(sol))

READ(1,*)natoms
READ(1,*)junk1
DO i=1,nsolute
  READ(1,*)ele(i),rdx(i),rdy(i),rdz(i)
END DO
DO j=1,sol
  READ(1,*)elesolv(j),rx(j),ry(j),rz(j)
END DO
!---------------------------------------------------------------------------
!calculating distance and neighbours
a=1
DO i=1,nsolute
   IF (ele(i)=='N') THEN
      IF (a==1) N1=i
      IF (a==2) N2=i
      a=2
   END IF
END DO


!----------------------------------------------------------------------------

!CALL pbccom(CMx,CMy,CMz,rdx,rdy,rdz,ele,elesolv,sol)
!------------------------------------------------------------------------------

!CALL rescale(rdx,rdy,rdz,rx,ry,rz,CMx,CMy,CMz,ele,idoi,sol,b)

! -----------------------------------------------------------------------
!OPEN(4,FILE='trial.xyz' ,STATUS='OLD',POSITION='APPEND')

!IF (idoi.ne.1) THEN
!WRITE(4,*)'45'
!WRITE(4,*)'water-elp'
!END IF
!DO i=1,nsolute
!WRITE(4,*)ele(i),rdx(i),rdy(i),rdz(i)
!END DO

CLOSE(4)
! -----------------------------------------------------------------------

!PRINT*,"rescaling with N1", CMx, CMy, CMz
!PRINT*,"rescaling with N1", rx(N1), ry(N1), rz(N1)
!PRINT*,"rescaling with N1", rx(N2), ry(N2), rz(N2)

!------------------------------------------------------------------------------
! rescale with pyrrole nitrogen
CMx=rdx(N1)
CMy=rdy(N1)
CMz=rdz(N1)

!CALL rescale(rdx,rdy,rdz,rx,ry,rz,CMx,CMy,CMz,ele,idoi,sol,b)
CALL whole(rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,idoi)

!PRINT*,"rescaling with N2", CMx, CMy, CMz
!PRINT*,"rescaling with N2", rdx(N1), rdy(N1), rdz(N1)
!PRINT*,"rescaling with N2", rdx(N2), rdy(N2), rdz(N2)

!-------------------------------------------------------------------------------
!pyrrole nitrogen
!WRITE(*,*) "!pyrrole nitrogen"

!CALL pyrrole1(index1,index2,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N1,idoi)

jarray = 0
jcount=0
dist_ref_O = 0.0d0

DO i=1,sol,3
    dist_ref_O = sqrt(rx(i)*rx(i)+ry(i)*ry(i)+rz(i)*rz(i))
    IF(dist_ref_O .gt. 6.0d0 .and. dist_ref_O .lt. real(b/2.d0)) THEN
      jcount=jcount+1
      jarray(jcount)=i
      !WRITE(*,*)"O index", i, dist_ref_O!,sqrt(rx(i)*rx(i)+ry(i)*ry(i)+rz(i)*rz(i))
    END IF
END DO
!WRITE(301,*)"No of water moleucles away from N2", jcount
!WRITE(301,*) (jarray(i), i = 1, jcount)


!index=index1
!ll=index1
!CALL writ(index,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,k)
!CALL angle(index,N1,N2,rdx,rdy,rdz,rx,ry,rz,sol,alpha1,theta2,hbbetaa1,hbgamma2)
!WRITE(*,*)"pyrrole1=",alpha1, hbbetaa1

!index=index2
!jj=index2
!CALL writ(index,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,k)
!CALL angle(index,N1,N2,rdx,rdy,rdz,rx,ry,rz,sol,alpha1,theta2,hbbetaa1,hbgamma2)
!WRITE(*,*)"pyrrole2=",alpha1, hbbetaa1 


ll = index1
icount=1
!CALL hbsolventdonor1(ll,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,icount,idoi)
ll = index1+1  ! With resepct to Hw1
!CALL hbsolventacceptor1(ll,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,icount,idoi)


!--------------------------------------------------------------------------------


!WRITE(*,*) "!pyridine nitrogen"
!------------------------------------------------------------------------------
! rescale with pyridine nitrogen
CMx=rdx(N2)
CMy=rdy(N2)
CMz=rdz(N2)
!CALL rescale(rdx,rdy,rdz,rx,ry,rz,CMx,CMy,CMz,ele,idoi,sol,b)

!PRINT*,"rescaling with N2", CMx, CMy, CMz
!PRINT*,"rescaling with N2", rdx(N1), rdy(N1), rdz(N1)
!PRINT*,"rescaling with N2", rdx(N2), rdy(N2), rdz(N2)

!CALL pyridine1(sindex1,sindex2,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,idoi)!count) 
karray = 0
dist_ref_O = 0.0d0
kcount=0
DO i=1,sol,3
dist_ref_O = sqrt(rx(i)*rx(i)+ry(i)*ry(i)+rz(i)*rz(i))
!IF(dist_ref_O .gt. 8.0d0 .and. dist_ref_O .lt. real(b/2.d0)) THEN
IF(dist_ref_O .gt. 6.5d0 ) THEN   ! Tryin to get Oxygen atoms which are 
!way from Pyridine N atom
kcount=kcount+1
karray(kcount)=i
!WRITE(*,*)"O index", i, dist_ref_O!,sqrt(rx(i)*rx(i)+ry(i)*ry(i)+rz(i)*rz(i))
END IF
END DO
!WRITE(302,*)"No of water moleucles away from N1 ", kcount
!WRITE(302,*) (karray(i), i = 1, kcount)

!print*,jcount, kcount
larray = 0
lcount = 0
DO i = 1 , kcount     ! array of pyridine nitrogen atoms
   DO j = 1, jcount   ! array of pyrrole nitroge atoms
     IF (karray(i) == jarray(j)) THEN
         lcount = lcount +1
         larray(lcount) = karray(i)
     END IF
   END DO
END DO
WRITE(303,*) (larray(i)+32, i = 1, lcount)



!index=index1
!ll = index1
!CALL writ(index,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,k)
!CALL angle(index,N1,N2,rdx,rdy,rdz,rx,ry,rz,sol,alpha1,theta2,hbbetaa1,hbgamma2)
!WRITE(*,*)"pyridine1=",theta2 , hbgamma2   

!index=index2
!jj = index2
!CALL writ(index,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,k)
!CALL angle(index,N1,N2,rdx,rdy,rdz,rx,ry,rz,sol,alpha1,theta2,hbbetaa1,hbgamma2)
!WRITE(*,*)"pyridine 2nd solvent",theta2, hbgamma2

ll = sindex1
icount=2
!CALL hbsolventdonor1(ll,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,icount,idoi)
ll = sindex1+1
!CALL hbsolventacceptor1(ll,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,icount,idoi)

!-------------------------------------------------------------------------------

!CALL whole(rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,idoi)





! Rescaling with resepct to Oxygen which are away from Nitrogen atoms.
dist_ref_O = 0.0d0
icount = 0
ksum = 0.0d0
cutoffr = 3.7d0

DO i = 1, lcount
      CMx = rx(larray(i))
      CMy = ry(larray(i))
      CMz = rz(larray(i))
      CALL rescalerefo(rx,ry,rz,CMx,CMy,CMz)
!PRINT*, rx(larray(i)), ry(larray(i)), rz(larray(i)) 
!      icount = 3

      ll = larray(i)
!      PRINT*, ll

!---------------------------------------------------------------------------
      !CALL hbsolventdonor1(ll,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,3,idoi)
!---------------------------------------------------------------------------
       jcount=0
   delx = 0.0d0
   dely = 0.0d0
   delz = 0.0d0
   DO j = 1, sol, 3  ! What it stores
     IF (j .ne. ll) THEN
        delx=rx(j)
        dely=ry(j)
        delz=rz(j)
        dist_ref_O = sqrt(delx*delx+dely*dely+delz*delz)
      IF(dist_ref_O .lt. cutoffr ) THEN
        jcount=jcount+1
      END IF
     END IF
   END DO
!print*, jcount
IF (jcount .LE. 1 ) THEN
!PRINT*, "HI"
ELSE 
   ALLOCATE(rref(jcount))
   rref(jcount) = 0.0d0
   jcount = 0
      DO j = 1, sol, 3
       IF (j .ne. ll) THEN
        delx=rx(j)
        dely=ry(j)
        delz=rz(j)
        dist_ref_O = sqrt(delx*delx+dely*dely+delz*delz)
         IF(dist_ref_O .lt. cutoffr ) THEN
           jcount=jcount+1
           rref(jcount)=dist_ref_O 
         END IF
       END IF
    END DO
!   PRINT*, rref
!
   CALL  Sort(rref(1:jcount), jcount)

   sum = 0.0d0 

  DO j = 1, jcount-1
    dist_ref_O = rref(j+1)-rref(j)
    sum = sum + dist_ref_O 
  END DO

   sum = sum/(jcount-1)

   isum = 0.0d0 
   delx = 0.0d0
   dely = 0.0d0
   delz = 0.0d0

   DO j = 1, jcount-1
       dist_ref_O = rref(j+1)-rref(j)
       isum = isum + (dist_ref_O-sum)*(dist_ref_O-sum)
   END DO

!
 DEALLOCATE(rref)

   icount=0
    IF (jcount > 0) THEN
        LSI=isum/(jcount-1)
      IF (LSI == 0.0d0) THEN
        icount=icount+1
      ELSE
!        rcount=rcount+1
!     OPEN(637,FILE='LSI-Info-BulkO.dat' ,STATUS='old',POSITION='APPEND')
!        WRITE(637,*) LSI 
!     CLOSE(637)
      END IF
    END IF

 END IF
!
END DO






!---------------------------------------------------------------------------
DEALLOCATE(rdx)
DEALLOCATE(rdy)
DEALLOCATE(rdz)
DEALLOCATE(rx)
DEALLOCATE(ry)
DEALLOCATE(rz)
DEALLOCATE(ele)
DEALLOCATE(elesolv)

END DO
CLOSE(1)

CALL cpu_time(finish)
!WRITE(*,*)"      "
WRITE(*,*)"time taken by cpu in s",finish-start































!----------------------------------------------------------------------------
!SUBROUTINES

CONTAINS

!-------------------------------------
!finding Centre of Mass

SUBROUTINE pbccom(CMx,CMy,CMz,rdx,rdy,rdz,ele,elesolv,sol) 
IMPLICIT NONE
INTEGER nsolute
INTEGER :: N,C,H,O
INTEGER :: i,j,sol
PARAMETER(nsolute=33)
REAL*8  :: cmx,cmy,cmz
REAL*8  :: sumx,sumy,sumz
REAL*8  :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
CHARACTER(len=10)::ele(nsolute),elesolv(sol)

C=0
N=0
H=0
O=0

DO i=1,nsolute
IF (ele(i)=='N') THEN
sumx=sumx+rdx(i)*14.0067d0
sumy=sumy+rdy(i)*14.0067d0
sumz=sumz+rdz(i)*14.0067d0
N=N+1
ELSEIF (ele(i)=='C') THEN
sumx=sumx+rdx(i)*12.0107d0
sumy=sumy+rdy(i)*12.0107d0
sumz=sumz+rdz(i)*12.0107d0
C=C+1
ELSEIF (ele(i)=='H') THEN
sumx=sumx+rdx(i)*1.0079d0
sumy=sumy+rdy(i)*1.0079d0
sumz=sumz+rdz(i)*1.0079d0
H=H+1
END IF
END DO

C=real(C)
N=real(N)
O=real(O)
H=real(H)

CMx=sumx/real(N*14.0067d0+C*12.0107d0+1.0079d0*H)!+16.d0*O)
CMy=sumy/real(N*14.0067d0+C*12.0107d0+1.0079d0*H)!+16.d0*O)
CMZ=sumz/real(N*14.0067d0+C*12.0107d0+1.0079d0*H)!+16.d0*O)



END SUBROUTINE pbccom







!---------------------------------------------------------
!rescaling according to the reference point CMx, CMy, CMz and APPlying periodic boundary

SUBROUTINE rescale(rdx,rdy,rdz,rx,ry,rz,CMx,CMy,CMz,ele,idoi,sol,b)
IMPLICIT NONE
INTEGER nsolute
PARAMETER(nsolute=33)
REAL*8    :: cmx,cmy,cmz
REAL*8    :: sumx,sumy,sumz
REAL*8    :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8    :: rx(sol),ry(sol),rz(sol),b
REAL*8    :: dist_ref_O
INTEGER   :: i,j,idoi,sol, count , array(100)
CHARACTER(len=10)::ele(nsolute)


DO i=1,nsolute
rdx(i)=rdx(i)-CMx
rdx(i)=rdx(i)-b*anint(rdx(i)/b)
rdy(i)=rdy(i)-CMy
rdy(i)=rdy(i)-b*anint(rdy(i)/b)
rdz(i)=rdz(i)-CMz
rdz(i)=rdz(i)-b*anint(rdz(i)/b)
END DO

DO i=1,sol
rx(i)=rx(i)-CMx
rx(i)=rx(i)-b*anint(rx(i)/b)
ry(i)=ry(i)-CMy
ry(i)=ry(i)-b*anint(ry(i)/b)
rz(i)=rz(i)-CMz
rz(i)=rz(i)-b*anint(rz(i)/b)
END DO

END SUBROUTINE rescale
!----------------------------------------------------------------------------








SUBROUTINE pyrrole1(index1,index2,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N1,idoi)

!
!            S_b
!           N___H    ----- Pyrrole nitrogen
!           |  /          d_(NH)    = S_b
!       S_c | :   S_a     d_(H--Od) = S_a
!           |/            d_(N--Od) = S_c
!           Od
!   alpha1 = /_C
! hbbetaa1 = /_B
! hbgamma2 = /_A

IMPLICIT NONE
INTEGER nsolute,idoi
PARAMETER(nsolute=33)
REAL*8  :: delx,dely,delz,S
REAL*8  :: rx(sol),ry(sol),rz(sol)
REAL*8  :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8  :: min1,min2,miz1,miz2
REAL*8  :: alpha1,hbbetaa1,hbgamma2
REAL*8  :: NH, NO
REAL*8  :: Sa, Sb, Sc
REAL*8  :: Sa1, Sb1, Sc1
REAL*8  :: alpha11
INTEGER :: i,j,sol,n1,index1,index2
CHARACTER(len=10)   ::  ele(nsolute),elesolv(sol)
REAL*8,DIMENSION(:),ALLOCATABLE :: R

ALLOCATE(R(sol))

DO j=1,sol,3
delx=rx(j)-rdx(N1+1)
dely=ry(j)-rdy(N1+1)
delz=rz(j)-rdz(N1+1)
R(j)=sqrt(delx**2+dely**2+delz**2)
END DO

miz1=R(1)
index1=1
DO i=4,sol,3
IF (R(i)<miz1) THEN
miz1=R(i)
index1=i
END IF
END DO


NH=sqrt(rdx(N1+1)*rdx(N1+1)+rdy(N1+1)*rdy(N1+1)+rdz(N1+1)*rdz(N1+1))
NO=sqrt(rx(index1)*rx(index1)+ry(index1)*ry(index1)+rz(index1)*rz(index1))

Sa = miz1  ! H--Od : Distance between pyrrole and Ow of nearest water
Sb = NH
Sc = NO
!CALL HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)

Sa1 = Sa
Sb1 = Sb
Sc1 = Sc
alpha11 = alpha1

!WRITE(101,*) index1+nsolute,Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2)
!WRITE(101,*) index1+nsolute,Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
!IF (hbgamma2 < 35.0d0 ) THEN
!WRITE(101,*) doi, Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
WRITE(201,*) idoi,Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(101,FILE='PTC_NU2.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(101,*) idoi,Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(101)
!END IF

R(index1)=1000

!IF(miz1==R(1)) THEN
!index2=1
!goto 15
!END IF

miz2=R(1)
index2=1


DO i=4,sol,3
!IF (miz1==R(i))THEN
!index2=i
!goto 15
!END IF

IF (R(i) < miz2) THEN
miz2=R(i)
index2=i
END IF
END DO

NO=sqrt(rx(index2)*rx(index2)+ry(index2)*ry(index2)+rz(index2)*rz(index2))
Sa = R(index2) !miz1
DEALLOCATE(R)
Sb = NH
Sc = NO
!CALL HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)

!WRITE(*,*) "pyrrole1",index1+nsolute,Sa,Sb,Sc,alpha1,hbgamma2

!15 WRITE(*,*) "pyrrole2",index2+nsolute,R(index2), NH, NO
!15 WRITE(102,*) "pyrrole2",index2+nsolute, Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb
!IF (hbgamma2 < 35.0d0 ) THEN
WRITE(202,*) idoi, Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(102,FILE='2PTC_NU2.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(102,*) idoi,Sa,Sb,Sc,alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(102)

!END IF

IF (alpha11 < alpha1 .and. Sc > Sc1) THEN
 WRITE(103,*) idoi, Sa,Sb,Sc,alpha1, Sa-Sb! Sb*cos(hbgamma2), cos(hbgamma2)
ELSE
 WRITE(103,*) idoi, Sa1,Sb1,Sc1,alpha11, Sa1-Sb1! Sb*cos(hbgamma2), cos(hbgamma2)
END IF


END SUBROUTINE pyrrole1

!--------------------------------------------------------------------------------------

!------------------------------------------------------------------------------
SUBROUTINE pyridine1(index1,index2,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,idoi)

!
!            S_b
!           Od__Hd           
!           |  /          d_(OdHd)  = S_b
!       S_c | :   S_a     d_(N--Hd) = S_a
!           |/            d_(N--Od) = S_c
!           N     ----- Pyridine nitrogen
!   alpha1 = /_C
! hbbetaa1 = /_B
! hbgamma2 = /_A


IMPLICIT NONE
INTEGER nsolute,idoi
PARAMETER(nsolute=33)
REAL*8  :: delx,dely,delz,S
REAL*8  :: rx(sol),ry(sol),rz(sol)
REAL*8  :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8  :: min1,min2,miz1,miz2
REAL*8  :: OH, NH
REAL*8  :: Sa, Sb, Sc
REAL*8  :: alpha1,hbbetaa1,hbgamma2
INTEGER :: i,j,sol,n2,index1,index2,index,ll
CHARACTER(len=10)::ele(nsolute),elesolv(sol)

REAL*8,DIMENSION(:),ALLOCATABLE :: R

ALLOCATE(R(sol))


DO j=1,sol
delx=rx(j)!-rdx(N2)
dely=ry(j)!-rdy(N2)
delz=rz(j)!-rdz(N2)

R(j)=sqrt(delx**2+dely**2+delz**2)

END DO
!--------------------------

miz1=R(1)
index1=1
DO i=2,sol,1
IF (R(i) < miz1) THEN
miz1=R(i)
index1=i
END IF
END DO

ll=index1

IF (elesolv(index1)=='O') THEN 
index1=index1
ELSEIF(elesolv(index1-1)=='O') THEN
index1=index1-1
ELSEIF(elesolv(index1-2)=='O') THEN
index1=index1-2
END IF
!NH=sqrt(rx(ll)*rx(ll)+ry(ll)*ry(ll)+rz(ll)*rz(ll))
OH = sqrt((rx(index1)-rx(ll))**2+(ry(index1)-ry(ll))**2+(rz(index1)-rz(ll))**2)
Sc = R(index1) ! N----O distance
Sb = OH        ! O-H Covalnet distance
Sa = R(ll)     ! N----H distance
!CALL HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)
!WRITE(111,*) "pyridine1",ll+nsolute, Sb, Sa, Sc, alpha1,hbgamma2, Sa-Sb
!WRITE(111,*) ll+nsolute, Sb, Sa, Sc, alpha1,hbgamma2, Sa-Sb
!IF (hbgamma2 < 35.0d0 ) THEN
WRITE(211,*) idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(111,FILE='PTC_NU1.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(111,*) idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(111)

!END IF
R(index1)=1000
R(index1+1)=1000
R(index1+2)=1000



!IF(miz1==R(1)) THEN
!index2=1
!goto 25
!END IF

miz2=R(1)
index2=1
DO i=2,sol,1
!IF(miz1==R(i)) THEN
!index2=i
!goto 25
!END IF
IF (R(i) < miz2) THEN
miz2=R(i)
index2=i
END IF
END DO

25 ll=index2

IF (elesolv(index2)=='O') THEN
index2=index2
ELSEIF(elesolv(index2-1)=='O') THEN
index2=index2-1
ELSEIF(elesolv(index2-2)=='O') THEN
index2=index2-2
END IF

!NH=sqrt(rx(ll)*rx(ll)+ry(ll)*ry(ll)+rz(ll)*rz(ll))
OH = sqrt((rx(index2)-rx(ll))**2+(ry(index2)-ry(ll))**2+(rz(index2)-rz(ll))**2)
Sc = R(index2) ! N----O distance
Sb = OH        ! O-H Covalnet distance
Sa = R(ll)     ! N----H distance
DEALLOCATE(R)
!CALL HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)
!WRITE(112,*) "pyridine 2nd solvent",ll+nsolute, Sb, Sa, Sc, alpha1,hbgamma2, Sa-Sb
!IF (hbgamma2 < 35.0d0 ) THEN
WRITE(212,*) idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(112,FILE='2PTC_NU1.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(112,*) idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(112)

!END IF
!WRITE(112,*) ll+nsolute, Sb, Sa, Sc, alpha1,hbgamma2, Sa-Sb
!WRITE(*,*) "pyridine 2nd solvent",ll+nsolute,R(ll), R(index2),NH



END SUBROUTINE pyridine1


!-------------------------------------------------------------------------------
SUBROUTINE hbsolventdonor1(index1,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2,icount,idoi)
!
!            S_b
!           O___H           
!           |  /          d_(OH)   = S_b
!       S_c | :   S_a     d_(O--H) = S_a
!           |/            d_(O--O) = S_c
!           Od     
!   alpha1 = /_C
! hbbetaa1 = /_B
! hbgamma2 = /_A
!
!    
IMPLICIT NONE
INTEGER nsolute,idoi
PARAMETER(nsolute=33)
REAL*8  :: delx,dely,delz,S
REAL*8  :: rx(sol),ry(sol),rz(sol)
REAL*8  :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8  :: min1,min2,miz1,miz2
REAL*8  :: NH
REAL*8  :: OOd, HOd, HOd1, HOd2, OH
REAL*8  :: del1hx,del1hy,del1hz,R1h(sol)
REAL*8  :: del2hx,del2hy,del2hz,R2h(sol)
REAL*8  :: alpha1,hbbetaa1,hbgamma2
REAL*8  :: Sa, Sb, Sc
INTEGER::i,j,sol,n2,index1,index2,index,ll, k, index3, icount, initialindex
CHARACTER(len=10)::ele(nsolute),elesolv(sol)

REAL*8,DIMENSION(:),ALLOCATABLE :: R

ALLOCATE(R(sol))


initialindex = index1
DO j=1,sol

    delx=rx(j)-rx(index1)
    dely=ry(j)-ry(index1)
    delz=rz(j)-rz(index1)

R(j)=sqrt(delx**2+dely**2+delz**2)

END DO

!--------------------------


R(index1) = 1000.d0
R(index1+1) = 1000.d0
R(index1+2) = 1000.d0

IF(index1==1) THEN 
k=4
   ELSE
k=1
END IF
miz1=R(k)
index3=1

DO i=k,sol,1
IF (R(i) < miz1) THEN
miz1=R(i)
index3=i
END IF
END DO
DEALLOCATE(R)
ll=index3

index=index1

IF (elesolv(index3)=='O') THEN 
!      print*,"H1"
        OOd = sqrt((rx(index3)-rx(index))**2 + (ry(index3)-ry(index))**2 + (rz(index3)-rz(index))**2)
       HOd1 = sqrt((rx(index3)-rx(index+1))**2 + (ry(index3)-ry(index+1))**2 + (rz(index3)-rz(index+1))**2)
       HOd2 = sqrt((rx(index3)-rx(index+2))**2 + (ry(index3)-ry(index+2))**2 + (rz(index3)-rz(index+2))**2)
    IF ( HOd1 < HOd2) THEN
       HOd = HOd1
       OH = sqrt((rx(index)-rx(index+1))**2 + (ry(index)-ry(index+1))**2 + (rz(index)-rz(index+1))**2)
      ELSE
       HOd = HOd2
       OH = sqrt((rx(index)-rx(index+2))**2 + (ry(index)-ry(index+2))**2 + (rz(index)-rz(index+2))**2)
    END IF
       index3=index3
ELSEIF(elesolv(index3-1)=='O') THEN
!      print*,"H2"
       HOd = sqrt((rx(index3)-rx(index))**2 + (ry(index3)-ry(index))**2 + (rz(index3)-rz(index))**2)
 index3=index3-1
       OOd = sqrt((rx(index3)-rx(index))**2 + (ry(index3)-ry(index))**2 + (rz(index3)-rz(index))**2)
       OH = sqrt((rx(ll)-rx(index3))**2 + (ry(ll)-ry(index3))**2 + (rz(ll)-rz(index3))**2)
ELSE!IF(elesolv(index3-2)=='O') THEN
     ! print*,"H3"
       HOd = sqrt((rx(index3)-rx(index))**2 + (ry(index3)-ry(index))**2 + (rz(index3)-rz(index))**2)
 index3=index3-2
       OOd = sqrt((rx(index3)-rx(index))**2 + (ry(index3)-ry(index))**2 + (rz(index3)-rz(index))**2)
       OH = sqrt((rx(ll)-rx(index3))**2 + (ry(ll)-ry(index3))**2 + (rz(ll)-rz(index3))**2)
END IF


Sa = HOd
Sb = OH
Sc = OOd

!CALL HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)
!OH=sqrt(rx(ll)*rx(ll)+ry(ll)*ry(ll)+rz(ll)*rz(ll))
IF (icount == 1) THEN
!WRITE(113,*) "HBSolvent at N1",ll+nsolute, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb
!WRITE(113,*) ll+nsolute, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb,  Sb*cos(hbgamma2), cos(hbgamma2)
WRITE(213,*) idoi, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb,  Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(113,FILE='OOPTC_NU2.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(113,*) idoi, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb,  Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE (113)
ELSE IF (icount == 2) THEN
!WRITE(114,*) "HBSolvent at N2",ll+nsolute, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb
!WRITE(114,*) ll+nsolute, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
WRITE(214,*) idoi, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(114,FILE='OOPTC_NU1.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(114,*) idoi, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb,  Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE (114)
ELSE IF (icount == 3) THEN
!WRITE(401,*) "HBSolvent at outer sphere",initialindex + nsolute, ll+nsolute, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb
!WRITE(401,*) initialindex + nsolute, ll+nsolute, Sa, Sb, Sc, alpha1, hbbetaa1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
WRITE(501,*) idoi, Sa, Sb, Sc, alpha1, hbbetaa1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(401,FILE='Bulk_PTC_NU-OOD.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
& WRITE(401,*) idoi, Sa, Sb, Sc, alpha1,hbbetaa1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(401)

ENDIF



END SUBROUTINE hbsolventdonor1
!-------------------------------------------------------------------------------








!-------------------------------------------------------------------------------
SUBROUTINE hbsolventacceptor1(index1,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,N2, icount,idoi)
!
!            S_b
!           O___H           
!           |  /          d_(OH)   = S_b
!       S_c | :   S_a     d_(O--H) = S_a
!           |/            d_(O--O) = S_c
!           Od     
!   alpha1 = /_C
! hbbetaa1 = /_B
! hbgamma2 = /_A
!
!    
!

IMPLICIT NONE
INTEGER nsolute,idoi
PARAMETER(nsolute=33)
!REAL*8::delx,dely,delz,R(sol),rx(sol),ry(sol),rz(sol),rdx(nsolute),rdy(nsolute),rdz(nsolute),min1,min2,miz1,miz2
REAL*8  :: delx,dely,delz,S
REAL*8  :: rx(sol),ry(sol),rz(sol)
REAL*8  :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8  :: min1,min2,miz1,miz2
REAL*8  :: NH
REAL*8  :: OOd, HOd, HOd1, HOd2, OH
REAL*8  :: del1hx,del1hy,del1hz,R1h(sol)
REAL*8  :: del2hx,del2hy,del2hz,R2h(sol)
REAL*8  :: alpha1,hbbetaa1,hbgamma2
REAL*8  :: Sa, Sb, Sc
INTEGER::i,j,sol,n2,index1,index2,index,ll, k, index3, COUNT
INTEGER :: icount
CHARACTER(len=10)::ele(nsolute),elesolv(sol)

REAL*8,DIMENSION(:),ALLOCATABLE :: R

ALLOCATE(R(sol))

COUNT = 1

90 IF (COUNT .ne. 1) THEN
index1 = index1+1
index=index1-2
ELSE
index1 = index1
index=index1-1
END IF
DO j=1,sol

    delx=rx(j)-rx(index1)
    dely=ry(j)-ry(index1)
    delz=rz(j)-rz(index1)

R(j)=sqrt(delx**2+dely**2+delz**2)

END DO

!--------------------------

IF (COUNT .ne. 1) THEN

!print*, 't1', R(index1-2), R(index1-1), R(index1+0)
R(index1-2) = 1000.d0
R(index1-1) = 1000.d0
R(index1+0) = 1000.d0
     IF(index1-1==1) THEN 
       k=4
      ELSE
       k=1
     END IF

ELSE

!print*, 't2', R(index1-1), R(index1+0), R(index1+1)
R(index1-1) = 1000.d0
R(index1+0) = 1000.d0
R(index1+1) = 1000.d0
     IF(index1-2==1) THEN 
       k=4
     ELSE
       k=1
     END IF

END IF


miz1=R(k)
index3=1

DO i=k,sol,3
     IF (R(i) < miz1) THEN
         miz1=R(i)
         index3=i
     END IF
END DO


ll=index3


       OOd = sqrt((rx(index3)-rx(index))**2 + (ry(index3)-ry(index))**2 + (rz(index3)-rz(index))**2)
       HOd = sqrt((rx(index3)-rx(index+1))**2 + (ry(index3)-ry(index+1))**2 + (rz(index3)-rz(index+1))**2)
       OH = sqrt((rx(index)-rx(index+1))**2 + (ry(index)-ry(index+1))**2 + (rz(index)-rz(index+1))**2)


Sa = HOd
Sb = OH
Sc = OOd

!CALL HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)
IF (icount == 1) THEN
IF (COUNT==1) THEN
   !WRITE(115,*) "HBSolvent acceptor  Hw1 at N1 ", COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb
   !WRITE(115,*)  COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb
!H1OPTC_NU2.dat
   WRITE(215,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(115,FILE='H1OPTC_NU2.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
&   WRITE(115,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(115)
ELSE
   !WRITE(116,*) "HBSolvent acceptor Hw2  at N1", COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb
   !WRITE(116,*)  COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb
   WRITE(216,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(116,FILE='H2OPTC_NU2.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
&   WRITE(116,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(116)
END IF
ELSE
IF (COUNT==1) THEN
   !WRITE(117,*) "HBSolvent acceptor Hw1  at N2", COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb
   !WRITE(117,*)  COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
   WRITE(217,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(117,FILE='H1OPTC_NU1.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
&   WRITE(117,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(117)
ELSE
   !WRITE(118,*) "HBSolvent acceptor Hw2  at N2", COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb
   !WRITE(118,*)  COUNT, ll+nsolute, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
   WRITE(218,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
OPEN(118,FILE='H2OPTC_NU1.dat' ,STATUS='old',POSITION='APPEND')
IF (Sc < 3.9d0 .and. Sa < 2.5d0 .and. alpha1 > 135.0d0 ) &
&   WRITE(118,*)   idoi, Sa, Sb, Sc, alpha1,hbgamma2, Sa-Sb, Sb*cos(hbgamma2), cos(hbgamma2)
CLOSE(118)
END IF
END IF

COUNT = COUNT+1

IF ( COUNT > 2 ) THEN
GO TO 99
ELSE
GO TO 90
ENDIF 
DEALLOCATE(R)

99 END SUBROUTINE hbsolventacceptor1
!-------------------------------------------------------------------------------










!--------------------------------------------------------------------------------------------
SUBROUTINE whole(rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,idoi)
IMPLICIT NONE
INTEGER nsolute
PARAMETER(nsolute=33)
REAL*8   :: rx(sol),ry(sol),rz(sol)
REAL*8   :: rdx(nsolute),rdy(nsolute),rdz(nsolute)
INTEGER  :: i,j,sol,idoi
CHARACTER(len=10)::ele(nsolute),elesolv(sol)

OPEN(2,FILE='POSITION.xyz' ,STATUS='OLD',POSITION='APPEND')

IF (idoi.ne.1) THEN
WRITE(2,*)sol+nsolute
WRITE(2,*)'14.861 14.861 14.861'
END IF
DO i=1,nsolute
WRITE(2,*)ele(i),rdx(i),rdy(i),rdz(i)
END DO
DO i=1,sol
WRITE(2,*)elesolv(i),rx(i),ry(i),rz(i)
END DO
CLOSE(2)



END SUBROUTINE whole




!------------------------------------------------------------------------------------------------

SUBROUTINE writ(index,rdx,rdy,rdz,rx,ry,rz,ele,elesolv,sol,k)
IMPLICIT NONE
INTEGER nsolute
PARAMETER(nsolute=33)
REAL*8::rx(sol),ry(sol),rz(sol),rdx(nsolute),rdy(nsolute),rdz(nsolute)
INTEGER::i,j,sol,index,k
CHARACTER(len=10)::ele(nsolute),elesolv(sol)


OPEN(4,FILE='trial.xyz' ,STATUS='OLD',POSITION='APPEND')

DO j=1,k
 WRITE(4,*)elesolv(index+j-1),rx(index+j-1),ry(index+j-1),rz(index+j-1)
END DO
CLOSE(4)
END SUBROUTINE writ







!------------------------------------------------------------------------------------------------
SUBROUTINE angle(index,N1,N2,rdx,rdy,rdz,rx,ry,rz,sol,alpha1,theta2,hbbetaa1,hbgamma2)
IMPLICIT NONE
INTEGER nsolute
PARAMETER(nsolute=33)
REAL*8::rx(sol),ry(sol),rz(sol),rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8::Sa,Sb,Sc,alpha1,theta2,hbbetaa1,hbgamma2,iii,jjj,pi
INTEGER::i,j,sol,index,k,n1,n2,ll
CHARACTER(len=10)::ele(nsolute),elesolv(sol)

pi=4*atan(1.d0)


Sa=sqrt((rdx(N1)-rdx(N1+1))**2+(rdy(N1)-rdy(N1+1))**2+(rdz(N1)-rdz(N1+1))**2)
Sb=sqrt((rdx(N1+1)-rx(index))**2+(rdy(N1+1)-ry(index))**2+(rdz(N1+1)-rz(index))**2)
Sc=sqrt((rdx(N1)-rx(index))**2+(rdy(N1)-ry(index))**2+(rdz(N1)-rz(index))**2)
alpha1=(acos((Sa**2+Sb**2-Sc**2)/(2*Sa*Sb)))*(180.d0/pi)
hbbetaa1=(acos((Sc**2+Sa**2-Sb**2)/(2*Sa*Sc)))*(180.d0/pi)


iii=sqrt((rdx(N2)-rx(index+1))**2+(rdy(N2)-ry(index+1))**2+(rdz(N2)-rz(index+1))**2)
jjj=sqrt((rdx(N2)-rx(index+2))**2+(rdy(N2)-ry(index+2))**2+(rdz(N2)-rz(index+2))**2)

IF (iii <= jjj)  THEN
ll=index+1
ELSE 
ll=index+2
END IF 

Sa=sqrt((rdx(N2)-rx(ll))**2+(rdy(N2)-ry(ll))**2+(rdz(N2)-rz(ll))**2)
Sb=sqrt((rx(ll)-rx(index))**2+(ry(ll)-ry(index))**2+(rz(ll)-rz(index))**2)
Sc=sqrt((rdx(N2)-rx(index))**2+(rdy(N2)-ry(index))**2+(rdz(N2)-rz(index))**2)
theta2=(acos((Sa**2+Sb**2-Sc**2)/(2*Sa*Sb)))*(180.d0/pi)
hbgamma2=(acos((Sc**2+Sb**2-Sa**2)/(2*Sc*Sb)))*(180.d0/pi)

END SUBROUTINE angle


!------------------------------------------------------------------------------------------------
SUBROUTINE HBANGLE(Sa,Sb,Sc,alpha1,hbbetaa1,hbgamma2,idoi)


!
!            S_b
!           O___H           
!           |  /          d_(OH)   = S_b
!       S_c | :   S_a     d_(O--H) = S_a
!           |/            d_(O--O) = S_c
!           Od     
!   alpha1 = /_C
! hbbetaa1 = /_B
! hbgamma2 = /_A



IMPLICIT NONE
INTEGER nsolute, idoi
!PARAMETER(nsolute=33)
!REAL*8::rx(sol),ry(sol),rz(sol),rdx(nsolute),rdy(nsolute),rdz(nsolute)
REAL*8::Sa,Sb,Sc,alpha1,theta2,hbbetaa1,hbgamma2,iii,jjj,pi
INTEGER::i,j,sol,index,k,n1,n2,ll
!CHARACTER(len=10)::ele(nsolute),elesolv(sol)
pi=4*atan(1.d0)

alpha1=(acos((Sa**2+Sb**2-Sc**2)/(2*Sa*Sb)))*(180.d0/pi)
hbbetaa1=(acos((Sc**2+Sa**2-Sb**2)/(2*Sa*Sc)))*(180.d0/pi)
hbgamma2=(acos((Sc**2+Sb**2-Sa**2)/(2*Sb*Sc)))*(180.d0/pi)

END SUBROUTINE HBANGLE
!----------------------------------------------------------------------------
!---------------------------------------------------------
!rescaling according to the reference and APPlying periodic boundary

SUBROUTINE rescalerefo(rx,ry,rz,CMx,CMy,CMz)
IMPLICIT NONE
INTEGER nsolute
PARAMETER(nsolute=33)
REAL*8::cmx,cmy,cmz,sumx,sumy,sumz,rdx(nsolute),rdy(nsolute),rdz(nsolute),rx(sol),ry(sol),rz(sol),b
REAL*8::dist_ref_O
INTEGER::i,j,idoi,sol, count , array(100)
CHARACTER(len=10)::ele(nsolute)

sol = 261                        !number of water molecules
b   = 14.861d0                     !box length of W-elp


DO i=1,sol
rx(i)=rx(i)-CMx
rx(i)=rx(i)-b*anint(rx(i)/b)
ry(i)=ry(i)-CMy
ry(i)=ry(i)-b*anint(ry(i)/b)
rz(i)=rz(i)-CMz
rz(i)=rz(i)-b*anint(rz(i)/b)
END DO


END SUBROUTINE rescalerefo
!----------------------------------------------------------------------------



!----------------------------------------------------------------------------
!SUBROUTINES

!CONTAINS




!---------------------------------------------------------
!rescaling according to the reference and APPlying periodic boundary

!SUBROUTINE rescalerefo(rx,ry,rz,CMx,CMy,CMz,sol,b)
!IMPLICIT NONE
!REAL*8       :: CMx,CMy,CMz
!REAL*8       :: rx(sol),ry(sol),rz(sol)
!REAL*8       :: b
!INTEGER      :: i,sol 



!DO i=1,sol
!rx(i)=rx(i)-CMx
!rx(i)=rx(i)-b*anint(rx(i)/b)
!ry(i)=ry(i)-CMy
!ry(i)=ry(i)-b*anint(ry(i)/b)
!rz(i)=rz(i)-CMz
!rz(i)=rz(i)-b*anint(rz(i)/b)
!END DO
!
!
!END SUBROUTINE rescalerefo
!----------------------------------------------------------------------------


   INTEGER FUNCTION  FindMinimum(x, Start, End)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                :: Start, End
      INTEGER                            :: Location
      INTEGER                            :: i
      REAL*8, DIMENSION(1:), INTENT(IN)  :: x
      REAL*8                             :: Minimum

      Minimum  = x(Start)                ! Assume the first is the min
      Location = Start                   ! Record its position
      DO i = Start+1, End                ! Start with next elements
         IF (x(i) < Minimum) THEN        ! If x(i) less than the min?
            Minimum  = x(i)              ! Yes, a new minimum found
            Location = i                 ! Record its position
         END IF
      END DO
      FindMinimum = Location             ! return the position
   END FUNCTION  FindMinimum

! --------------------------------------------------------------------
! SUBROUTINE  Swap():
!    This subroutine swaps the values of its two formal arguments.
! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      REAL*8, INTENT(INOUT) :: a, b
      REAL*8                :: Temp
      Temp = 0
      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

! --------------------------------------------------------------------
! SUBROUTINE  Sort():
!    This subroutine receives an array x() and sorts it into ascending
! order.
! --------------------------------------------------------------------

   SUBROUTINE  Sort(x, Size)
      IMPLICIT  NONE
      REAL*8, DIMENSION(1:), INTENT(INOUT)  :: x
      INTEGER, INTENT(IN)                   :: Size
      INTEGER                               :: i
      INTEGER                               :: Location
      DO i = 1, Size-1                          ! except for the last
         Location = FindMinimum(x, i, Size)     ! find min from this to last
         CALL  Swap(x(i), x(Location))          ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort




END PROGRAM ELPTICINE




