Subroutine uGbk_forward_RK4(uGbk, hGGk, NG, Nocc, Nk, dt)
  Implicit none
  Complex(kind(0d0)), parameter :: zI=(0.0d0, 1.0d0)
  Integer, intent(in) :: NG, Nocc, Nk
  Double precision, intent(in) :: dt
  Complex(kind(0d0)), intent(inout) :: uGbk(1:Nk,1:Nocc,1:NG)
  Complex(kind(0d0)), intent(in) :: hGGk(1:Nk,1:NG,1:NG)
  Integer :: ik, ig, ib
  Complex(kind(0d0)) :: k1(1:Nocc,1:NG), k2(1:Nocc,1:NG), k3(1:Nocc,1:NG), k4(1:Nocc,1:NG)

  Do ik = 1, Nk
    k1 = u_h2hu(uGbk(ik,:,:)              , hGGk(ik,:,:))/zI
    k2 = u_h2hu(uGbk(ik,:,:) + 0.5d0*dt*k1, hGGk(ik,:,:))/zI
    k3 = u_h2hu(uGbk(ik,:,:) + 0.5d0*dt*k2, hGGk(ik,:,:))/zI
    k4 = u_h2hu(uGbk(ik,:,:) + dt*k3      , hGGk(ik,:,:))/zI
    uGbk(ik,:,:) = uGbk(ik,:,:) + (k1 + 2.0d0*k2 + 2.0d0*k3 + k4)*dt/6.0d0
  End Do
  Return
Contains
  Function u_h2hu(u, h) result(hu)
    Implicit none
    Complex(kind(0d0)), intent(in) :: u(1:Nocc, 1:NG), h(1:NG, 1:NG)
    Complex(kind(0d0)) :: hu(1:Nocc, 1:NG)
    Integer :: i

!    hu(1:Nocc, 1:NG) = 0.d0
!    Do ig = 1,NG
!      Do ib = 1,Nocc
!        Do i = 1,NG
!          hu(ib,ig) = hu(ib,ig) + h(i,ig)*u(ib,i)
!        End Do
!      End Do 
!    End Do
    hu = matmul(u,h)
  End Function u_h2hu
End Subroutine uGbk_forward_RK4
!==========================================================================================
Subroutine uGbk_forward_exp(uGbk, hGGk, NG, Nocc, Nk, dt)
  Implicit none
  Complex(kind(0d0)), parameter :: zI=(0.0d0, 1.0d0)
  Integer, intent(in) :: NG, Nocc, Nk
  Double precision, intent(in) :: dt
  Complex(kind(0d0)), intent(inout) :: uGbk(1:Nk,1:Nocc,1:NG)
  Complex(kind(0d0)), intent(in) :: hGGk(1:Nk,1:NG,1:NG)
  Integer :: ik, ig, ib
  Complex(kind(0d0)) :: U(1:NG,1:NG)

  Do ik = 1, Nk
    U = hdt2U(hGGk(ik,:,:)*dt)
    uGbk(ik,:,:) = U_u2Uu(U, uGbk(ik,:,:))
  End Do
  Return
Contains
  !==
  Function hdt2U(hdt) result(U)
    Implicit none
    Complex(kind(0d0)), intent(in) :: hdt(1:NG, 1:NG)
    Complex(kind(0d0)) :: U(1:NG, 1:NG)
    Complex(kind(0d0)) :: coef(1:NG, 1:NG)
    Double precision :: eigs(1:NG)
    Integer :: i,j,k

    U = (0.d0, 0.d0)
    Call eigh(NG,transpose(hdt),eigs,coef) !Transpose is for the conversion from Row-major to Column-major
    Do i = 1,NG
      Do J = 1,NG
        Do k = 1,NG
          U(j,k) = U(j,k) + exp(-zI*eigs(i))*coef(j,i)*conjg(coef(k,i))
        End Do
      End Do
    End Do
    U = transpose(U)                      !Transpose is for the conversion from Column-major to Row-major 
  End Function hdt2U
  !==
  Function U_u2Uu(Unitary, uorb) result(Unitaryuorb)
    Implicit none
    Complex(kind(0d0)), intent(in) :: uorb(1:Nocc, 1:NG), Unitary(1:NG, 1:NG)
    Complex(kind(0d0)) :: Unitaryuorb(1:Nocc, 1:NG)
    Integer :: i 

!    Uu(1:Nocc, 1:NG) = 0.d0
!    Do ig = 1,NG
!      Do ib = 1,Nocc
!        Do i = 1,NG
!          Uu(ib,jg) = Uu(ib,ig) + U(i,ig)*u(ib,i)
!        End Do
!      End Do 
!    End Do
    Unitaryuorb = matmul(uorb,Unitary)
  End Function U_u2Uu
End Subroutine uGbk_forward_exp
!==========================================================================================
Subroutine eigh(N,H,E,V)
  Implicit none
  Integer, intent(in) :: N
  Complex(kind(0d0)), intent(in) :: H(1:N, 1:N)
  Double precision, intent(out) :: E(1:N)
  Complex(kind(0d0)), intent(out) :: V(1:N, 1:N)
!For ZHEEV
  Integer :: LWORK_EV 
  Complex(kind(0d0)) :: WORK_EV(2*(2*N-1)) !The argument is just LWORK_EV
  Double precision :: RWORK_EV(3*N - 2)
  Integer :: INFO_EV
  LWORK_EV = 2*(2*N-1)

  V = H
  Call ZHEEV('V','U',N,V,N,E,WORK_EV,LWORK_EV,RWORK_EV,INFO_EV)
    
  Return
End Subroutine eigh
