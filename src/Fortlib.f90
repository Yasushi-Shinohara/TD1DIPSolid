Subroutine uGbk_forward_RK4(uGbk, hGGk, NG, Nocc, Nk, dt)
  Implicit none
  Complex(8), parameter :: zI=(0.0d0, 1.0d0)
  Integer(8), intent(in) :: NG, Nocc, Nk
  Real(8), intent(in) :: dt
  Complex(8), intent(inout) :: uGbk(1:Nk,1:Nocc,1:NG)
  Complex(8), intent(in) :: hGGk(1:Nk,1:NG,1:NG)
  Integer :: ik, ig, ib
  Complex(8) :: k1(1:Nocc,1:NG), k2(1:Nocc,1:NG), k3(1:Nocc,1:NG), k4(1:Nocc,1:NG)

!  Write(*,*) NG, Nocc, Nk, dt
!  Stop
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
    Complex(8), intent(in) :: u(1:Nocc, 1:NG), h(1:NG, 1:NG)
    Complex(8) :: hu(1:Nocc, 1:NG)
    Integer(8) :: i

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
