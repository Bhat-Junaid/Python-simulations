program spectral_analysis
    implicit none
    integer :: i 
    integer, parameter :: n = 861
    real ,dimension(n) :: fexp ,ftheo_wi1 ,ftheo_wi2 , omega, omega2, wi1, wi2, x, x2, x3, x4, y, yx, yx2 
    real ::  sum1, sum2, sum3, sum4, sum5, sum6, sum7, sum8, sum9, sum10, sum11, sum12, sum13, sum14, sum15, sum16 ,sum17, S1, S2
    real ::  avg1_x, avg1_x2, avg1_x3, avg1_x4, avg1_y, avg1_yx, avg1_yx2, avg2_1, avg2_x, avg2_x2, avg2_x3, avg2_x4, avg2_y, wm1
    real :: avg2_yx, avg2_yx2, omega_mean, omega2_mean, omega_mean2, sigma , a20, a21, a22, a10, a11, a12, pi, gamma1, gamma2, wm2   
    
    pi = 4.0*ATAN(1.0)
    
    ! Reading the spectrum values stored in spectre_1.txt
    open(unit=5,file='spectre_5.txt') ! open the file 
    sum1 = 0. 
    sum2 = 0.
    do i=1,n
        read(5,*) fexp(i) ! read intensity values and store them into a vector
        omega(i) = 2280. + 0.01*(i-1) ! Omega vector with a constant time step and starting value of 2280
        omega2(i) = omega(i)**2 
        sum1 = sum1 + omega(i) ! Sum the omega terms
        sum2 = sum2 + omega2(i) ! Sum of the omega square terms
        write(12,*) omega(i)
    end do
    close(unit=5) ! close the file 
     
    omega_mean = sum1/float(n) ! The mean value of omega 
    omega2_mean = sum2/float(n) ! The mean value of the square of omega 
    omega_mean2 = omega_mean**2 ! the square of the mean value of omega 
    sigma = sqrt(omega2_mean - omega_mean2)! in order to calculate the value of 
    
    print*,"omega_mean: ", omega_mean
    print*, "omega_mean2: ", omega_mean2
    print*, "Sigma: ", sigma
    
! To calculate the various values x, x^2, x^3, x^4, y, xy, yx^2
    do i=1, n
        x(i)= (omega(i)-omega_mean)/sigma 
        y(i) = 1/ (fexp(i))
        yx(i) = x(i)*y(i)
        yx2(i) = yx(i)*x(i)
        x2(i) = x(i)*x(i)
        x3(i) = x2(i)*x(i)
        x4(i) = x2(i)*x2(i)
        wi1(i) = 1 ! Without weight function 
        wi2(i) = fexp(i)*fexp(i) ! calculate the weights w = fexp^2
        
    end do 
    
!Initializing various sums to find the means of the above quantities    
    sum3 = 0. ; sum4 = 0. ;sum5 = 0. ;sum6 = 0. ;sum7 = 0. ;sum8 = 0. ;sum9 = 0. ;sum10 = 0. ;sum11 = 0. ;sum12 = 0. ;sum13 = 0. 
    sum14 = 0. ;sum15 = 0. ;sum16 = 0.
    
    do i=1,n
        ! Summation of various terms corresponding to w=1
        sum3 = sum3 + wi1(i)*x(i) ! sum of x values
        sum4 = sum4 + wi1(i)*x2(i) ! sum of x^2 values
        sum5 = sum5 + wi1(i)*x3(i) ! sum of x^3 values
        sum6 = sum6 + wi1(i)*x4(i) ! sum of x^4 values
        sum7 = sum7 + wi1(i)*y(i) ! sum of y values
        sum8 = sum8 + wi1(i)*yx(i) ! sum of y*x values
        sum9 = sum9 + wi1(i)*yx2(i) ! sum of y*x^2 values
        ! Summation of various terms corresponding to W =(fext)^2
        sum10 = sum10 + wi2(i)*x(i) 
        sum11 = sum11 + wi2(i)*x2(i)
        sum12 = sum12 + wi2(i)*x3(i)
        sum13 = sum13 + wi2(i)*x4(i)
        sum14 = sum14 + wi2(i)*y(i)
        sum15 = sum15 + wi2(i)*yx(i)
        sum16 = sum16 + wi2(i)*yx2(i)

    end do
    
    ! mean values corresponding to W = 1
    avg1_x = sum3/float(n) ! mean of x
    avg1_x2 = sum4/float(n) ! mean of x^2
    avg1_x3 = sum5/float(n) ! mean of x^3
    avg1_x4 = sum6/float(n) ! mean  of x^4
    avg1_y = sum7/float(n) ! mean of y
    avg1_yx = sum8/float(n) ! mean of x*y
    avg1_yx2 = sum9/float(n) ! mean of x^2*y
    
    ! mean values corresponding to W =(fext)^2
    avg2_x = sum10/float(n)
    avg2_x2 = sum11/float(n)
    avg2_x3 = sum12/float(n)
    avg2_x4 = sum13/float(n)
    avg2_y = sum14/float(n)
    avg2_yx = sum15/float(n)
    avg2_yx2 = sum16/float(n)

    
    ! coefficients values a corresponding to W = 1, In this case <1> = 1
    a10 = ( avg1_y*avg1_x2*avg1_x4 + avg1_yx*avg1_x3*avg1_x2 + avg1_x*avg1_x3*avg1_yx2 - avg1_yx2*avg1_x2*avg1_x2 &
    &    - avg1_y*avg1_x3*avg1_x3 - avg1_x*avg1_yx*avg1_x4) / (avg1_x2*avg1_x4 + 2*avg1_x*avg1_x2*avg1_x3 - avg1_x2**3 &
    &    - avg1_x3*avg1_x3 - avg1_x*avg1_x*avg1_x4)
    print*, "a10 = ", a10 

    a11 = ( avg1_yx*avg1_x4 + avg1_x2*avg1_x3*avg1_y + avg1_x*avg1_x2*avg1_yx2 - avg1_x2*avg1_x2*avg1_yx & 
    &    - avg1_x*avg1_x4*avg1_y - avg1_x3*avg1_yx2 ) / ( avg1_x2*avg1_x4 + 2*avg1_x*avg1_x2*avg1_x3 &
    &    - avg1_x2*avg1_x2*avg1_x2 - avg1_x3*avg1_x3 - avg1_x*avg1_x*avg1_x4 )
    print*, "a11 = ", a11

    a12 = ( avg1_x2*avg1_yx2 + avg1_x*avg1_x3*avg1_y + avg1_x*avg1_x3*avg1_yx - avg1_x2*avg1_x2*avg1_y   &
    &    - avg1_x*avg1_x*avg1_yx2 - avg1_x3*avg1_yx ) / ( avg1_x2*avg1_x4 + 2*avg1_x*avg1_x2*avg1_x3 &
    &    - avg1_x2*avg1_x2*avg1_x2 - avg1_x3*avg1_x3 - avg1_x*avg1_x*avg1_x4 )
    print*, "a12 = ", a12
    
    ! coefficients values a corresponding to W = (fexp)^2
    ! calculation of <1> in the case of W = fexp^2 :
    sum17 = 0
    do i=1,n 
        sum17 = sum17 + wi2(i) 
    end do 
    avg2_1 = sum17/float(n) 

    a20 = ( avg2_y*avg2_x2*avg2_x4 + avg2_yx*avg2_x3*avg2_x2 + avg2_x*avg2_x3*avg2_yx2 - avg2_yx2*avg2_x2*avg2_x2 &
    &    - avg2_y*avg2_x3*avg2_x3 - avg2_x*avg2_yx*avg2_x4) / (avg2_1*avg2_x2*avg2_x4 + 2*avg2_x*avg2_x2*avg2_x3 - avg2_x2**3 &
    &    - avg2_1*avg2_x3*avg2_x3 - avg2_x*avg2_x*avg2_x4)
    print*, "a20 = ", a20
    a21 = ( avg2_1*avg2_yx*avg2_x4 + avg2_x2*avg2_x3*avg2_y + avg2_x*avg2_x2*avg2_yx2 - avg2_x2*avg2_x2*avg2_yx & 
    &    - avg2_x*avg2_x4*avg2_y - avg2_1*avg2_x3*avg2_yx2 ) / ( avg2_1*avg2_x2*avg2_x4 + 2*avg2_x*avg2_x2*avg2_x3 &
    &    - avg2_x2*avg2_x2*avg2_x2 - avg2_1*avg2_x3*avg2_x3 - avg2_x*avg2_x*avg2_x4 )
    print*, "a21 = ", a21
    a22 = ( avg2_1*avg2_x2*avg2_yx2 + avg2_x*avg2_x3*avg2_y + avg2_x*avg2_x3*avg2_yx - avg2_x2*avg2_x2*avg2_y   &
    &    - avg2_x*avg2_x*avg2_yx2 - avg2_1*avg2_x3*avg2_yx ) / ( avg2_1*avg2_x2*avg2_x4 + 2*avg2_x*avg2_x2*avg2_x3 &
    &    - avg2_x2*avg2_x2*avg2_x2 - avg2_1*avg2_x3*avg2_x3 - avg2_x*avg2_x*avg2_x4 )
    print*, "a22 = ", a22


    S1 = pi*sigma/sqrt(a10*a12 - a11*a11*0.25) ! calculation of the parameter S for W = 1
    S2 = pi*sigma/sqrt(a20*a22 - a21*a21*0.25) ! calculation of the parameter S for W = (fexp)^2
    gamma1 = Sigma*sqrt(a10/a12 - a11*a11/4/a12/a12) ! calculation of the gamma for W = 1
    gamma2 = Sigma*sqrt(a20/a22 - a21*a21/4/a22/a22) ! calculation of the gamma for W = (fexp)^2
    wm1 = omega_mean - Sigma*a11/2/a12 ! calculation of the maximum wavenumber for W = 1
    wm2 = omega_mean - Sigma*a21/2/a22 ! calculation of the maximum wavenumber for W = (fexp)^2
    print*, "S corresponding to (W = 1) = ", S1 
    print*, "gamma corresponding to (W = 1) = ", gamma1 
    print*, "wm corresponding to (W = 1) = ", wm1 
    print*, "S corresponding to (W = fexp^2) = ", S2  
    print*, "gamma corresponding to (W = fexp^2) = ", gamma2
    print*, "wm corresponding to (W = fexp^2) = ", wm2
    
    ! Calculation of theoritical intensities corresponding to w = 1 and storing them in a file
    ! We will use python to plot these intensities 
    open(unit=20,file='theo_wi1_spectre_5.txt') ! open a new file 
    do i= 1, n
        ftheo_wi1(i) = S1*gamma1/((omega(i)-wm1)*(omega(i)-wm1)+gamma1*gamma1)/pi
        write(20,*) ftheo_wi1(i) ! put the value in the file
    end do
    close(20) 
    
    ! Calulation of intensities with w = Fexp^2   
    open(unit=20,file='theo_wi2_spectre_5.txt') 
    do i= 1, n
        ftheo_wi2(i) = S2*gamma2/((omega(i)-wm2)*(omega(i)-wm2)+gamma2*gamma2)/pi
        write(20,*) ftheo_wi2(i)  
    end do
    close(20) 

end program spectral_analysis

! In the program, i am getting different a0, a1 and a2 values because of some error in the standard calculation but 
! as you run the program you will see that i am correct values of S, gamma and wm, i did not understand why this was happening.
!one explanation is that due to different standard deviaition, the values of a0,a1, a2 change in such a way that it has no 
!effect the results
! I calculated the standard deviation in python and got correct values of a0, a1 and a2
!Python code to calculate the standard deviation so that i got correct values of a0,a1,a2
!(which is unneccesary because i am getting correct values of parameters)

!import numpy as np
!omega = np.loadtxt(r"C:\Academics\Numerical Methods\omega_values.txt")
!mean_1= np.mean(omega)
!w = np.zeros(len(omega))
!for i in range(len(omega)):
!    w[i] = (omega[i])**2
!print()
!mean = np.mean(w)
!std = np.sqrt(mean - (mean_1)**2)
!std1 = np.std(omega)
!print(mean_1, mean, std, std1)
