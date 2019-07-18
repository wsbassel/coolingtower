function [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
    // function for estimation of enthalpy and density of steam

          R=4.6151;
          R_j=R/10;
 
if state <2 
ro=p/(R*t_k);
else
    ro=.99;
end
	am=zeros(10,7);

        am(1,:) = [29.492937 -5.198586  6.8335354   -0.1564104   -6.3972405   -3.9661401  -0.69048554];

		am(2,:) =[ -132.13917 7.7779182  -26.149751  -0.72546108  26.409282  15.453061  2.7407416];

		am(3,:)= [274.64632	 -33.301902  65.326396  -9.2734289  -47.740374  -29.14247  -5.102807];

		am(4,:) =[ -360.93828  -16.254622  -26.181978  4.312584  56.32313  29.568796  3.9636085];

		am(5,1) =342.18431;
		am(5,2)=-177.3107;


		am(6,1)=-244.50042;
		am(6,2)= 127.48742;



		am(7,1)=155.18535;
		am(7,2)=137.46153;
		



		am(8,1) = 5.9728487;
		am(8,2) = 155.97836;
		
		


		am(9,:)=[ -410.30848  337.3118  -137.46618  6.7874983  136.87317  79.84797  13.041253];



		am(10,:) =[ -416.0586    -209.88866  -733.96848  10.401717  645.8188  399.1757  71.531353];
        
 
roa=zeros(7,1);
taua=zeros(7,1);
		roa(1)=.634;
		taua(1)=1.544912;
	
	for j=2:7

	roa(j)= 1.0;
	taua(j)= 2.5;
    end

		tauc = 1.544912;
		em = 4.8;

	tau=1000.0/t_k;


	if tau == 2.5   then   //& tau<2.500001;
	tau = 2.5001;
    end
		



//	                              %pow(e_nap,m_e);
	E=4.8;
     
          for iter=1:30
	for j=1:7
        con_Q1(j)=(tau-taua(j))^(j-2);
            for i=1:8
                con_Q2(i)=am(i,j)*(ro-roa(j))^(i-1);
            end
            con_Q2(9)=exp(-E*ro)*am(9,j);
            con_Q2(10)=exp(-E*ro)*am(10,j)*ro;
            con_Q3(j)=con_Q1(j)*sum(con_Q2);
    end
    Q=(tau-tauc)*sum(con_Q3);
    
                        
                       for j=1:7
        con_dQ1(j)=(tau-taua(j))^(j-2);
            for i=1:8
                con_dQ2(i)=am(i,j)*(i-1)*(ro-roa(j))^(i-2);
            end
            con_dQ2(9)=exp(-E*ro)*(-E)*am(9,j);
            con_dQ2(10)=exp(-E*ro)*am(10,j)+exp(-E*ro)*(-E)*am(10,j)*ro;
            con_dQ3(j)=con_dQ1(j)*sum(con_dQ2);
    end
    dQ_dro=(tau-tauc)*sum(con_dQ3);
   
                      
                       for j=1:7
        con_d2Q1(j)=(tau-taua(j))^(j-2);
            for i=1:8
                con_d2Q2(i)=am(i,j)*(i-1)*(i-2)*(ro-roa(j))^(i-3);
            end
            con_d2Q2(9)=exp(-E*ro)*(-E)*(-E)*am(9,j);
            con_d2Q2(10)=exp(-E*ro)*(-E)*am(10,j)+exp(-E*ro)*(-E)*am(10,j)+exp(-E*ro)*(-E)*(-E)*am(10,j)*ro;
            con_d2Q3(j)=con_d2Q1(j)*sum(con_d2Q2);
    end
    d2Q_dro=(tau-tauc)*sum(con_d2Q3);
        
         
          
		f=-p/(R*t_k)+ro+ro^2*Q+ro^3*dQ_dro;
        df_dro=1+ro^2*dQ_dro+2*ro*Q+ro^3*d2Q_dro+3*ro^2*dQ_dro;
      dro=f/df_dro;
      ro=ro-dro/2;
          end
          ro1=ro;
          //%================dQ_dtau===============
			for j=1:7
        con_dQ_dtau1(j)=(j-2)*(tau-taua(j))^(j-3);
            for i=1:8
                con_dQ_dtau2(i)=am(i,j)*(ro-roa(j))^(i-1);
            end
            con_dQ_dtau2(9)=exp(-E*ro)*am(9,j);
            con_dQ_dtau2(10)=exp(-E*ro)*am(10,j)*ro;
            con_dQ_dtau3(j)=con_dQ_dtau1(j)*sum(con_dQ_dtau2);
    end
    dQ_dtau=sum(con_Q3)+(tau-tauc)*sum(con_dQ_dtau3);
    
    //%====================dcon_dQ_dtau3===================
    for j=1:7
        dcon_dQ_dtau1(j)=(j-3)*(j-2)*(tau-taua(j))^(j-4);
            for i=1:8
                dcon_dQ_dtau2(i)=am(i,j)*(ro-roa(j))^(i-1);
            end
            dcon_dQ_dtau2(9)=exp(-E*ro)*am(9,j);
            dcon_dQ_dtau2(10)=exp(-E*ro)*am(10,j)*ro;
            dcon_dQ_dtau3(j)=dcon_dQ_dtau1(j)*sum(dcon_dQ_dtau2);
    end
    dcon_dQ_dtau3s=sum(dcon_dQ_dtau3);
    //%======================================
    //%===========================dsum(con_q3=================
    for j=1:7
        dcon_Q1(j)=(j-2)*(tau-taua(j))^(j-3);
            for i=1:8
                dcon_Q2(i)=am(i,j)*(ro-roa(j))^(i-1);
            end
            dcon_Q2(9)=exp(-E*ro)*am(9,j);
            dcon_Q2(10)=exp(-E*ro)*am(10,j)*ro;
            dcon_Q3(j)=dcon_Q1(j)*sum(dcon_Q2);
    end
    dsum_con_Q3=sum(dcon_Q3);
    //%==================================
    d2Q_dtau2=dsum_con_Q3+(tau-tauc)*dcon_dQ_dtau3s+sum(con_dQ_dtau3);
    
    
    
     c_epsi=[1857.065  3229.12 -419.465 36.6649 -20.5516 4.85233 46 -1011.249];
        epsi=zeros(8,1);
        for i=1:6
            epsi(i)=c_epsi(i)/tau^(i-1);
        end
          epsi(7)=c_epsi(7)*log(t_k);
          epsi(8)=c_epsi(8)*log(t_k)/tau;
          
          epsi_zero=sum(epsi);
      for i=1:6
            depsi_tau(i)=-(i-1)* c_epsi(i)*tau^(-(i-1)-1);
        end
          depsi_tau(7)=-c_epsi(7)/tau;                      //  %(1/(1000/tau)*(-1000)/tau^2);
          depsi_tau(8)=c_epsi(8)*(-1/tau^2-log(t_k)/tau^2);
          
          depsi_tau_zero=sum(depsi_tau);
          //%==========================depsi2_dtau2=================
           for i=1:6
            d2epsi_tau(i)=-(i-1)*(-i)* c_epsi(i)*tau^(-i-1);
        end
          d2epsi_tau(7)=c_epsi(7)/tau^2;                       // %(1/(1000/tau)*(-1000)/tau^2);
          d2epsi_tau(8)=c_epsi(8)*(2/tau^3+1/tau^3+2*log(t_k)/tau^3);
          
          d2epsi_tau_zero=sum(d2epsi_tau);
          //%=============================
         
          
          
           u1=R_j*t_k*ro*tau*dQ_dtau+tau*depsi_tau_zero+epsi_zero;
          
              
          h1=u1+R_j*t_k*(1+ro*Q+ro^2*dQ_dro);
          s1=-R_j*(log(ro)+ro*Q-ro*tau*dQ_dtau)+depsi_tau_zero*1000/(t_k^2);
// %============================repetit com .1 c mais======

endfunction
//p=1
//t_k=373.15
//state=1

  //  [ ro , u1 , h1, s1] =wageeh_keenan (p, t_k, state)
   // disp (ro)
   // disp (u1)
    //disp (h1)
    //disp (s1)
