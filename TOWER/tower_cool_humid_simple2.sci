function [Z]=tower_cool_humid_simple (w_ar_ratio ,m_dot_w  , np ,lamda , p_barm  ,apro ,kmax ,tw_in, tar_in , hum_rel_in )
   le = 1
    
           t=tar_in;
           hum_rel=hum_rel_in;
           [ twb ]=wet_bulb(t ,hum_rel);
           t_wet_bulb=twb
           //tw_in = 80;
        disp 'inlet wet bulb temperature'
        disp (twb)
           // kmax=11;
        Dt_w=(tw_in-twb-apro)/kmax
                  // m_dot_w =6  //1.66      //ate 6.66
                 //  w_ar_ratio=1.6
                 //  tw=zeros(kmax+1)
                    m_dot_ar =m_dot_w/w_ar_ratio
           ratio_agua_ar = m_dot_w / m_dot_ar;
                //  p_barm = 1;
           r_j = 287;
           r_j_w = 461.51;
           cpar = 1;
           cpw = 4.18;
          // lamda = .417;                //%309;
          // np =.47;    
           
            t_k = tar_in+273.15;
           [ p ] = saturation(t_k );
           //%ttar_in = tt
           //%ftp_ar_in = ftp
           p_tar_in = p*hum_rel_in;
           disp 'vapor pressure at inlet air temperature'
           disp (p_tar_in)
           //%dftp_dtt_ar_in = dftp_dtt
           //%dftp_dt_ar_in = dftp_dt
           p=p_tar_in
         state=1
          [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
      // [ ro ] = steam_kenan( p ,t_k , state)  // ro kg/l
           ro_vap_tar_in =ro*1000                         //p_tar_in * 100000 / r_j_w / t_k;
           
   //ro_vap_tar_in_r= p_tar_in * 100000 / r_j_w / t_k;
  
  
           ro_ar_tar_in_seco = ((p_barm -p_tar_in ) * 100000 / (r_j * (tar_in + 273.15)));
      
           hum_ar_in = (ro_vap_tar_in ) / (ro_ar_tar_in_seco+ro_vap_tar_in);// kg de vapor por kg de ar seco
           m_dot_ar_seco=m_dot_ar*ro_ar_tar_in_seco/(ro_ar_tar_in_seco+ro_vap_tar_in)
   disp 'specific humidity at inlet air  temperature'
   disp (hum_ar_in)
           h_vap_tar_in = h1
           
          // dry_wet_f=ro_ar_tar_in/(ro_ar_tar_in+ro_vap_tar_in)
           hmist =(ro_vap_tar_in * h_vap_tar_in + ro_ar_tar_in_seco *cpar * tar_in)/(ro_vap_tar_in + ro_ar_tar_in_seco); // entalpia por kg de ar seco
           hg_in = hmist
    disp 'inlet enthalpy of air'
    disp (hg_in)
    
       m_dot_ar_new(1)=m_dot_ar
        m_dot_w_new(1)=m_dot_w
      
        Qtot_div(1)= m_dot_w_new(1)*cpw*Dt_w
        t_w(1)=tw_in-Dt_w*(kmax)
        t_ar(1)=tar_in
        hum_vap_t(1)=hum_ar_in
        T_k(1)=t_w(1)+273.15
        t_k=T_k(1)
         [ p ] = saturation(t_k );
          
      state=1
       [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
                ro_vap_t_w(1) =ro*1000            //p_tw(k) * 100000 / r_j_w / ttw(k);

                ro_ar_t_w(1) = (p_barm- p ) * 100000 / r_j / t_k;
                H_s(1)=(h1*ro_vap_t_w(1)+ro_ar_t_w(1)*cpar*t_w(1) )/(ro_vap_t_w(1)+ro_ar_t_w(1))

              f_s(1)=ro_vap_t_w(1)/(ro_vap_t_w(1)+ ro_ar_t_w(1))
              
              H_g(1)=hg_in
              f_g(1)=hum_vap_t(1)
              
              dfg(1)= m_dot_w/m_dot_ar*cpw*Dt_w/(H_s(1)-H_g(1))/(1-f_s(1))*(f_s(1)-f_g(1))
             
              D_ip= lamda * (m_dot_w / m_dot_ar) ^ (-np)
              Dz(1)=cpw*Dt_w/(H_s(1)-H_g(1))/D_ip
            
              H_g(2)=H_g(1)+ Qtot_div(1)/m_dot_ar_new(1)
              f_g(2)=f_g(1)+dfg(1)
             // t_ar(2)=(H_g(2)-f_g(2)*2501.6)/cpar
              //disp t_ar(2)
              //disp (t_ar(2))
              // H_g(2)mass vapor/)masvapo+massar)
             
             hum_vap_t(2)=f_g(2)         ///(1-f_g(2))   //hum_vap_t2=kg vapor/kg air moist
              
              
              //t_ar(2)=(H_g(2)-hum_vap_t(2)*2501.6)/(hum_vap_t(2)*1.827+(1-hum_vap_t(2))*cpar)
             
              m_dot_ar_new(2)=(1+hum_vap_t(1))*m_dot_ar
              m_dot_w_new(2)=m_dot_w+hum_vap_t(1)*m_dot_ar
              
              Qtot_div(2)=m_dot_w_new(2)*cpw*Dt_w
              
              
              for k =2 :(kmax)
                  t_w(k)=t_w(k-1)+Dt_w
        //t_ar(k)=tar_in
        //hum_ar(k)=hum_vap_t(k)
        t_k=t_w(k)+273.15
         T_k(k)=t_w(k)+273.15
        t_k=T_k(k)
         [ p ] = saturation(t_k );
          //      %ttw(k) = tt
           //     %ftp_w_in = ftp
                
             //   %dftp_dtt_w_in = dftp_dtt
              //  %dftp_dt_w_in = dftp_dt

  
   // [ ro ] = steam_kenan( p ,t_k , state)  // ro kg/l
   
  //p=p_tar_in
      state=1
       [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
                ro_vap_t_w(k) =ro*1000            //p_tw(k) * 100000 / r_j_w / ttw(k);

                ro_ar_t_w(k) = (p_barm- p ) * 100000 / r_j / t_k;
               
                H_s(k)=(h1*ro_vap_t_w(k)+ro_ar_t_w(k)*cpar*t_w(k))/(ro_vap_t_w(k)+ro_ar_t_w(k)) 

              f_s(k)=ro_vap_t_w(k)/(ro_vap_t_w(k)+ ro_ar_t_w(k))
              
              //f_g(k)=hum_ar(k)
             H_g(k)=H_g(k-1)+ Qtot_div(k-1)/m_dot_ar_new(k-1)
              dfg(k)= m_dot_w_new(k)/m_dot_ar_new(k)*cpw*Dt_w/(H_s(k)-H_g(k))/(1-f_s(k))*(f_s(k)-f_g(k))
              d_Hs_Hg(k)=H_s(k)-H_g(k)
              if d_Hs_Hg(k)< 0 then
                disp  'PROCESSING STOPPED DUE TO INOPUT DATA NOT SUITABLE- H_S IS LESS THAN h_G  '
                  break
                  end
              // dfg(k)= m_dot_w_new(k-1)/m_dot_ar_new(k-1)*Qtot_div(k-1)/(H_s(k)-H_g(k))/(1-f_s(k))*(f_s(k)-f_g(k))
               // f_g(k+1)=f_g(k)+dfg(k)
             // H_g(k)=hg_in
              D_ip= lamda * (m_dot_w_new(k-1) / m_dot_ar_new(k-1)) ^ (-np)
              Dz(k)=cpw*Dt_w/(H_s(k)-H_g(k))/D_ip
              //disp Dz(k)
              //disp (Dz(k))
             // H_g(k+1)=H_g(k)+ Qtot_div(k-1)/m_dot_ar_new(k-1)
              // H_g(2)=H_g(1)+ Qtot_div(1)/m_dot_ar_new(1)
              f_g(k+1)=f_g(k)+dfg(k)
             
             
             // t_ar(2)=(H_g(2)-f_g(2)*2501.6)/cpar
             // disp t_ar(2)
              //disp (t_ar(2))
              // H_g(2)mass vapor/)masvapo+massar)
             
             hum_vap_t(k+1)=f_g(k+1)/(1-f_g(k+1))   //hum_vap_t2=kg vapor/kg air moist
              
              
              //t_ar(k+1)=(H_g(k+1)-hum_vap_t(k+1)*2501.6)/(hum_vap_t(k+1)*1.827+(1-hum_vap_t(k+1))*cpar)
              
              //disp t_ar(k+1)
              //disp (t_ar(k+1))
              m_dot_ar_new(k+1)=(1+hum_vap_t(k))*m_dot_ar
              m_dot_w_new(k+1)=m_dot_w+(hum_vap_t(k)*m_dot_ar)
                Qtot_div(k)=m_dot_w_new(k)*cpw*Dt_w
              end
              
             // for k=1: kmax-1
             //     D_hs_hg(k)=H_s(k)-H_g(k)
              //end
             // disp D_hs_hf
              //disp(D_hs_hg)
              //m_dot_ar_new(k)=m_dot_ar_new(k-1)+dfg(k)*m_dot_ar                                     //(1+hum_vap_t(k+1))*m_dot_ar
              //m_dot_w_new(k)=m_dot_w_new(k-1)+dfg(k)*m_dot_ar                             //m_dot_w*(1+hum_vap_t(k+1)*m_dot_ar)
              //Qtot_div(k+1)=m_dot_w_new(k)*cpw*Dt_w
             // end
             // for k=kmax-1 : kmax+1
             //========================
             t_w_inlet  = t_w(kmax)+Dt_w
            
             //==========================
        //t_ar(k)=tar_in
        //hum_ar(k)=hum_vap_t(k)
       // t_k=t_w_exit+273.15
       //===================================
     //   T_k_exit=t_w_exit+273.15
      // t_k=T_k_exit
        // [ p ] = saturation(t_k );
        
   
  //p=p_tar_in
     // state=1
      // [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
      

        //       ro_ar_t_w(exit) = (p_barm- p ) * 100000 / r_j / t_k;
               
              
           H_g_exit=H_g(kmax)+ Qtot_div(kmax)/m_dot_ar_new(kmax)
            
             
   m_dot_ar_new_exit=(1+hum_vap_t(kmax))*m_dot_ar
   m_dot_w_new_inlet =m_dot_w+(hum_vap_t(kmax)*m_dot_ar)
              
              
             //===================================================================== 
             // pressure loss
             //Np=np*2*g/v^2*dw/dair
             
             Z=sum(Dz)
             //====================================================
             ro_ar_in=( ro_vap_tar_in + ro_ar_tar_in_seco)
           v=m_dot_ar/(ro_ar_in)
           cons_dp=2*9.81/v^2*1000/ro_ar_in
           Np=3.401*m_dot_w^.59*m_dot_ar^-.542
           np=1000*Np/cons_dp*Z
          
               disp 'exit air specific humidity'
              disp (f_g(kmax+1))
               
             disp 'exit water temperature'
             disp (t_w(1))
               disp 'inlet water temperature'
              disp (t_w_inlet)
               disp  'exit air enthalpy'
              disp ( H_g_exit)
              disp   'exit air mass flow rate'
              disp (  m_dot_ar_new_exit)
              disp   'inlet water mass flow rate'
              disp (  m_dot_w_new_inlet)
             //===================================================================== 
             
             disp 'packing height meter' 
             disp (Z)
        
              Q_Ztot=sum(Qtot_div)
              disp 'heat transfer load kW/square meter'
              disp(Q_Ztot)
              
              t_ar_out=(H_g_exit-2501.5*f_g(kmax+1))/(cpar*(1-f_g(kmax))+1.83*f_g(kmax+1))
              disp 'exit air temperature'
              disp (t_ar_out)
              
               disp 'pressure loss mm water' 
           disp (np)
             
    endfunction
    
////////////////////////////////////////////////////////////////////////////////////
//tic()
//getd ("C:\Users\wsb\Documents\share\");
//tar_in=30
//hum_rel_in=.2
//[Z]=tower_cool_humid_simple( tar_in , hum_rel_in )

    
    
