function [Z]=tower_cool_humid_simple (w_ar_ratio ,m_dot_w  , np ,lamda , p_barm  ,apro ,kmax ,tw_in, tar_in , hum_rel_in )
   le = 1
    
           t=tar_in;  // inlet temperature of air degrees centigrade
           hum_rel=hum_rel_in;   // relative humidity 
           [ twb ]=wet_bulb(t ,hum_rel);
           t_wet_bulb=twb     // wet bulb temperature 
           
        disp 'inlet wet bulb temperature'
        disp (twb)
          
        Dt_w=(tw_in-twb-apro)/(kmax-1)  // tw_in  : inlet  water temperature
                                    // twb : wet bulb temperature 
                                    //apro : approach degrees centigrade
                                    // kmax ; number of packing divisions +1      
                  
                 
                 
                   
           m_dot_ar =m_dot_w/w_ar_ratio // m_dot_ar: air mass rate of flow kg/m^2/sec
                                        // m_dot_w : water mass rate of flow kg/m^2/sec
                                        //w_ar_ratio : ratio of water to air mass rates of flow   
                    
           //ratio_agua_ar = m_dot_w / m_dot_ar;
                ratio_agua_ar =w_ar_ratio 
           r_j = 287;      //  air gas constant j/kg/K
           r_j_w = 461.51;  // water vapor gas constant j/kg/K
           cpar = 1;         // cpar specific heat of air kj/kg/K
           cpw = 4.18;       // specific heat of water kj/kg/K 
          
           
           
            t_k = tar_in+273.15;  //t_k : air inlet temperature degrees K 
           [ p ] = saturation(t_k );  // p :saturation pressure bar
           
           
           p_tar_in = p*hum_rel_in; // p_tar_in : partial pressure of vapor at inlet
           disp 'vapor pressure at inlet air temperature'
           disp (p_tar_in)
           
           
           p=p_tar_in
         state=1    // state=1 means estamate thermpdynamic properties for vaapor phase
          
          [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)// function to calculate enthalpy and density
      
           ro_vap_tar_in =ro*1000 // ro_vap_tar_in : density of vapor phase kg/m^3
          
  
           ro_ar_tar_in_seco = ((p_barm -p_tar_in ) * 100000 / (r_j * (tar_in + 273.15)));
                                       // ro_ar_tar_in_seco : density of dry air at inlet kg/m^3
      
           hum_ar_in = (ro_vap_tar_in ) / (ro_ar_tar_in_seco+ro_vap_tar_in);
                                // hum_ar_in : specific humidity of air at inlet
           m_dot_ar_seco=m_dot_ar*ro_ar_tar_in_seco/(ro_ar_tar_in_seco+ro_vap_tar_in)
   disp 'specific humidity at inlet air  temperature'
   disp (hum_ar_in)
           h_vap_tar_in = h1  // enthalpy of air at inlet tempoerature
           
          
           hmist =(ro_vap_tar_in * h_vap_tar_in + ro_ar_tar_in_seco *cpar * tar_in)/(ro_vap_tar_in + ro_ar_tar_in_seco); 
           hg_in = hmist// h_mist enthalpy of humid air at inlet kj/kg 
    disp 'inlet enthalpy of air'
    disp (hg_in)
    //--------------------------------------------------------------------------
    //calculation of the height and the change of humidity of airof 1st division 
         
    
    
       m_dot_ar_new(1)=m_dot_ar   // mass rate of flow of air INLET of the 1st division kg/m^2/sec
        m_dot_w_new(1)=m_dot_w    // mass rate of flow of water at EXIT of the 1st division
        Qtot_div(1)= m_dot_w_new(1)*cpw*Dt_w // energy transfer from water to air of the 1st division
        t_w(1)=tw_in-Dt_w*(kmax-1) // exit water temperature requerid  degrees centigrade 
        t_ar(1)=tar_in           // inlet air temperature degrees centigrade of the 1st division
        hum_vap_t(1)=hum_ar_in   // specific humidity of inlet air of 1st division
        
        
        T_k(1)=t_w(1)+273.15
        t_k=T_k(1) // absolute temperature of water at exit of 1st division
         [ p ] = saturation(t_k );  // saturation pressure 
          
      state=1
       [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state) // function to find enthalpy and density 
                ro_vap_t_w(1) =ro*1000   // density of water vapor at temperature of water at exit                         

                ro_ar_t_w(1) = (p_barm- p ) * 100000 / r_j / t_k;
                H_s(1)=(h1*ro_vap_t_w(1)+ro_ar_t_w(1)*cpar*t_w(1) )/(ro_vap_t_w(1)+ro_ar_t_w(1))
                //H_s(1): enthalpy of saturated vapo at temperature of water
              f_s(1)=ro_vap_t_w(1)/(ro_vap_t_w(1)+ ro_ar_t_w(1))
              // f_s(1) : specific humidoty of saturated vapor at temperatyre of water
          
              H_g(1)=hg_in // H_g(1) enthalpy of humid air at air inlet temperature 
              f_g(1)=hum_vap_t(1)// f_g(1): specific humidity of air at inlet temoerature
              
              dfg(1)= m_dot_w/m_dot_ar*cpw*Dt_w/(H_s(1)-H_g(1))/(1-f_s(1))*(f_s(1)-f_g(1))
                                  // dfg : the increment of specific humidity ofb air at exit  of division
                                  
                                  
              D_ip= lamda * (m_dot_w / m_dot_ar) ^ (-np)  // D_ip parameter of merkel equation
              Dz(1)=cpw*Dt_w/(H_s(1)-H_g(1))/D_ip         // height of packing division
            
              H_g(2)=H_g(1)+ Qtot_div(1)/m_dot_ar_new(1)   // enthalpy of humid air at the inlet of the second division
              f_g(2)=f_g(1)+dfg(1)                         // specific humidity at the entrance of the second division
             
             hum_vap_t(2)=f_g(2)         
              
              m_dot_ar_new(2)=(1+hum_vap_t(2))*m_dot_ar    //  mass rate of flow of air at INLET of the 2nd division kg/m^2/sec
              kst division kg/m^2/sec
              m_dot_w_new(2)=m_dot_w+hum_vap_t(2)*m_dot_ar  // mass rate of flow of water at the exit of the 2nd division
              
              //Qtot_div(2)=m_dot_w_new(2)*cpw*Dt_w         // energy transfer of 2nd division
              
              
              for k =2 :(kmax-1)
                  t_w(k)=t_w(k-1)+Dt_w
        
        
        t_k=t_w(k)+273.15
         T_k(k)=t_w(k)+273.15
        t_k=T_k(k)
         [ p ] = saturation(t_k );
         
      state=1
       [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
                ro_vap_t_w(k) =ro*1000           // density of saturated vapor at water temperature at exit kth fivision
                ro_ar_t_w(k) = (p_barm- p ) * 100000 / r_j / t_k; // density of dray air at warer temperature at exit kth division
               
                H_s(k)=(h1*ro_vap_t_w(k)+ro_ar_t_w(k)*cpar*t_w(k))/(ro_vap_t_w(k)+ro_ar_t_w(k)) 
                      // enthalpy of vapor at water temperature at exit

              f_s(k)=ro_vap_t_w(k)/(ro_vap_t_w(k)+ ro_ar_t_w(k))
                         // specific humidity at water temperature of kth division at exit 
              
             H_g(k)=H_g(k-1)+ Qtot_div(k-1)/m_dot_ar_new(k-1) 
             // enthalpy of humid air at air temperature at inlet of kth division
              dfg(k)= m_dot_w_new(k)/m_dot_ar_new(k)*cpw*Dt_w/(H_s(k)-H_g(k))/(1-f_s(k))*(f_s(k)-f_g(k))
                             // increment of specific humidity of kth division  
              d_Hs_Hg(k)=H_s(k)-H_g(k)
              if d_Hs_Hg(k)< 0 then
                disp  'PROCESSING STOPPED DUE TO INOPUT DATA NOT SUITABLE- H_S IS LESS THAN h_G  '
                  break
                  end
              
               
             
              D_ip= lamda * (m_dot_w_new(k-1) / m_dot_ar_new(k-1)) ^ (-np)
              Dz(k)=cpw*Dt_w/(H_s(k)-H_g(k))/D_ip  // height of kth division
              
              
             
              
              f_g(k+1)=f_g(k)+dfg(k) // specific humidity of air at air temperature at inlet of (k+1) division
             
             
             
             
              
            
             
             hum_vap_t(k+1)=f_g(k+1)/(1-f_g(k+1))// humidity of air at inlet of (k+1) th division
              
              
              
              
              
              
              m_dot_ar_new(k+1)=(1+hum_vap_t(k+1))*m_dot_ar
              // mass flow rate of humid air at inlet of (k+1)divisiom until k=kmax-2
              
              m_dot_w_new(k+1)=m_dot_w+(hum_vap_t(k+1)*m_dot_ar)
               // mass flow rate of water at exit of (k+1)divisiom until k=kmax-2
                Qtot_div(k)=m_dot_w_new(k)*cpw*Dt_w
              end
              
             
            
             t_w_inlet  = t_w(kmax-1)+Dt_w
            
           
              
           H_g_exit=H_g(kmax-1)+ Qtot_div(kmax-1)/m_dot_ar_new(kmax)
            
             
   m_dot_ar_new_exit= m_dot_ar_new(kmax) // exit humid air mass rate of flow 
   m_dot_w_new_inlet =m_dot_w+(hum_vap_t(kmax)*m_dot_ar)// inlet water mass of flow 
              
              
           
             
             Z=sum(Dz)
             
             ro_ar_in=( ro_vap_tar_in + ro_ar_tar_in_seco)
           v=m_dot_ar/(ro_ar_in)
           cons_dp=2*9.81/v^2*1000/ro_ar_in
           Np=3.401*m_dot_w^.59*m_dot_ar^-.542
           np=1000*Np/cons_dp*Z
          
               disp 'exit air specific humidity'
              disp (f_g(kmax))
               
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
             
             
             disp 'packing height meter' 
             disp (Z)
        
              Q_Ztot=sum(Qtot_div)
              disp 'heat transfer load kW/square meter'
              disp(Q_Ztot)
              
              t_ar_out=(H_g_exit-2501.5*f_g(kmax))/(cpar*(1-f_g(kmax))+1.83*f_g(kmax))
              disp 'exit air temperature'
              disp (t_ar_out)
              
               disp 'pressure loss mm water' 
           disp (np)
             
    endfunction
    


    
    
