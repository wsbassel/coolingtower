function [ wbt , hg ] = wet_bulb(t , hum_rel );
   
        t_k=t+273.15
      //  p_t=p;
        p_barm= 1;
        r_j = 287;
        r_j_w = 461.51;
        cpar = 1;
        [ p ] = saturation(t_k )
   p_vap = p
     p_t = p;
     
   state=1    
    [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
    
     ro_vap_t = ro*1000
    
         ro_ar_t = ((p_barm -p_vap) * 100000 / (r_j * (t + 273.15)));
          hum_ar_t = (ro_vap_t * hum_rel) / ro_ar_t;
           h_vap_t = h1                  //2501.6+ 1.827 * t ;
           hmist = hum_ar_t * h_vap_t + cpar * t;
        hg=hmist;
        
       
        
        t11=t;
        pt1(1)=p_t;
        
        if hum_rel< 1
           
               t_k=t11+273.15
        
        
        [ p ] = saturation(t_k )
   p_vap=p
   state=1     //vapor
    [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
     ro_vap_t = ro*1000
    
   
         ro_ar_t = ((p_barm-p_vap ) * 100000 / (r_j * (t11 + 273.15)));
          hum_ar_t = (ro_vap_t * hum_rel) / ro_ar_t;
           h_vap_t = h1                  //2501.6+ 1.827 * t ;
           hgl1 = hum_ar_t * h_vap_t + cpar * t11;
     
          end
          //hgl=zeros(10,1)
               for k =1: 100
              //------------------------------
                t1(k)=t11-.1*k;
                t_k=t1(k)+273.15;
                [ p ] = saturation(t_k );
               
                
                 [ p ] = saturation(t_k )
   p_vap=p
   state=1     //vapor
    [ ro, u1 , h1, s1] =wageeh_keenan (p, t_k, state)
     ro_vap_t = ro*1000
    
   
         ro_ar_t = ((p_barm -p_vap) * 100000 / (r_j * (t1(k) + 273.15)));
          hum_ar_t = (ro_vap_t ) / ro_ar_t;
           h_vap_t = h1                  //2501.6+ 1.827 * t ;
          hgl(k) = hum_ar_t * h_vap_t + cpar * t1(k);
               // ro_vap_t = p_t1(k) * 100000 / r_j_w / t_k;


              //  hmist = hum_ar_t * h_vap_t + cpar * t1(k);
dhgl(k)=hg-hgl(k)
 end   
 //disp (dhgl)        
[minv,mpos]=min(abs(dhgl))


wbt=t1(mpos)

 
               
   
     
    endfunction
    
    
//tic();
//getd ("C:\Users\wsb\Documents\Code\WetBulb\func\");
//t=40;
//hum_rel=.6
//[ wbt , hg ] = wet_bulb(t , hum_rel );
 

// disp (wbt);
// disp (hg);
// t = toc();
//disp(t);
