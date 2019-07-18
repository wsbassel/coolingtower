tic();
getd ("share");
getd ("TOWER")
disp 'tar_in is inlet air temperature ='
tar_in=30
disp (tar_in)
disp "hum rel_in is  inlet relative humidity ="
hum_rel_in=.2
disp (hum_rel_in)
disp 'inlet water temperature'
tw_in=39
disp(tw_in)
disp 'kmax=number of packing divisions'

kmax=11
disp (kmax)
disp 'apro =approach in degrees centigrades'
apro=5
disp (apro)
disp 'p_barm = atmospheric pressure'
p_barm=1
disp (p_barm)
disp 'lamda = packing parameter'
lamda=.341
disp (lamda)
disp ' np= packing parameter'
np=.57
disp (np)
disp ' m_dor_w water= mass rate of flow kg/square meter/sec'
m_dot_w  =6
disp (m_dot_w)
disp 'w_ar_ratio = ratio o mass rate of flow of water to mass rate of flow of air'
w_ar_ratio=1.6
disp (w_ar_ratio)
[Z]=tower_cool_humid_simple (w_ar_ratio ,m_dot_w  , np ,lamda , p_barm  ,apro ,kmax ,tw_in, tar_in , hum_rel_in )
t=toc()
disp 'CPU time seconds'
disp (t)
