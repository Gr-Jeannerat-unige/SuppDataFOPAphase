clear all
pul_dur=10e-6;
angle_pulse=90/180*pi;%deg
ampli_hz=(angle_pulse/pul_dur)/(2*pi);
disp(['pulse amplitude : ' num2str(ampli_hz) ' Hz'])



inca=1;
list_angl=(0:inca:360)*pi/180;
hold on
offsset_first_null=sqrt(15)/(4*pul_dur);
disp(['offset_first_null : ' num2str(offsset_first_null) ' Hz ' ])
count_main=0;
start_pt=[0 -1 0];
start_pt_p=start_pt;
inc_store=1;
values=0:ampli_hz/50:ampli_hz*6;%for 12
tab_off=zeros(1,size(values,2));
for loop_offset=values
    
    nu_eff=sqrt(loop_offset*loop_offset+ampli_hz*ampli_hz);
       tab_off(inc_store)=nu_eff;

    %plot field vector
    tilt_angle=atan(ampli_hz/loop_offset)
    
    pos_mag=[0 0 1];
    pos_mag_p=pos_mag;
    field=[0 0 1];
    increment_tilt=pi/100000;
    inc=0;
    inc_sto=1;
    for til_tim=0:increment_tilt:pi/2
        di=cross([sin(tilt_angle) 0 cos(tilt_angle) ],pos_mag);
        di=di/norm(di);
        
        pos_mag=pos_mag+di*increment_tilt;
      
          
                %    plot3([pos_mag_p(1,1) pos_mag(1,1)],[pos_mag_p(1,2) pos_mag(1,2)],[pos_mag_p(1,3) pos_mag(1,3)],'r-')
        
            pos_mag_p=pos_mag;
            inc=inc+1;
            
      
        
    end
   
    start_pt_p=start_pt;
    
    
    
    start_pt=[pos_mag_p(1,1) pos_mag_p(1,2)  pos_mag(1,3)];
    store_traj(inc_store,1)=norm(start_pt(1,1:2));
    disp_part_mag=store_traj(inc_store,1)
    store_phe(inc_store,1)=atan(start_pt(1,1)/start_pt(1,2));
    
    inc_store=inc_store+1;
    % plot3([start_pt_p(1,1) start_pt(1,1)],[start_pt_p(1,2) start_pt(1,2)],[start_pt_p(1,3) start_pt(1,3)],'g-')
    
end
max(max(store_traj))
cur_f=figure(1);clf;
plot(values/ampli_hz,store_traj,'k-','linewidth',1.5);hold on
plot(-values/ampli_hz,store_traj,'k-','linewidth',1.5)
[a b]=min(abs(store_traj-0.9));txt=[num2str(values(1,b)/ampli_hz) 'x' num2str(ampli_hz) '=' num2str(values(1,b)) ' Hz']
title(['function of nu(eff)/w1 ' txt])
axis([-max(values/ampli_hz) max(values/ampli_hz) 0 1.0001])
plot(3*[-1 1],[0.9 0.9],'k:');


        
              
cur_f.Units='centimeters';
cur_f.Position=[1 6 4 4];
  

      orient landscape

print('-depsc','-tiff','-r600',[ 'Fig_demo_profile_offset.eps']);%here

figure(2);clf;
plot(values/ampli_hz,store_phe,'k-','linewidth',1.5);hold on
plot(-values/ampli_hz,store_phe,'k-','linewidth',1.5)
title('function of w eff/w1')
%print('-depsc','-tiff','-r600',[ 'Fig_demo_profile_phase_err.eps']);%here
