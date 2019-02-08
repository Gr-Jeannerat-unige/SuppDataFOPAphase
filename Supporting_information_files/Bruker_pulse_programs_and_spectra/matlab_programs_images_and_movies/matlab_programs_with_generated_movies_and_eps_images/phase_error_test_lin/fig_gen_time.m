function fig_gen_time(factor,angle_deg,angle_view)
angle_deg
tim=[0 52];
if nargin==0
    factor=1;
    angle_deg=90;
end
figure(1)
pul_dur=10e-6;
full_dur=30e-6;
tau_corr=pul_dur*2/pi;
time_ser=0:0.01*1e-6:full_dur-(pul_dur-tau_corr);



angle_pulse=90/180*pi;%deg
ampli_hz=(angle_pulse/pul_dur)/(2*pi);
%disp(['pulse amplitude : ' num2str(ampli_hz) ' Hz'])



cur_f=figure(1);clf;
inca=1;

list_angl=(0:inca:360)*pi/180;
hold on
offsset_first_null=sqrt(15)/(4*pul_dur);
%disp(['offset_first_null : ' num2str(offsset_first_null) ' Hz ' ])
count_main=0;
start_pt=[0 -1 0];
start_pt_p=start_pt;
inc_store=1;
for loop_offset=factor*ampli_hz% not a loop anymore
    % if angle_deg>89.999
    %   plot(time_ser+(pul_dur-tau_corr),,'k-');hold on
    % plot(time_ser+(pul_dur-tau_corr),atan(sin(time_ser*loop_offset*2*pi)/cos(time_ser*loop_offset*2*pi)),'k--');hold on
    
    %  end
end
this_inc=0;

for loop_offset=factor*ampli_hz% not a loop anymore
    this_inc=this_inc+1;
    if count_main==0
        disp_on=1;
        count_main=count_main+1;
        
    else
        disp_on=0;
        
        count_main=count_main+1;
        if count_main==(5)
            count_main=0;
        end
    end
    nu_eff=sqrt(loop_offset*loop_offset+ampli_hz*ampli_hz);
    
    %plot field vector
    tilt_angle=atan((ampli_hz/loop_offset));
    if tilt_angle<0, tilt_angle=tilt_angle+pi;end
    % if disp_on
    % disp(['Offset : ' num2str(loop_offset) ' Hz w_eff=' num2str(nu_eff) ' Hz'])
    how_much_further=2;
    
    %   plot3([sin(tilt_angle)],[ 0],[ cos(tilt_angle)],'ko','MarkerSize',5)
    
    % end
    pos_mag=[0 0 1];
    pos_mag_p=pos_mag;
    field=[0 0 1];
    increment_tilt=pi/100000;
    inc=0;
    inc_sto=1;
    for til_tim=0:increment_tilt:angle_deg/180*pi
        di=cross([sin((tilt_angle)) 0 cos((tilt_angle)) ],pos_mag);
        di=di/norm(di);
        
        pos_mag=pos_mag+di*increment_tilt;
        if inc==0
            
            stor_tr(inc_sto,:)=pos_mag;inc_sto=inc_sto+1;
            pos_mag_p=pos_mag;
            inc=inc+1;
            
        else
            inc=inc+1;
            if inc==1000
                inc=0;
            end
        end
        
    end
    st1(this_inc)=power(power(stor_tr(end,1),2)+power(stor_tr(end,2),2),0.5);
    stx(this_inc)=stor_tr(end,1);
    sty(this_inc)=stor_tr(end,2);
    st2(this_inc)=180/pi*atan(stor_tr(end,1)./stor_tr(end,2));
    %  if disp_on
    %%%    plot3(stor_tr(:,1),stor_tr(:,2),stor_tr(:,3),'k:','linewidth',1.25)
    
    tim2=([1:size(stor_tr,1)]-1)/50*pul_dur;
    
    %  plot(tim2,stor_tr(:,1),'r--','linewidth',2,'color',[ 1 1 1]*0.5);hold on
 %   plot(tim2,180/pi*atan(stor_tr(:,1)/stor_tr(:,2)),'k','linewidth',2,'color',[ 1 1 1]*0.5)
    %  plot(tim2,stor_tr(:,3),'b:','linewidth',2,'color',[ 1 1 1]*0.5)
    axis([ 0 pul_dur  -180 180   ])
    
    % end
    
    %% plot3([0 pos_mag_p(1,1) ],[0 pos_mag_p(1,2) ],[0 pos_mag(1,3)])
    
    %  arrow([0 0 0 ],[pos_mag_p(1,1)  pos_mag_p(1,2)  pos_mag(1,3)],'linewidth',1.5)
    %  end
    drawnow
    start_pt_p=start_pt;
    start_pt=[pos_mag_p(1,1) pos_mag_p(1,2)  pos_mag(1,3)];
    store_traj(inc_store,:)=start_pt;inc_store=inc_store+1;
    %     figure(3);
    %     plot(180/pi*atan(stor_tr(:,1)./stor_tr(:,2)))
    %     print('-depsc','-tiff','-r600',[ 'lineraity_phase_' num2str(angle_deg) 'deg.eps']);%here
    %
    %     figure(1)
end
mi=min(min(factor));
ma=max(max(factor));
if size(factor,2)==1
    txt_ti=[' ' num2str(mi) ' x B1'];
else
    if mi==-ma
        txt_ti=[' +/- ' num2str(ma) ' x B1'];
        
    else
        txt_ti=[' (' num2str(mi) ':' num2str(ma) ') x B1'];
    end
end
%text(0.4,0.4,-0.9,[txt_ti '  ' num2str(angle_deg) ' deg.'])
%axis off
%if angle_deg==90

% print('-depsc','-tiff','-r600',[ 'Final_time_domain_sim.eps']);%here
figure(2);clf;
plot(factor,st2)

mid=((size(st2,2)-1)/2);
slope= (st2(mid+1)-st2(mid))/(factor(mid+1)-factor(mid));
slope1=(st2(mid+2)-st2(mid+1))/(factor(mid+2)-factor(mid+1));
tmp=(mid):(size(st1,2)-1);
tmpcc=0;
clear stoz;
for mid=tmp
    tmpcc=tmpcc+1;
   slope=(st2(mid+1)-st2(mid))/(factor(mid+1)-factor(mid));
        text(factor(mid),st2(mid),[ num2str(slope) ' ' num2str(factor(mid)) ])
stoz(tmpcc,1)=factor(mid);
stoz(tmpcc,2)=slope;
end
        title([' sp1=' num2str(slope) ' sp2=' num2str(slope1)  ])

        print('-depsc','-tiff','-r600',[ 'linearity_phase_' num2str(angle_deg) '_deg.eps']);%here
figure(22);clf;
size(factor)
size(st1)
size(st2)

plot(factor,st1,'k-');
[a b]=min(abs(st1-0.9));text(factor(b),st1(b),[num2str(factor(b)) ' '  num2str(st1(b)) ' ' num2str(st1(b+1))]);

hold on
plot(factor,stx,'k--')

plot(factor,sty,'k:')
plot(factor,sty.*sty,'c:')
tot=sty.*sty;
[a b]=min(abs(sty+0.9));text(factor(b),sty(b),[num2str(factor(b)) ' '  num2str(sty(b)) ' ' num2str(sty(b+1))]);
[a b]=min(abs(sty-0.9));text(factor(b),sty(b),[num2str(factor(b)) ' '  num2str(sty(b)) ' ' num2str(sty(b+1))]);
[a b]=min(abs(tot-0.9));text(factor(b),tot(b),[num2str(factor(b)) ' '  num2str(tot(b)) ' ' num2str(tot(b+1))]);

drawnow
        print('-depsc','-tiff','-r600',[ 'linearity_amplitude_' num2str(angle_deg) '_deg.eps']);%here
figure(3);clf;
plot(stoz(:,1),stoz(:,2))
ylim([-180 180])
        print('-depsc','-tiff','-r600',[ 'linearity_slope_phase_' num2str(angle_deg) 'deg_full.eps']);%here
ylim([-58 -40])
        print('-depsc','-tiff','-r600',[ 'linearity_slope_phase_' num2str(angle_deg) 'deg_zoom.eps']);%here

figure(1)
% end
end
