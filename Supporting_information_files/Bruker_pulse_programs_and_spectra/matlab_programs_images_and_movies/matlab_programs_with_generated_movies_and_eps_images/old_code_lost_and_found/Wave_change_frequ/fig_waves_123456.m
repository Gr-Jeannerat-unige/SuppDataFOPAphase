clear all
larmor=1;
nu0=larmor;
figure(3);clf;

si_pu=1.5;
b1=1/(4*si_pu)

si_follow=6;
si_full=8;
si_before_pu=4;
t=si_before_pu;
time=-(si_full-si_follow):0.01:si_follow;
time=0:0.001:si_follow;
time2=si_before_pu:0.001:si_pu+si_before_pu;
time3=[0:0.001:si_before_pu ];
time4=[  si_pu+si_before_pu:0.001:si_follow];
for shifts=[0:1*b1:2*b1]
nu1=larmor+shifts;

phi=2*pi*((nu1-nu0)/1)*t;

 cur_f=figure(3);
 signal2=cos(2*pi*nu0*time2);
 signal=cos(2*pi*nu0*time);
 signal1=cos(2*pi*nu1*time);
        plot(time,signal,'k-');hold on      
        plot(time,signal1,'k--');hold on      
       % plot(time,signal1-signal,'r--');hold on   
        plot([min(time) max(time)],[0 0],'k-')
        off_ver=2.5;
        signal=cos(2*pi*nu0*time);
 signal22=cos(2*pi*nu1*time3-phi);
 signal222=cos(2*pi*nu1*time4-phi);
 signal11=cos(2*pi*nu1*time2-phi);
         plot(time,signal-off_ver,'k-','color',0.7*[1 1 1]);hold on      
         plot(time2,signal2-off_ver,'k-');hold on      
        plot(time3,signal22-off_ver,'k--','color',0.7*[1 1 1]);hold on      
        plot(time4,signal222-off_ver,'k--','color',0.7*[1 1 1]);hold on      
        plot(time2,signal11-off_ver,'k--');hold on      
       % plot(time,signal1-signal,'r--');hold on   
        plot([min(time) max(time)],[0 0]-off_ver,'k-')
        plot( si_before_pu*[1 1]       ,[-1 1.1]-off_ver,'k-')
        plot((si_before_pu+si_pu)*[1 1],[-1 1.1]-off_ver,'k-')
cur_f.Units='centimeters';
cur_f.Position=[20 6 12 6];
end

      orient landscape

print('-depsc','-tiff','-r600',[ 'Fig_wave123456.eps']);%here
