clear all
detector_phase=[0];
main_ratio=1;

%for mainlooop=[-2:4]
for mainlooop=[0:1]
    %for mainlooop=[2]
    
    phi_an_orig=20;
    
    phase_in=0;
    dist_in_hz=0;
    erro_in_deg=0;
    add_text='';
    sw=5;
    dw=1/(1*sw);
    disp(['Dwell time = ' num2str(dw) '/n with n=1 or 2'])
    %t=0:dw:1000;
    %t=0:dw:10000-dw;
    ffac=1000;
    th=0:dw/(2*ffac):2+dw/(2*ffac);
    
    t=th(1:2:end);
    lb=1/5;%1/2
    listpos=1:ffac:size(t,2);
    
    
    %nu=[2.01321321 ];
    amp=[1 0.1 0.1 0.1 0.5];
    height=[0.99 0.9 0.5 ];
    height=[0.99 0.9 ];
    heightfirst=[ 0.5 0.1 0.02  ];
    heightfirstd=[ 0.1 0.02 ];
    
    %%sim td complex points
    % fid_sim= (main_ratio)*exp(j*2*pi*((t )*nu)).*exp(-t*2*pi*lb);
    %%seq td*2 real points
    fid_full=(0)*th;
    nu=[5 -5 -2 1 9.5];
    nu=[-1.13];
    nu=[-2];%from -sw/2 to sw/2
    nurref=-sw/2;
    if mainlooop==1
        ti_txt='Sim. y';
        
        nuref=[nurref]*0;
        fid_ref=fid_full+exp(0*pi/2+1i*2*pi*((th)*nuref));%.*exp(-th*2*pi*lb);
        fid_ref=fid_ref(1:2:end);% 1, 3, 5...
    end
    if mainlooop==0
        ti_txt='Sim. x';
        nuref=[nurref]*0;
        fid_ref=fid_full+exp(1i*(pi/2+2*pi*((th)*nuref)));%.*exp(-th*2*pi*lb);
        fid_ref=fid_ref(1:2:end);% 1, 3, 5...
    end
    
    if mainlooop==2
        ti_txt='Redfield ';
        
        nuref=[nurref]*1;
        fid_ref=fid_full+exp(1i*(1*pi/2+2*pi*((th)*nuref)));%.*exp(-th*2*pi*lb);
    end
    if mainlooop==3
        ti_txt='Redfield unwinded';
        
        nuref=[0]*0000000;
        nu=nu-[ nurref];
        fid_ref=fid_full+exp(1i*(1*pi/2+2*pi*((th)*nuref)));%.*exp(-th*2*pi*lb);
    end
    if mainlooop==4
        ti_txt='Single channel -0.125 Hz shift';
        
        nuref=-0.125;
        fid_ref=fid_full+exp(1i*(1*pi/2+2*pi*((th)*nuref)));%.*exp(-th*2*pi*lb);
    end
    if mainlooop==-1
        ti_txt='Sim y -0.125 Hz shift';
        
        nuref=-0.125;
        fid_ref=fid_full+exp(1i*(000*pi/2+2*pi*((th)*nuref)));%.*exp(-th*2*pi*lb);
        fid_ref=fid_ref(1:2:end);% 1, 3, 5...
        
    end
    if mainlooop==-2
        ti_txt='Sim x -0.125 Hz shift';
        
        nuref=-0.125;
        fid_ref=fid_full+exp(1i*(1*pi/2+2*pi*((th)*nuref)));%.*exp(-th*2*pi*lb);
        fid_ref=fid_ref(1:2:end);% 1, 3, 5...
        
    end
    nu_shifted=nu+sw/2;%list of frequ all positive for "qf"
    for loo_spi=1:size(nu,2)
        fid_full=fid_full+amp(1,loo_spi)*exp(1i*2*pi*((th)*nu(1,loo_spi))).*exp(-th*2*pi*lb);
    end
    fid_full=(fid_full*exp(1i*detector_phase/180*pi));
    fid_full_shfited=0*fid_full;
    for loo_spi=1:size(nu_shifted,2)
        fid_full_shfited=fid_full_shfited+amp(1,loo_spi)*exp(1i*2*pi*((th)*nu_shifted(1,loo_spi))).*exp(-th*2*pi*lb);
    end
    fid_full_shfited=real(fid_full_shfited*exp(1i*detector_phase/180*pi));
    
    fid_sim=fid_full(1:2:end);% 1, 3, 5...
    fid_seq= 0*fid_full;
    %  figure(1111);clf;plot(real(phase_in));dfasdf
    
    %mamip to get seq data algebra (phase 90)
    %     skip_detail=1;
    %     if skip_detail==0,
    %         for cursor=1:size(fid_full,2),
    %             if mod(cursor-1,4)==0,
    %                 fid_seq(cursor)=+real(fid_full(1,cursor));
    %             end
    %             if mod(cursor-1,4)==1,
    %                 fid_seq(cursor)=+imag(fid_full(1,cursor));
    %             end
    %             if mod(cursor-1,4)==2,
    %                 fid_seq(cursor)=-real(fid_full(1,cursor));
    %             end
    %             if mod(cursor-1,4)==3,
    %                 fid_seq(cursor)=-imag(fid_full(1,cursor));
    %             end
    %         end
    %     else
    %same as above, but using phase
    phi_an=(-90-phi_an_orig)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
    phases_pt=0:size(fid_seq,2)-1;
    phases=phi_an*phases_pt;
    fid_seq_shifted=real(fid_full.*exp(-1i*phases));
    
    phi_an=(-90-0)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
    phases_pt=0:size(fid_seq,2)-1;
    phases=phi_an*phases_pt;
    fid_seq=real(fid_full.*exp(-1i*phases));
    
    
    fid_seq_qf=real(fid_full_shfited);
    figure(2);clf;
    % plot3(t,imag(fid_sim),real(fid_sim),'b-','linewidth',1.5);hold on
    if mainlooop==0
        plot3(t,imag(fid_sim),-real(fid_sim),'r-','linewidth',1.5);hold on
    end
   if mainlooop==1
        plot3(t,-imag(fid_sim),real(fid_sim),'r-','linewidth',1.5);hold on
    end
  
plot3([min(t) max(t) ],[0 0],[0 0],'k-');hold on
if mainlooop>=2
    plot3(th,imag(fid_ref),real(fid_ref),'k:');hold on
    for_quad_pt=1:ffac/2:size(th,2)/2;
    fid_ref=fid_ref(1:2:end);
    
else
    plot3(t,imag(fid_ref),real(fid_ref),'k:');hold on
    for_quad_pt=1:ffac:size(th,2)/2;
    
end
%
% for for_quad_pt2=for_quad_pt
%     for_quad_pt2
%     plot3([for_quad_pt2],[,for_quad_pt2*0],[for_quad_pt2*0],'.')
%         plot3([t(for_quad_pt2) t(for_quad_pt2)],[0 imag(fid_sim(for_quad_pt2))],[0 real(fid_sim(for_quad_pt2))],'k-')
%
%     plot3([t(for_quad_pt2) t(for_quad_pt2)],[0 imag(fid_sim(for_quad_pt2))],[0 t(for_quad_pt2)*0],'g-')
%     plot3([t(for_quad_pt2) t(for_quad_pt2)],[0 t(for_quad_pt2)*0],[0 real(fid_sim(for_quad_pt2))],'r-')
% end
sp=imag(fid_ref()).*imag(fid_sim())+real(fid_ref()).*real(fid_sim());

plot3(t,sp.* imag(fid_ref),sp.* real(fid_ref),'k-','linewidth',1.5)

for for_quad_pt2=for_quad_pt
    plot3([t(for_quad_pt2) t(for_quad_pt2)],[-imag(fid_ref(for_quad_pt2))*000000 imag(fid_ref(for_quad_pt2))],[-real(fid_ref(for_quad_pt2))*000000 real(fid_ref(for_quad_pt2))],'k:')
    sp=imag(fid_ref(for_quad_pt2))*imag(fid_sim(for_quad_pt2))+real(fid_ref(for_quad_pt2))*real(fid_sim(for_quad_pt2));
    v1=imag(fid_ref(for_quad_pt2))*imag(fid_ref(for_quad_pt2))+real(fid_ref(for_quad_pt2))*real(fid_ref(for_quad_pt2));
    v2=imag(fid_sim(for_quad_pt2))*imag(fid_sim(for_quad_pt2))+real(fid_sim(for_quad_pt2))*real(fid_sim(for_quad_pt2));
    if sp<0 sii=-1;else sii=1;end
    %  sp=sii*sqrt(sp*sii)/(sqrt(v1)*sqrt(v2));
    sp=(sp)/((v1));
    
    angr=acos(sp)*sii;
    angi=cos(angr);
    anga=sin(angr);
    %%% plot3([t(for_quad_pt2) t(for_quad_pt2)],sp*[0 imag(fid_ref(for_quad_pt2))],sp*[0 real(fid_ref(for_quad_pt2))],'k-')
    if mainlooop==0
        
        plot3([t(for_quad_pt2) t(for_quad_pt2)],[imag(fid_sim(for_quad_pt2)) sp*imag(fid_ref(for_quad_pt2))],-[real(fid_sim(for_quad_pt2)) sp*real(fid_ref(for_quad_pt2))],'k:')
    end
    if mainlooop==1
        
        plot3([t(for_quad_pt2) t(for_quad_pt2)],-[imag(fid_sim(for_quad_pt2)) sp*imag(fid_ref(for_quad_pt2))],[real(fid_sim(for_quad_pt2)) sp*real(fid_ref(for_quad_pt2))],'k:')
    end
    % plot3([t(for_quad_pt2) t(for_quad_pt2)],anga*[0 imag(fid_ref(for_quad_pt2))],anga*[0 real(fid_ref(for_quad_pt2))],'r--')
    
    %
    %     plot3([t(for_quad_pt2) t(for_quad_pt2)],[0 imag(fid_sim(for_quad_pt2))],[0 t(for_quad_pt2)*0],'g-')
    %     plot3([t(for_quad_pt2) t(for_quad_pt2)],[0 t(for_quad_pt2)*0],[0 real(fid_sim(for_quad_pt2))],'r-')
    
end
view([-24 24])
%view([0 0])
axis([0 max(t) -1 1 -1 1]);
pbaspect([2 1 1])



%
%
% fid_sim(1)=fid_sim(1)/2;
% fid_seq(1)=fid_seq(1)/2;
% fid_seq_shifted(1)=fid_seq_shifted(1)/2;
% fid_seq_qf(1)=fid_seq_qf(1)/2;
%
%
%
% spectrum0=fftshift(fft((fid_sim)))*exp(-1i*detector_phase/180*pi);
% spectrum1=(fft((fid_seq)))*exp(-1i*detector_phase/180*pi);
% spectrumqf=(fft((fid_seq_qf)))*exp(-1i*detector_phase/180*pi);
% spectrum1_shifted=(fft((fid_seq_shifted)))*exp(-1i*detector_phase/180*pi);
% %keep second part
% %spectrum1=spectrum1(:,size(spectrum1,2)/2+1:end);
%
% phi_an=2*(90-0)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
% phases_pt=0:size(fid_sim,2)-1;
% phases=phi_an*phases_pt;
% fid_simt=(fid_sim.*exp(-j*phases));
% spectrum2=(fft((fid_simt)))*exp(-1i*detector_phase/180*pi);
%
%
% phi_an=2*(90-phi_an_orig)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
% phases_pt=0:size(fid_sim,2)-1;
% phases=phi_an*phases_pt;
% fid_sim=(fid_sim.*exp(-j*phases));
% spectrum2_shifted=(fft((fid_sim)))*exp(-1i*detector_phase/180*pi);
%
%
% % top=max(max(abs(spectrum)));
% % spectrum=spectrum/top;
% incsw=sw/size(t,2);
%
% scale=-sw/2+incsw/2:incsw :sw/2-incsw/2;
% scale1=-sw/2+incsw/2:incsw :3*sw/2-incsw/2;
%
%
%
% where_projy=sw/2;
% where_projx=-1;
%
% %plot3(where_projy+0*scale,imag(spectrum),real(spectrum),'k-'); hold on
% off_fig0=-18;
% off_fig1=2*off_fig0;
% off_fig2=3*off_fig0;
% off_fig3=4*off_fig0;
% off_fig4=5*off_fig0;
% off_figm=off_fig0;
% spacce=off_fig0;
% off_fig0=0;
%
% cur_f=figure(1);clf;
%
% for loop=1:6,
%     if loop==2,
%         plot(10*[1 1],spacce*(loop-1)+[0 -spacce-3],'k-'); hold on
%     else
%         plot(0*[1 1],spacce*(loop-1)+[0 -spacce-3],'k-'); hold on
%
%     end
% end
% plot(scale,off_fig0+real(spectrum0),'k-'); hold on
%
% pos_x_text=-9.5;pos_y=15;
%
% text(pos_x_text,off_fig0+pos_y,'sim')
% plot(scale1,off_fig1+real(spectrum1),'k-'); hold on
% text(pos_x_text,off_fig1+pos_y,'seq')
%
% plot(scale1,off_figm+real(spectrumqf),'k-'); hold on
% text(pos_x_text,off_figm+pos_y,'QF')
%
% plot(scale1,off_fig3+real(spectrum1_shifted),'k-'); hold on
% text(pos_x_text,off_fig3+pos_y,'Shifted seq')
%
% plot(scale,off_fig2+real(spectrum2),'k-'); hold on
% text(pos_x_text,off_fig2+pos_y,'sim-TPPI')
%
%
% plot(scale,off_fig4+real(spectrum2_shifted),'k-'); hold on
% text(pos_x_text,off_fig4+pos_y,'Shifted sim-TPPI')


title(ti_txt,'interpreter','none')
set(gcf,'color','w');

set(gca,'ytick',[])
ylabel('Ix');
set(gca,'ztick',[])
zlabel('Iy');
xlabel('t [s]');

cur_f.Units='centimeters';
cur_f.Position=[30 30 16 10];
cur_f.Position=[30 30 8 5];


orient landscape
print('-depsc','-tiff','-r600',[ 'Fig_demo_red_'  num2str(mainlooop) '.eps']);%here

end
