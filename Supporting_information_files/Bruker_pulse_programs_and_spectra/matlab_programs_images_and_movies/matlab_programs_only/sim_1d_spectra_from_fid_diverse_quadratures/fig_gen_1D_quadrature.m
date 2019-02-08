clear all
inc_fig=10;
for detector_phase=[0 45 90]
    for aliase=[0 1]
        for first_tp_non_zero=[0 2]%1 is for non corrected phase error
            main_ratio=1;
            
            mainlooop=0;
            
            phi_an_orig=20;
            
            phase_in=0;
            dist_in_hz=0;
            erro_in_deg=0;
            add_text='';
            sw=20;
            dw=1/sw;
            
            %t=0:dw:1000;
            %t=0:dw:10000-dw;
            th=0:dw/2:10000-dw/2;
            if first_tp_non_zero
                th=th+ dw/2;
            end
            t=th(1:2:end);
            lb=1/5;%1/2
            if aliase
                nu=[5 -5 -2 1 12.5];addtxt='aliasing_';addtxt2='aliased peak at 12.5 Hz ';
            else
                nu=[5 -5 -2 1 9.5];addtxt='';addtxt2='';
            end
            nu_shifted=nu+sw/2;%list of frequ all positive for "qf"
            %nu=[2.01321321 ];
            amp=[1 0.1 0.1 0.1 0.5];
            height=[0.99 0.9 0.5 ];
            height=[0.99 0.9 ];
            heightfirst=[ 0.5 0.1 0.02  ];
            heightfirstd=[ 0.1 0.02 ];
            
            %%sim td complex points
            % fid_sim= (main_ratio)*exp(j*2*pi*((t )*nu)).*exp(-t*2*pi*lb);
            %%seq td*2 real points
            
            %initialized
            fid_full        =zeros(size(th));%(0)*th;
            fid_full_shfited=zeros(size(th));%0*fid_full;

            %generation of the fid / loop over signals 
            for loo_spi=1:size(nu,2)
                fid_full        =fid_full        +amp(1,loo_spi)*exp(1i*2*pi*((th)*        nu(1,loo_spi))).*exp(-th*2*pi*lb);
                fid_full_shfited=fid_full_shfited+amp(1,loo_spi)*exp(1i*2*pi*((th)*nu_shifted(1,loo_spi))).*exp(-th*2*pi*lb);
            end
           
            %apply phase of detector
            fid_full        =(fid_full        *exp(1i*detector_phase/180*pi));
            fid_full_shfited=(fid_full_shfited*exp(1i*detector_phase/180*pi));

            fid_sim=fid_full(1:2:end);% 1, 3, 5...
% %             fid_seq= 0*fid_full;
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
            phases_pt=0:size(fid_full,2)-1;
            phases=phi_an*phases_pt;
            fid_seq_shifted=real(fid_full.*exp(-1i*phases));
            
            %phase correcton for sequ (redfield)
            phi_an=(-90-0)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
            phases_pt=0:size(fid_full,2)-1;
            phases=phi_an*phases_pt;
            fid_seq=real(fid_full.*exp(-1i*phases));
            
% % % %              %phase correcton for sequ (True sequ...)
% % % %              CH1=real(fid_full(1:2:end));
% % % %              CH2=    (fid_full(2:2:end));
% % % %             phi_an=(-90-0)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
% % % %             phases_pt=0:size(CH2,2)-1;
% % % %            phases_pt=phases_pt+0;
% % % %             phases=phi_an*phases_pt*1;
% % % %             CH2=imag(CH2.*exp(-1i*phases));
% % % %                         fid_sim=CH1*0+1*CH2;

            fid_seq_qf=real(fid_full_shfited);
% % %             figure(2);clf;
% % %             size(t)
% % %             
% % %             
% % %             plot3(t,real(fid_sim),imag(fid_sim))
% % %             axis([0 2 -1 1 -1 1])
% % %             view([-12 25])
         %   figure(3);clf
            
            if first_tp_non_zero==0
                fid_sim(1)=fid_sim(1)/2;
                fid_seq(1)=fid_seq(1)/2;
                fid_seq_shifted(1)=fid_seq_shifted(1)/2;
                fid_seq_qf(1)=fid_seq_qf(1)/2;
            end
            spectrum0=fftshift(fft((fid_sim)))*exp(-1i*detector_phase/180*pi);
            spectrum1=(fft((fid_seq)))*exp(-1i*detector_phase/180*pi);
            spectrumqf=(fft((fid_seq_qf)))*exp(-1i*detector_phase/180*pi);
            spectrum1_shifted=(fft((fid_seq_shifted)))*exp(-1i*detector_phase/180*pi);
            %keep second part
            %spectrum1=spectrum1(:,size(spectrum1,2)/2+1:end);
            
            phi_an=2*(90-0)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
            phases_pt=0:size(fid_sim,2)-1;
            phases=phi_an*phases_pt;
            fid_simt=(fid_sim.*exp(-1i*phases));
            spectrum2=(fft((fid_simt)))*exp(-1i*detector_phase/180*pi);
            
            
            phi_an=2*(90-phi_an_orig)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
            phases_pt=0:size(fid_sim,2)-1;
            phases=phi_an*phases_pt;
            fid_sim=(fid_sim.*exp(-1i*phases));
            spectrum2_shifted=(fft((fid_sim)))*exp(-1i*detector_phase/180*pi);
            
            if first_tp_non_zero>1
                ph1cor=-180;
                phc0cor=90;
                phi_an=1*(ph1cor)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
                phi_an=phi_an/size(spectrum0,2);
                phases_pt=0:size(spectrum0,2)-1;
                phases=phc0cor/180*pi+1*phi_an*phases_pt;
                
                spectrum0=spectrum0.*exp(1i*phases);
                spectrum2=spectrum2.*exp(1i*phases);
                spectrum2_shifted=spectrum2_shifted.*exp(1i*(phases+phi_an_orig/180*pi));
                
                
                ph1cor=-180*2;
                phc0cor=0;
                phi_an=1*(ph1cor)/180*pi;%normally for qseq, tppi, refield: phi_an=pi/2;
                phi_an=phi_an/size(spectrum1,2);
                phases_pt=0:size(spectrum1,2)-1;
                phases=phc0cor/180*pi+1*phi_an*phases_pt;
                
                
                spectrum1=spectrum1                .*exp(1i*(phases+pi/2));
                spectrum1_shifted=spectrum1_shifted.*exp(1i*(phases+pi/2));
                spectrumqf=spectrumqf.*exp(1i*phases);
            end
            
            % top=max(max(abs(spectrum)));
            % spectrum=spectrum/top;
            incsw=sw/size(t,2);
            
            scale=-sw/2+incsw/2:incsw :sw/2-incsw/2;
            scale1=-sw/2+incsw/2:incsw :3*sw/2-incsw/2;
            
            
            
            where_projy=sw/2;
            where_projx=-1;
            
            %plot3(where_projy+0*scale,imag(spectrum),real(spectrum),'k-'); hold on
            off_fig0=-18;
            off_fig1=2*off_fig0;
            off_fig2=3*off_fig0;
            off_fig3=4*off_fig0;
            off_fig4=5*off_fig0;
            off_figm=off_fig0;
            spacce=off_fig0;
            off_fig0=0;
            
            cur_f=figure(inc_fig);clf;
            for loop=1:6
                plot([-10 10],spacce*(loop-1)+[0 0],'k-'); hold on
                
                if loop==2
                    plot(10*[1 1],spacce*(loop-1)+[0 -spacce-3],'k-'); hold on
                else
                    plot(0*[1 1],spacce*(loop-1)+[0 -spacce-3],'k-'); hold on
                    
                end
            end
            
            linewidth=1.5;
            
            plot(scale,off_fig0+real(spectrum0),'k-','Linewidth',linewidth); hold on
            
            pos_x_text=-9.5;pos_y=15;
            
            text(pos_x_text,off_fig0+pos_y,'sim')
            plot(scale1,off_fig1+real(spectrum1),'k-','Linewidth',linewidth); hold on
            text(pos_x_text,off_fig1+pos_y,'seq')
            
            plot(scale1,off_figm+real(spectrumqf),'k-','Linewidth',linewidth); hold on
            text(pos_x_text,off_figm+pos_y,'QF')
            
            plot(scale1,off_fig3+real(spectrum1_shifted),'k-','Linewidth',linewidth); hold on
            text(pos_x_text,off_fig3+pos_y,'Shifted seq')
            %
            plot(scale,off_fig2+real(spectrum2),'k-','Linewidth',linewidth); hold on
            text(pos_x_text,off_fig2+pos_y,'sim-TPPI')
            
            
            plot(scale,off_fig4+real(spectrum2_shifted),'k-','Linewidth',linewidth); hold on
            text(pos_x_text,off_fig4+pos_y,'Shifted sim-TPPI')
            
            
            title([ addtxt2 'detector phase=' num2str(detector_phase)],'interpreter','none')
            if first_tp_non_zero
                title([ addtxt2 'detector phase=' num2str(detector_phase) ' First t1inc=0.5xSW'],'interpreter','none')
                
            end
            
            cur_f.Units='centimeters';
            tm=mod(inc_fig-10,4);
            cur_f.Position=[1+tm*16 1+10*(inc_fig-10-tm)/4 16 10];
            
            
            orient landscape
            print('-depsc','-tiff','-r600',[ 'Fig_demo_quad_1D_det_phase_' addtxt num2str(detector_phase) ' .eps']);%here
            if first_tp_non_zero
                print('-depsc','-tiff','-r600',[ 'Fig_demo_quad_1D_det_phase_' addtxt num2str(detector_phase) '_first_p1_non_zero' num2str(first_tp_non_zero) '.eps']);%here
            end
            inc_fig=inc_fig+1;
           % afsd
        end
    end
end
