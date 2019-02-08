%demo phase
clear all
disp('you may have to enlage the main figure for proper rendering')
%set(groot,'defaultLineLineWidth',1);%less than 1 does not works
%  set(0,'defaultlinelinewidth',0.5);
%  set(0,'defaultaxeslinewidth',0.5);
%  set(0,'defaultpatchlinewidth',0.5);
% include in plot 1 the phase
ma_font_size=16;


thick_line_widht=2;
sw=20;
dw=1/sw;
factor_more_pt=1;
t=0:dw/factor_more_pt:100;%to
lb=1;
height=[0.99 0.9 sqrt(2)/2 0.5 ];%where shows width
height=[0.95 sqrt(2)/2 0.5  ];%where shows width and spiral lines
%height=[ sqrt(2) ];

heightfirst=[ 0.5 0.1 0.02  ];
heightfirst=[ 0.5   ];
heightfirstd=[ 0.1 0.02 ];
heightfirstd=[0.25];
%for anglet=[0 -2*360/22],
for list_of_cases=[13]% :12
    figure(200+list_of_cases);clf;hold on;       

   set_pos(ma_font_size);
    
    list_frequencies=[00 ];
    this_title=['Pure Lorentzian LW = ' num2str(lb) ];
    
    if list_of_cases==1
        list_frequencies=[2 ];
        
        
    end
   
    
    if list_of_cases>1 & list_of_cases<4,
        list_frequencies=[-2  2];
        this_title=[this_title ' J = 5 Hz'];
        
    end
    list_phases=[0 0];
     if list_of_cases==13,
        list_frequencies=[1 ];%2
        
         list_phases=[ -30];
        this_title=[this_title ' phase: 30 deg.'];
    end
    if list_of_cases==3,
        list_phases=[0 -30];
        this_title=[this_title ' phase: 30 deg.'];
    end
    
    if list_of_cases==4,
        list_frequencies=[-3 -2 2 3];
        list_phases=0*list_frequencies;
        this_title=[this_title ' J = 5, 1 Hz'];
        
    end
    if list_of_cases==5,
        J_coup=0.6;
        % list_frequencies=[-2.5-J_coup -2.5+J_coup 2.5-J_coup 2.5+J_coup];
        list_frequencies=[-J_coup/2 J_coup/2];
        list_phases=0*list_frequencies;
        
        this_title=[this_title ' d: J = ' num2str(J_coup) ' Hz'];
    end
    if list_of_cases==6,
        J_coup=0.6;
        % list_frequencies=[-2.5-J_coup -2.5+J_coup 2.5-J_coup 2.5+J_coup];
        list_frequencies=[-J_coup 0 0 J_coup];
        list_phases=0*list_frequencies;
        this_title=[this_title ' t: J = ' num2str(J_coup) ' Hz'];
        
    end
    if list_of_cases<7,
        this_title=[this_title ' (no window function)'];
        
    end
    window=1+0*t;
    if list_of_cases==7,
        line_br=1;
        this_title=[this_title 'Exp. mult LB = ' num2str(line_br) ' Hz'];
        window=exp(-t*2*pi*line_br/2);
        plot(t,window,'k-','DisplayName','Window fn.')
        
    end
    if list_of_cases==8,%lorentz and gauss
        line_br=0;
        window=exp(-t*2*pi*line_br/2);
        %sigma=2/2.355;mean=0;
        lw_g=1;
        sigma=0.5/(sqrt(2*log(2))*lw_g);mean=0;
        gauss=1/sqrt(2*pi)/sigma*exp(-(t-mean).^2/2/sigma/sigma);
        %      figure(555);clf;plot(t,gauss);
        %      asdfafd
        window=window.*gauss;
        this_title=[this_title 'Exp. mult LB = ' num2str(line_br) ' Hz / Gauss FWHM ' num2str(lw_g)];
        plot(t,window,'k-','DisplayName','Window fn.')
    end
    
    if list_of_cases==9,%lorentz and gauss
        line_br=-1.2;
        window=exp(-t*2*pi*line_br/2);
        %sigma=2/2.355;mean=0;
        lw_g=0.8;
        sigma=0.5/(sqrt(2*log(2))*lw_g);mean=0;
        gauss=1/sqrt(2*pi)/sigma*exp(-(t-mean).^2/2/sigma/sigma);
        %      figure(555);clf;plot(t,gauss);
        %      asdfafd
        %figure(555);clf;plot(t,window);
        % plot(t,window,'k.')
       % plot(t,gauss,'k:','DisplayName','Gauss')
        
        window=window.*gauss;
        plot(t,window/max(max(window)),'k-','DisplayName','Window fn.')
        
        this_title=[this_title 'Exp. mult LB = ' num2str(line_br) ' Hz / Gauss FWHM ' num2str(lw_g)];
        
    end
    aq_cos=.5;
    
    if list_of_cases==11,%cos
        power_cos=1;
        factor_cos=(pi/2)/aq_cos;
        window=power(cos(t*factor_cos),power_cos);
        window=window.*(t<aq_cos);
        plot(t,window,'k-','DisplayName','Window fn.')
        this_title=[this_title 'cos power ' num2str(power_cos) '. Truncated after t = ' num2str(aq_cos) ' s'];
        
        if power_cos==1;
            this_title=[this_title 'cos(0 to pi/2). Truncated after t = ' num2str(aq_cos) ' s'];
        end
        if power_cos==2;
            this_title=[this_title 'cos(0 to pi/2) sqared. Truncated after t = ' num2str(aq_cos) ' s'];
        end
    end
    
    if list_of_cases==12,%cos square
        power_cos=2;
        factor_cos=(pi/2)/aq_cos;
        window=power(cos(t*factor_cos),power_cos);
        window=window.*(t<aq_cos);
        this_title=[this_title 'cos power ' num2str(power_cos) '. Truncated after t = ' num2str(aq_cos) ' s'];
        
        if power_cos==1;
            this_title=[this_title 'cos(0 to pi/2). Truncated after t = ' num2str(aq_cos) ' s'];
        end
        if power_cos==2;
            this_title=[this_title 'cos(0 to pi/2) sqared. Truncated after t = ' num2str(aq_cos) ' s'];
        end
        
        plot(t,window,'k-','DisplayName','Window fn.')
        
    end
    if list_of_cases==10,%truncation
        aq_trunc=aq_cos;
        window=window.*(t<aq_trunc);
        this_title=[this_title 'Truncated after t = ' num2str(aq_cos) ' s'];
        %plot(t,window,'k-')
        
    end
    %     anglet =0;
    %
    %     if list_of_cases==2,
    %         anglet =-2*360/22;
    %     end
    
    
    % calculate shape
    fid=0*t;%initialize fid
    for loop_list_signals=1:size(list_frequencies,2),
        anglet=0*90+list_phases(1,loop_list_signals);
        phase=cos(anglet/360*2*pi)+j*sin(anglet/360*2*pi);%as complex number
        phase=phase/(abs(phase));
        %fid=(phase)*exp(j*2*pi*((t)*nu)).*exp(-t*2*pi*lb);
        freq_mod=exp(j*2*pi*list_frequencies(1,loop_list_signals)*t);
        
        % figure(111);clf;plot(real(freq_mod),'g-');;hold on; plot(imag(freq_mod),'r-') ;   figure(list_of_cases*3-2)
        
        fid=fid+(phase)*exp(-t*2*pi*lb/2).*freq_mod;
        %fid=fid+(phase)*exp(-t*2*pi*lb);
    end
    
    plot(t,real(fid),'k--','DisplayName','Crude FID');
    
    fid=fid.*window;
    %
    %
    %
    %fid=fid-offset;
    plot(t,real(fid),'k-','LineWidth',thick_line_widht,'DisplayName','After proc.');
    max_showed_in_sec=1;
    xlim([0 max_showed_in_sec]);
    % ylim([0 1.2]);
    other_method=1;
    if other_method,
        offset=0;
        fid(1)=fid(1)/2;
    else
        offset=sum(sum(abs(fid)))*(1/2*pi)*(1/sw);%%% important !!
    end
    spectrum=fftshift(fft((fid)));
    take=size(spectrum,2)/factor_more_pt;middle=(size(spectrum,2)/2);
    %     middle-take/2+1
    %     middle+take/2
    %     size(spectrum)
    spectrum=spectrum(1,round(middle-take/2+1:middle+take/2));
    % spectrum=conv(spectrum,window,'same');
    spectrum=spectrum-offset;
    top=max(max(abs(spectrum)));
    spectrum=spectrum/top;
    incsw=sw/size(spectrum,2);
    
    scale=-sw/2+incsw/2:incsw :sw/2-incsw/2;
    title(this_title)
    plot([0 max_showed_in_sec ],[0 0],'k-')
    
  %  legend('show')
    set_pos(ma_font_size);
          subplot_0=gca;

%    print('-depsc','-tiff','-r2400',[ './fig_phase_demo0_' num2str(list_of_cases) '.eps']);
    figure(list_of_cases*3-2)
    
    clf
set_pos(ma_font_size);
    plot([-sw/2 sw/2],[0 0 ],'k-');hold on
    plot(scale,abs(spectrum),'k-','LineWidth',thick_line_widht);hold on
    plot(scale,real(spectrum),'k-');
    plot(scale,imag(spectrum),'k--');
    
    %% draw lines and comments on line with
    %for real part
    for loop=1:size(heightfirst,2),%plot line width at different heights
        what=real(spectrum);
        tmp=(what(1,1:size(what,2)-1)-heightfirst(1,loop)).*(what(1,2:size(what,2)-0)-heightfirst(1,loop));
        list_i= find(tmp < 0);
        inc=0;a=0;b=0;ii=i;
        for i=list_i,
            
            a=scale(1,i)*0.5+0.5*scale(1,i+1);
            b=real(what(1,i))*0.5+0.5*real(what(1,i+1));
            inc=inc+1;
            if inc==2;
                % plot([a aa] ,[ b bb],'k-') ;
                plot([a aa] ,[ b b],'k-') ;
                plot([a ] ,[ b ],'kx') ;
                plot([ aa] ,[  b],'kx') ;
                text_hz=[num2str( abs( scale(1,i)-scale(1,ii)),'%.2f') ' Hz'];
                text(aa-4,bb,['Abs.:' num2str(heightfirst(1,loop)*100,'%.0f') '% ' text_hz],'FontSize',ma_font_size);
                store_width_from=aa;
                store_width_to=a;
                
            end
            aa=a;bb=b;ii=i;
        end
    end
    %same for imag part
    for loop=1:size(heightfirstd,2),%plot line disp
        what=imag(spectrum);
        tmp=(what(1,1:size(what,2)-1)-heightfirstd(1,loop)).*(what(1,2:size(what,2)-0)-heightfirstd(1,loop));
        list_i= find(tmp < 0);
        inc=0;a=0;b=0;ii=i;
        for i=list_i,
            
            a=scale(1,i)*0.5+0.5*scale(1,i+1);
            b=real(what(1,i))*0.5+0.5*real(what(1,i+1));
            inc=inc+1;
            if inc==2;
                %  plot([a aa] ,[ b bb],'k--') ;
                plot([a aa] ,[ b b],'k:') ;
                
                text(aa-4,bb,['Disp.:' num2str(heightfirstd(1,loop)*100,'%.0f') '% ' num2str( abs( scale(1,i)-scale(1,ii)),'%.1f') ' Hz'],'FontSize',ma_font_size);
            end
            aa=a;bb=b;ii=i;
            
        end
    end
    %same for abs part
    for loop=1:size(heightfirst,2),%plot line disp
        what=abs(spectrum);
        tmp=(what(1,1:size(what,2)-1)-heightfirst(1,loop)).*(what(1,2:size(what,2)-0)-heightfirst(1,loop));
        list_i= find(tmp < 0);
        inc=0;a=0;b=0;ii=i;
        for i=list_i,
            
            a=scale(1,i)*0.5+0.5*scale(1,i+1);
            b=real(what(1,i))*0.5+0.5*real(what(1,i+1));
            inc=inc+1;
            if inc==2;
                %  plot([a aa] ,[ b bb],'k--') ;
                plot([a aa] ,[ b b],'k--','LineWidth',thick_line_widht) ;
                text(a+0.5,bb,['Magn.: ' num2str( abs( scale(1,i)-scale(1,ii)),'%.1f') ' Hz'],'FontSize',ma_font_size);
            end
            aa=a;bb=b;ii=i;
            
        end
    end
    title(this_title)
     set_pos(ma_font_size);
          subplot_1=gca;

  %  print('-depsc','-tiff','-r2400',[ './fig_phase_demo1_' num2str(list_of_cases) '.eps']);
    
    %plot phase diagram
    figure(100+list_of_cases);clf
set_pos(ma_font_size); 
    min_for_bold=0.25;
    % %      data_to_plot=( 180+mod( angle(spectrum-0.5)/pi*180-180,360)-360);
    % %     plot(scale,data_to_plot,'g-');hold on
    data_to_plot=( 180+mod( angle(spectrum)/pi*180-180,360)-360);
    plot(scale,data_to_plot,'k--');hold on
    
    list_from  =find(  abs(spectrum(1,1:end)<=min_for_bold).*[abs(spectrum(1,2:end)>min_for_bold) 0]);
    list_to=find(  abs(spectrum(1,1:end)>=min_for_bold).*[abs(spectrum(1,2:end)<min_for_bold) 1]);
    size(list_to)
    min(size(list_from,2),size(list_to,2))
    for loop_bold_plot=1:(min(size(list_from,2),size(list_to,2))),
        tmp_dat=data_to_plot(1,list_from(1,loop_bold_plot):list_to(1,loop_bold_plot));
        tmp_sca=scale(1,list_from(1,loop_bold_plot):list_to(1,loop_bold_plot));
        plot(tmp_sca,data_to_plot(1,list_from(1,loop_bold_plot):list_to(1,loop_bold_plot)),'k-','LineWidth',thick_line_widht);
        
        
        %     vt1=max(tmp_dat);
        %     vt2=min(tmp_dat);
        if loop_bold_plot==1,
            [del tmpu1]=min(abs(tmp_sca-   store_width_to));
            vt2=tmp_dat(1,tmpu1);
            [del tmpu2]=min(abs(tmp_sca-   store_width_from));
            vt1=tmp_dat(1,tmpu2);
           % plot(store_width_from*[1 1],[vt1 vt2],'k--');
            plot(store_width_to*[1 1],[vt1 vt2],'k:');
            plot(store_width_from*[1 ],[vt1 ],'kx');
            plot(store_width_to*[1 ],[vt1 ],'kx');
            plot([store_width_from store_width_to],[vt1 vt1],'k-');
            text(store_width_to,vt1,[' ' text_hz ' half of disp.'],'FontSize',ma_font_size);
            tmda=-(vt1 - vt2)/(store_width_to-store_width_from);
            text(store_width_to,vt1/2,['  Slope: ' num2str(tmda,2) ' deg/Hz'],'FontSize',ma_font_size);
            %  diagonal
            plot([store_width_from store_width_to],[vt1  vt2],'k:');
            
        end
    end
    %position theoretical phases
    
    for loop_list_signals=1:size(list_frequencies,2),
                [del tmpu1]=min(abs(scale-   list_frequencies(1,loop_list_signals)));

        plot(scale(1,tmpu1)*[1 1], [ data_to_plot(1,tmpu1) list_phases(1,loop_list_signals) ] ,'k-')
                mid_pos=0.5*( data_to_plot(1,tmpu1) +list_phases(1,loop_list_signals));
             val_heret=round(   (data_to_plot(1,tmpu1) -list_phases(1,loop_list_signals) ),1);
text(scale(1,tmpu1),mid_pos,[' Error: ' num2str(val_heret)  '  deg.'],'FontSize',ma_font_size)
        plot(list_frequencies(1,loop_list_signals), list_phases(1,loop_list_signals)  ,'ko')
        plot(scale(1,tmpu1),data_to_plot(1,tmpu1) ,'k+')
    end
    
    
    title(['Phase in deg. (Broad line where S > ' num2str(min_for_bold) ')']);
    plot([min(scale) max(scale) ],[0 0],'k-')
    
     set_pos(ma_font_size);
     subplot_4=gca;

    %     print('-depsc','-tiff','-r2400',[ './fig_phase_demo4_' num2str(list_of_cases) '.eps']);

    figure(list_of_cases*3-1)
    
    clf
  %  set_pos(ma_font_size); 
    where_projy=sw/2;
    where_projx=-1;
       extra=sw/20;

    % plot3(where_projy+0*scale,imag(spectrum),real(spectrum),'k-'); hold on
    plot3(extra+where_projy+0*scale,imag(spectrum),real(spectrum),'b-','LineWidth',thick_line_widht); hold on
    plot3(scale,1+0*imag(spectrum),real(spectrum),'b-','LineWidth',thick_line_widht); hold on
    plot3(scale,imag(spectrum),real(spectrum),'r-','LineWidth',thick_line_widht); hold on
   % plot3(scale,0*imag(spectrum),1*real(spectrum),'LineWidth',thick_line_widht); hold on
    plot3(scale,1*imag(spectrum),0*real(spectrum),'k-'); hold on
    plot3(sw/2*[-1 1],[ 0 0],[ 0 0],'k-')
    plot3(sw/2*[-1 1],1+[ 0 0],[ 0 0],'k-')
  %  plot3(sw/2*[1 1],[ 1 -1],[ 0 0],'k-')
    % set(gca,'YDir','Reverse')
    
    %axis([ where_projy*[-1 1] -1 1 0 1])
    %% draw text and phase differences
    for loop=1:size(height,2),
        what=abs(spectrum);
        tmp=(what(1,1:size(what,2)-1)-height(1,loop)).*(what(1,2:size(what,2)-0)-height(1,loop));
        list_i= find(tmp < 0);
        [a b]=max(what);
        [a c]=max(abs(spectrum));
        list_i=[list_i b c];
        for i=list_i,
            plot3(scale(1,i)*[1 1],[ imag(spectrum(1,i)) 0],[ real(spectrum(1,i)) 0],'k:')
            
            av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))   *0.5 + 0.5* angle(real(spectrum(1,i+1))+j*imag(spectrum(1,i+1)));
            delta=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  - angle(real(spectrum(1,i+1))+j*imag(spectrum(1,i+1)));
            %text(where_projy*[1],[ imag(spectrum(1,i)) ],[ real(spectrum(1,i)) ],[num2str(360/(2*pi)*av_al,'%.1f') '(' num2str(360/(2*pi)*delta,'%.1f')  ')']);
            addtext=[num2str(height(1,loop)*100,'%.0f') '% '] ;
            signal_me=0;
            if i==c,
                addtext='max(abs) ';
                av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  ;
                signal_me=1;
            end
            if i==b,
                addtext='max(real) ';
                av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  ;
                
            end
            if i==b && i==c,
                addtext='max(real & abs) ';
                av_al=angle(real(spectrum(1,i))+j*imag(spectrum(1,i)))  ;
                signal_me=1;
            end
            if signal_me,
                %        plot3(where_projy*[1 1],[ imag(spectrum(1,i)) 0],[ real(spectrum(1,i)) 0],'k:','LineWidth',thick_line_widht)
            else
                %       plot3(where_projy*[1 1],[ imag(spectrum(1,i)) 0],[ real(spectrum(1,i)) 0],'k:')
            end
            
            %   text(where_projy*[1],[ imag(spectrum(1,i)) ],[ real(spectrum(1,i)) ],[addtext num2str(360/(2*pi)*av_al,'%.1f') ]);
            %  text(scale(1,i)*[1],[ imag(spectrum(1,i)) ],[ real(spectrum(1,i)) ],[addtext num2str(360/(2*pi)*av_al,'%.1f') ]);
        end
    end
    plot3(extra+[sw/2 sw/2],[0 0 ],[0 1.05 ],'k');
  %  plot3([-sw/2 -sw/2],[-1 1 ],[0 0 ],'k');
  %  plot3([-sw/2 sw/2],[-1 -1 ],[0 0 ],'k');
   % plot3(cos(0:pi/100:2*pi)/2*0+sw/2,cos(0:pi/100:2*pi)/2,sin(0:pi/100:2*pi)/2+0.5,'k--')
       pos_ref=1;

    plot3(extra+pos_ref*cos(0:pi/100:2*pi)/2*0+sw/2,pos_ref*cos(0:pi/100:2*pi),pos_ref*sin(0:pi/100:2*pi),'k-');%,'LineWidth',thick_line_widht)
    sta=30;Len=pos_ref+[0 0.05];Lenq=[0 min(min(Len))];
    wl=max(max(Len))+0.1;
    pha=pi/2;
    for lo_a=0:sta:360-sta
        plot3(extra+[1 1]*sw/2,Len*cos(pha+lo_a/180*pi),Len*sin(pha+lo_a/180*pi),'k-');%,'LineWidth',thick_line_widht)
        plot3(extra+[1 1]*sw/2,Lenq*cos(pha+lo_a/180*pi),Lenq*sin(pha+lo_a/180*pi),'k:');%,'LineWidth',thick_line_widht)
        text(extra+sw/2,wl*cos(pha+lo_a/180*pi),wl*sin(pha+lo_a/180*pi),num2str(lo_a),'HorizontalAlignment','center','VerticalAlignment','middle')
        if lo_a==30
            Len2=[0 1];
            plot3(extra+[1 1]*sw/2,Len2*cos(pha+lo_a/180*pi),Len2*sin(pha+lo_a/180*pi),'k-')
            
        end
    end
    %view([-24,22])
    view([-24,22])
    view([-20,15])
    %%%% set_pos(ma_font_size);
    subplot_2=gca;
    axis off
    
    this_ax=gca;
pos_ax=[1 1 25 20];set(this_ax,'Unit','centimeters');set(this_ax,'OuterPosition',pos_ax);
    print('-painters','-depsc2',[ './dispa.eps']);
    %    print('-depsc','-tiff','-r2400',[ './fig_phase_demo2_' num2str(list_of_cases) '.eps']);
    
    %%
    figure(list_of_cases*3);clf;
set_pos(ma_font_size);
   % plot(cos(0:pi/100:2*pi)/2,sin(0:pi/100:2*pi)/2+0.5,'k--')
    plot(cos(0:pi/100:2*pi),sin(0:pi/100:2*pi),'k--')
    hold on
    
    plot(imag(spectrum),real(spectrum),'k-','LineWidth',thick_line_widht)
    plot( imag(spectrum(1,1)),real(spectrum(1,1)),'ko');
    %list of extrema:
    tmp=abs(spectrum);
    tmp=[0 (tmp(3:end)<tmp(2:end-1)).*(tmp(1:end-2)<tmp(2:end-1)) 0];
    tmp=find(tmp);
    plot(imag(spectrum(1,tmp)),real(spectrum(1,tmp)),'kx');
    for loop_over_extrema=1:size(tmp,2),
        text(imag(spectrum(1,tmp(1,loop_over_extrema))),real(spectrum(1,tmp(1,loop_over_extrema))),[num2str(scale(1,tmp(1,loop_over_extrema))) ' Hz'],'FontSize',ma_font_size);
    end
    
    plot([0 0], [min(real(spectrum)) max(real(spectrum))],'k-');
    plot( [min(imag(spectrum)) max(imag(spectrum))],[0 0],'k-');
    axis equal
    set_pos(ma_font_size);
box off
subplot_3=gca;
%    print('-depsc','-tiff','-r2400',[ './fig_phase_demo3_' num2str(list_of_cases) '.eps']);
    
    %% assemble into full_figure
   if (list_of_cases==1) || (list_of_cases==13),
   
   main_fig=figure(1001);clf;
       set(main_fig,'PaperType','a4');%20x29 cm
         set(main_fig,'PaperUnit','centimeter');
         set(main_fig,'PaperOrientation','portrait');
     set(main_fig,'PaperPosition',[0 0 28 40]);% to adjust..
     set(main_fig,'Unit','centimeter');
     set(main_fig,'Position',[0 0 28 40]);
     drawnow
%     axis off
f=0.65;
   width=3.5*2*f;
   space_h=0.6*2*f;
   height=4.0*2*f;
   space_v=0.5*2*f;
   test_position='Position';
  % test_position='Position';
%             coorr=[0 0];
%             pos=[(height+space_h)*(coorr(1,1)) (height+space_v)*(coorr(1,2)) width height]
%           cur_obj1=  copyobj(subplot_0,main_fig); set(cur_obj1,'Unit','centimeters');set(cur_obj1,test_position,pos);set(main_fig, 'currentaxes', cur_obj1);drawnow;
        
            coorr=[0 1];
            pos=[2+(height+space_h)*(coorr(1,1)) 2+(height+space_v)*(coorr(1,2)) width+5.5 height]
          cur_obj2=  copyobj(subplot_1,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
        
            coorr=[0 0];
            pos=[2+(height+space_h)*(coorr(1,1)) 2+(height+space_v)*(coorr(1,2)) width+5.5 height]
          cur_obj2=  copyobj(subplot_2,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
        
            coorr=[1 0];
            pos=[2+(height+space_h)*(coorr(1,1))+5.5+0.5 2+(height+space_v)*(coorr(1,2))+0.5 width-0.5 height-0.5]
          cur_obj2=  copyobj(subplot_3,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
          
            coorr=[1 1];
            pos=[2+(height+space_h)*(coorr(1,1))+5.5 2+(height+space_v)*(coorr(1,2)) width height]
          cur_obj2=  copyobj(subplot_4,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
          
          
  %         print('-depsc','-tiff','-r2400',[ './fig_phase_demo_assembled_' num2str(list_of_cases) '.eps']);
 
    end
    if list_of_cases>=7 & list_of_cases<13 ,
%         if list_of_cases==7
%             main_fig=figure(1001);clf;
%         else
%              main_fig=figure(1001);
%         end
          main_fig=figure(1001);clf;
          coorr=[0 1];
%         if list_of_cases==10,coorr=[0 0];end
%         if list_of_cases==7,coorr=[0 1];end
%         if list_of_cases==8,coorr=[0 2];end
%         if list_of_cases==9,coorr=[0 3];end
%         if list_of_cases==11,coorr=[0 4];end
%         if list_of_cases==12,coorr=[0 5];end
        width=3.2*2;
        space_h=0.8*2;
        height=3.0*2;
        space_v=0.5*2;
        test_position='Position';
        % test_position='Position';
        coorr(1,1)=0;
        pos=[1.9+(width+space_h)*(coorr(1,1)) (height+space_v)*(coorr(1,2)) width-2 height];
        cur_obj1=  copyobj(subplot_0,main_fig); set(cur_obj1,'Unit','centimeters');set(cur_obj1,test_position,pos);set(main_fig, 'currentaxes', cur_obj1);drawnow;
        
        coorr(1,1)=1;
        pos=[2+(width+space_h)*(coorr(1,1))-2 (height+space_v)*(coorr(1,2)) width height];
        cur_obj2=  copyobj(subplot_1,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
%         
%         coorr(1,1)=2;
%         pos=[2+(width+space_h)*(coorr(1,1)) (height+space_v)*(coorr(1,2)) width+2 height];
%         cur_obj2=  copyobj(subplot_2,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
%         
        coorr(1,1)=2;
        pos=[2+(width+space_h)*(coorr(1,1))-2 (height+space_v)*(coorr(1,2))+0.3 width-2 height];
        cur_obj2=  copyobj(subplot_3,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
        
        coorr(1,1)=3;
        pos=[2+(width+space_h)*(coorr(1,1))-4+0.5 (height+space_v)*(coorr(1,2)) width-1.5 height];
        cur_obj2=  copyobj(subplot_4,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
      %  if list_of_cases==12
      orient landscape
 %           print('-depsc','-tiff','-r2400',[ './fig_phase_demo_assembled_' num2str(list_of_cases) '.eps']);
       % end
        
        
    end
    
    
     if list_of_cases>=2 & list_of_cases<7 ,
%         if list_of_cases==7
%             main_fig=figure(1001);clf;
%         else
%              main_fig=figure(1001);
%         end
          main_fig=figure(1001);clf;
          coorr=[0 1];
%         if list_of_cases==10,coorr=[0 0];end
%         if list_of_cases==7,coorr=[0 1];end
%         if list_of_cases==8,coorr=[0 2];end
%         if list_of_cases==9,coorr=[0 3];end
%         if list_of_cases==11,coorr=[0 4];end
%         if list_of_cases==12,coorr=[0 5];end
        width=3.2*2;
        space_h=0.8*2;
        height=3.0*2;
        space_v=0.5*2;
        test_position='Position';
        % test_position='Position';
%         coorr(1,1)=0;
%         pos=[1.9+(width+space_h)*(coorr(1,1)) (height+space_v)*(coorr(1,2)) width-2 height];
%         cur_obj1=  copyobj(subplot_0,main_fig); set(cur_obj1,'Unit','centimeters');set(cur_obj1,test_position,pos);set(main_fig, 'currentaxes', cur_obj1);drawnow;
%         
        coorr(1,1)=0;
        pos=[2+(width+space_h)*(coorr(1,1)) (height+space_v)*(coorr(1,2)) width-0.5 height];
        cur_obj2=  copyobj(subplot_1,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
         
         coorr(1,1)=1;
         pos=[2+(width+space_h)*(coorr(1,1))-0.5 (height+space_v)*(coorr(1,2)) width+0.5 height];
         cur_obj2=  copyobj(subplot_2,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
%          
%         coorr(1,1)=2;
%         pos=[2+(width+space_h)*(coorr(1,1)) (height+space_v)*(coorr(1,2)) width-2 height];
%         cur_obj2=  copyobj(subplot_3,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
%         
        coorr(1,1)=2;
        pos=[2+(width+space_h)*(coorr(1,1))+0.5 (height+space_v)*(coorr(1,2)) width-1.5 height];
        cur_obj2=  copyobj(subplot_4,main_fig); set(cur_obj2,'Unit','centimeters');set(cur_obj2,test_position,pos);set(main_fig, 'currentaxes', cur_obj2);drawnow;
      %  if list_of_cases==12
      orient landscape
   %         print('-depsc','-tiff','-r2400',[ './fig_phase_demo_assembled_' num2str(list_of_cases) '.eps']);
       % end
        
        
    end
    
end

%set(groot,'defaultLineLineWidth','remove');
