clear all
%for main_ratio=[ 0.01 1 0.5 0.2 0.1],
for loop_type=[ 0 1]%1
for mainlooop=-2
    % for mainlooop=-3
    list_factor=mainlooop+0.01;
    list_an=[90 ];
    if (mainlooop==-1)
        list_factor=[0 1 2 3 4 5 6];
        
    end
    if mainlooop==-2 %not working does not know why
        list_factor=[-0.50:0.1:0.5];
       % list_factor=[-0.60:0.2:0.6];
        list_factor=[-1:0.2:1];
        list_factor=[-0.5:0.1:0.5];
        list_an=[90 180 270 360];
        list_an=[360 270 180 90];
    end
    if mainlooop==-3 %not working does not know why
        list_factor=[-6:0.2:6];
        list_an=[90 ];
    end
    for main_ratio=list_an
        
        figure(mainlooop+4)
        Z = peaks;
        surf(Z)
        axis tight manual
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        if loop_type
            extra_text='sphere';
        else
            extra_text='time';
            
        end
        
        writerObj = VideoWriter(['mov' num2str(mainlooop) '_' num2str(main_ratio) '_' extra_text '.avi']);
        open(writerObj);
        
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        
        list_of_values= 0.00:1:main_ratio+1;
        loops = 360;
        loops = 1;
        
        F(size(list_of_values,2)+1) = struct('cdata',[],'colormap',[]);
        counter=1;
        for loopj = list_of_values
            % X = sin(loopj*pi/10)*Z;
            %surf(X,Z)
            clf
            
            anglet= loopj-1;
            % phase_a=cos(anglet/360*2*pi)+j*sin(anglet/360*2*pi);%as complex number
            %  phase_a=phase_a/(abs(phase_a));
            %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            %  [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            
         %   if main_ratio==270
          %      fig_gen_spheres(list_factor,loopj,[152 19]);
         %   else
         if  loop_type
                fig_gen_spheres(list_factor,loopj)
                
         else
                             fig_gen_time(list_factor,loopj);
                             
         end
         %  end
         set(gcf,'color','w');
         drawnow;
         tmp_frame = getframe;
         if counter==1
             si=size(tmp_frame.cdata);
%              if loop_type
%                  si=si-2;
%              end
         end
         if loop_type
             
             tmp_frame.cdata=tmp_frame.cdata(1:si(1,1),1:si(1,2),1:si(1,3));
         else
            
             tmp_frame.cdata=tmp_frame.cdata(1:si(1,1),1:si(1,2),1:si(1,3));
             
         end
         
         
         
         F(counter)=tmp_frame;
         writeVideo(writerObj,F(counter));
            counter=counter+1;
            
            
        end
        
        close(writerObj);
        %movie(F,2)
                 if loop_type

        print('-depsc','-tiff','-r600',[ 'Final_time_domain_sphere_' num2str(main_ratio) '.eps']);%here
                 else
                         print('-depsc','-tiff','-r600',[ 'Final_time_domain_time_' num2str(main_ratio) '.eps']);%here
    
        
                 end
    end
    %     figure(111)
    %     clf
    %     plot(store_dis_in_hz(:,1),store_erro_in_deg(:,1),'b-');
    %     hold on
    %     plot(store_dis_in_hz(:,2),store_erro_in_deg(:,2),'r-');
    %     print('-depsc','-tiff','-r600',[ 'Phase_error_nearby_small_signals' num2str(main_ratio)  '.eps']);%here
    
end   
end
