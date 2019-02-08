clear all
%for main_ratio=[ 0.01 1 0.5 0.2 0.1],
for main_ratio=[ 1],
    for mainlooop=0:2,
        figure(mainlooop+1)
        Z = peaks;
        surf(Z)
        axis tight manual
        ax = gca;
        ax.NextPlot = 'replaceChildren';
        
        
        writerObj = VideoWriter(['mov' num2str(mainlooop) '.avi']);
        open(writerObj);
        
        axis tight
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        
        
        loops = 360;
     %   F(loops) = struct('cdata',[],'colormap',[]);
        counter=1;
        for loopj = 1:3:89,
            % X = sin(loopj*pi/10)*Z;
            %surf(X,Z)
            clf
            
            anglet= loopj-1;
            phase_a=cos(anglet/360*2*pi)+j*sin(anglet/360*2*pi);%as complex number
            phase_a=phase_a/(abs(phase_a));
            [dist_in_hz erro_in_deg]=shap_fn3d(loopj-1,mainlooop,main_ratio);
            
            drawnow;
            F(loopj) = getframe;
            writeVideo(writerObj,F(loopj));
            
        end
        
        close(writerObj);
        %movie(F,2)
    end
%     figure(111)
%     clf
%     plot(store_dis_in_hz(:,1),store_erro_in_deg(:,1),'b-');
%     hold on
%     plot(store_dis_in_hz(:,2),store_erro_in_deg(:,2),'r-');
%     print('-depsc','-tiff','-r600',[ 'Phase_error_nearby_small_signals' num2str(main_ratio)  '.eps']);%here
    
    
end
