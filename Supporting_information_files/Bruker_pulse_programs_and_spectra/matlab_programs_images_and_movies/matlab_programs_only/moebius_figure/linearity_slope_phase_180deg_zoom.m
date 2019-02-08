
a=4;
b=3;
max=1+2;
figure(1);clf;
for l1=1:a
    for l2=1:b
        
        switch l1
            case 2
                beta=max*pi;
            case 3
                beta=2*max*pi;
            case 4
                beta=2.53*max*pi;
            otherwise
                beta=0;
        end
        ratio=2;
        if l2==1
            ratio=4;
        end
        % subplot(a,b,(l1-1)*b+l2)
                sp=0.02;

        wi=(1-2.3*sp)/2.2;he=(1-2*sp)/4;
        wi2=wi;
        he2=he;
        addme=0;
        if l2==3
           wi2=wi/3.5;  
           %he2=he2/2;
           addme=he/6;
                      ratio=1;

        end
        %wi=3;he=2;
        pos=[sp+(l2-1)*(wi) sp+(4-l1)*(he)+addme wi2 he2];
        subplot('Position',pos)
        
        incx=1/256;
        x=(-max/2):incx:max/2;
        zo=beta*x/max;
        
        y=cos(x*2*pi);
        z=sin(x*2*pi);
        
        
        
        iic=0.3;
        
        is=iic/2;
        if l2==1
            is=is/2;
            
        end
        for l3=-iic:is:iic
            %             zo=((l3+iic)/2*iic)*2*pi;
            if l2==3
                x=0*x;
            end
            x2=x;
            y2=y;
            z2=z;
            if l2==1
                x2=0;
                y2=x;
                % zo=beta*y/max;
                
                z2=0;
                %  zo=beta*y/max;
                
            end
            
            if l2==1
                s=1*l3*sin(1*zo);
                
                x2=x2+l3*cos(1*zo);
                
                y2=y2;
                z2=z2+s;
            else
                s=1*l3*sin(1*zo);
                
                x2=x2+l3*cos(1*zo)/ratio;
                
                y2=y2-s.*y;
                z2=z2+s.*z;
            end
            lww=2;
            main_col=[0.5 0.5 0.5];
            if iic>=(abs(l3)+0.001)
                if abs(l3)>0.01
                    plot3(x2,y2,z2,'color',main_col,'linewidth',lww);hold on
                else
                    plot3(x2,y2,z2,'color',main_col,'linewidth',lww);hold on
                    
                end
            else
                if l3>0
                    blue='b-';
                    if l2==3
                        if l1==2
                            blue='m';
                        end
                    end
                    plot3(x2,y2,z2,blue,'linewidth',lww);hold on
                    stx=x2;
                    sty=y2;
                    stz=z2;
                else
                    %                     if l2==3
                    %                         plot3(x2,y2,z2,'r-');hold on
                    %                     else
                    plot3(x2,y2,z2,'r-','linewidth',lww);hold on
                    
                    %                     end
                    enx=x2;
                    eny=y2;
                    enz=z2;
                end
            end
        end
        list=1:64:size(enx,2);
        for loperpendicular=list
            %             if loperpendicular==(size(list,2)-1)/2
            %                 col='k';al=2;
            %             else
            col=main_col;al=lww;
            %             end
            plot3([stx(1,loperpendicular) enx(end,loperpendicular)],[sty(1,loperpendicular) eny(end,loperpendicular)],[stz(1,loperpendicular) enz(end,loperpendicular)] ,'color',col,'linewidth',al);hold on
            
        end
        list=(size(enx,2)-1)/2+1;
        for loperpendicular=list
            %             if loperpendicular==(size(list,2)-1)/2
            col='k';al=lww*1.5;
            %             else
            %               col=main_col;al=1;
            %             end
            plot3([stx(1,loperpendicular) enx(end,loperpendicular)],[sty(1,loperpendicular) eny(end,loperpendicular)],[stz(1,loperpendicular) enz(end,loperpendicular)] ,'color',col,'linewidth',al);hold on
            
        end
        %          x=[x fliplr(x)+iic x(1,1)];
        %         y=[y fliplr(y) y(1,1)];
        %         z=[z fliplr(z) z(1,1)];
        %         plot3(x,y,z,'b-'); hold on
        tk=1.3;
        
        if l2==2
            axis([-max/2-tk+1 max/2+tk-1 -tk tk -tk tk])
            
        xticks([-1.5 -1 -0.5 0 0.5 1 1.5])
        yticks([])
        zticks([])
                pbaspect([ratio 1 1])
         axis('off')

        end
        if l2==3
            axis([-0.4 0.4  -tk tk -tk tk])
        xticks([ ])
        yticks([])
        zticks([])
         axis('off')
                 pbaspect([0.3 1 1])

        end
        tk2=iic;
        view([-25 31])
        
        if l2==1
            pbaspect([1 ratio 1])
            view([-40 31])
            
            axis([-tk2 tk2  -max/2 max/2 -tk2 tk2])
              xticks([])
        yticks([-1.5 -1 -0.5 0 0.5 1 1.5])
        zticks([])
        end
%         if l2==3
%             view([-37 31])
%                     view([-15 31])
% 
%         end
       % axis('off')
      
        %  print('-depsc',['./moebius_' num2str(l1) '_' num2str(l2) '.eps']);%%DEVEL
        
    end
end
print('-depsc',['./Fig_Smoebius.eps']);%%DEVEL
print('-painters','-depsc',['./Fig_Smoebius_new_rendering.eps']);%%DEVEL


%print('-fillpage','-dpdf',['./moebius.pdf']);%%DEVEL
