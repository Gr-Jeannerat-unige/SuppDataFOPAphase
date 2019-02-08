function set_pos(ma_font_size)
%%set(groot,'defaultLineLineWidth',3);
this_ax=gca;
pos_ax=2*[0 0 4 3];set(this_ax,'Unit','centimeters');set(this_ax,'OuterPosition',pos_ax);


set(this_ax,'FontSize',ma_font_size);
% 
%  set(0,'defaultlinelinewidth',0.5);
%  set(0,'defaultaxeslinewidth',0.5);
%  set(0,'defaultpatchlinewidth',0.5);
%set(groot,'defaultLineLineWidth','remove')
end