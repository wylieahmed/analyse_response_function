function [a1,a2, fr]=get_response_AOD_2beads(data_in,f,xy_slopes,xy_k,cal,rate)

%first I want to calculate the displacement. For this I need to calculate
%the position of each bead in µm anf the substract them
%first bead

for i=1:2
  data=squeeze(data_in(i,:,:));
  x(i,:)=data(1,:)/cal(5)*1e-6-1./(xy_slopes(1).*1e6).*data(4,:)./data(8,:);
  y(i,:)=data(2,:)/cal(6)*1e-6-1./(xy_slopes(2).*1e6).*data(6,:)./data(8,:);
  Fx(i,:)=xy_k(1)./(xy_slopes(1)*1e6).*data(4,:)./data(8,:);
  Fy(i,:)=xy_k(2)./(xy_slopes(2)*1e6).*data(6,:)./data(8,:);
end

%now we get the difference in x and y
dx=abs(x(1,:)-x(2,:));
dy=abs(y(1,:)-y(2,:));

p=length(dx);

alpha_x1=fft(dx)./fft(Fx(1,:));
alpha_y1=fft(dy)./fft(Fy(1,:));
alpha_x2=fft(dx)./fft(Fx(2,:));
alpha_y2=fft(dy)./fft(Fy(2,:));

fr=rate/p*([0:p/2]);

[m,p_i]=min(abs(fr-f)) 
a1=[alpha_x1(p_i),alpha_y1(p_i)];
a2=[alpha_x2(p_i),alpha_y2(p_i)];