%Test file for test-benching

x_size=400;
y_size=400;
z_steps=15;
sx_size=512;
sy_size=512;
lambda=1;
z=200;
fig=6;

[v_in,x,y]=slit(x_size,y_size,10,20,sx_size,sy_size);
%[v_in,x,y]=doubleslit(x_size,y_size,4,8,0,20,sx_size,sy_size);
vx_out = zeros (1,sx_size);
vx_out2 = zeros (1,sx_size);
zx = linspace(1,z*2,2*z_steps);
xx = linspace(-x_size,x_size,sx_size);
[Z,X] = meshgrid(xx,zx);
figure(fig);

%[v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,1000,lambda,'lens');
%v_in = v_out;
%mesh(abs(v_out));
%view(0,90);
%axis tight;
%title('HARAMBLErGH');
%drawnow;
%v_out=v_in;

for z = linspace(0,z,z_steps)
%for z = [7000]
    [v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,z,lambda,'biDFT');
    if z==0
        vx_out = [v_out(sx_size/2,:)];
    else
        vx_out = [vx_out;v_out(sx_size/2,:)];
    end
    
    subplot(1,2,1)
    mesh(abs(v_out));
    view(0,90);
    axis tight;
    title(sprintf('z=%d',z));  
    subplot(1,2,2)
    mesh(angle(v_out));
    view(0,90);
    axis tight;
    title(sprintf('z=%d',z));
    drawnow;    
    %v_in=v_out;
end

%Change field at the end of the propagation to new input
v_in = v_out;

 [v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,1000,lambda,'lens');
 figure(fig+1);
 vx_out = [vx_out;v_out(sx_size/2,:)];
 v_in = v_out;
 subplot(1,2,1)
 mesh(abs(v_out));
 view(0,90);
 axis tight;
 title('field after lens');
 subplot(1,2,2)
 mesh(angle(v_out));
 view(0,90);
 axis tight;
 title('field after lens');
 drawnow;

%propagate again at z'

for zp = linspace(0,z,z_steps)
%for z = [7000]
    [v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,zp,lambda,'biDFT');
   vx_out = [vx_out;v_out(sx_size/2,:)];
    subplot(1,2,1)
     mesh(abs(v_out));
     view(0,90);
     axis tight;
     title('field after lens');
     subplot(1,2,2)
     mesh(angle(v_out));
     view(0,90);
     axis tight;
    title(sprintf('zp=%d',zp));
    drawnow;   
end

 figure(fig+2);
    mesh(X,Z,abs(vx_out(2:end,:)))
    view(0,90)
    xlabel('propagation z-coordinates (in um)') % y-axis label
    ylabel('|E| [a.u.]') % x-axis label
    axis tight
    %colormap jet
    drawnow