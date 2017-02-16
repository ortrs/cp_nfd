% Test file for test-benching
% All units are normalized to [um]

x_size=50;      %[um] size of simulation area in x
y_size=50;      %[um] size of simulation area in y
sx_size=512;    %[px] number of pixels in x direction 
sy_size=512;    %[px] number of pixels in y direction
z_size=20;      %[n/a] number of z-direction iterations
lambda=1;       %[um] excitation field wavelength
end_vector=10;  %[n/a] wavelength multiplier for end value simulation
                % e.g. end_vector=10 -> 10*lambda = 10 um

[v_in,x,y]=slit(x_size,y_size,4,8,sx_size,sy_size);
%Uncomment this and comment slit for double slit
%[v_in,x,y]=doubleslit(x_size,y_size,4,4,0,200,sx_size,sy_size);

%ZX-plane vector definition
%z=logspace(0.1*lambda,10*lambda,z_size);
zx = linspace(0.01*lambda,end_vector*lambda,z_size);
xx = linspace(-x_size,x_size,sx_size);
[Z,X] = meshgrid(xx,zx);

%Part one!
%Discrete calculation of fields, change z for specific fields

% % for z = [0*lambda,0.5*lambda,1*lambda,5*lambda,50*lambda,5000*lambda]
% %     [v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,z,lambda,'biDFT');
% %     figure(2);
% %     subplot(1,2,1)
% %     mesh(x,y,abs(v_out))
% %     view(0,90)
% %     axis tight
% %     colorbar
% %     title(sprintf('Field Amplitude xy z=%.1f\\lambda',z));
% %     subplot(1,2,2)
% %     plot(x(1,:),abs(v_out(sx_size/2,:)))
% %     axis tight
% %     title(sprintf('Field Amplitude x-cut z=%.1f\\lambda',z));
% %     drawnow
%GIF-storing
    %figure(2);
% % % %     filename = 'diffall_10mm_all.gif';
% % % %     frame = getframe(1);
% % % %     im = frame2im(frame);
% % % %    [imind,cm] = rgb2ind(im,256);
% % % %    if z == 1
% % % %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
% % % %    else
% % % %          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
% % % %    end
% % % %end

%Part two!
%Continuous field propagation along z, XY-plane definitions given above

for z = linspace(0.1*lambda,10*lambda,z_size)
%for z = 0
    [v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,z,lambda,'biDFT');
    figure(1);
    subplot(2,2,1)
    mesh(x,y,abs(v_out))
    view(0,90)
    axis tight
    colorbar
    title(sprintf('Field Amplitude xy z=%.3f\\lambda',z));
    subplot(2,2,2)
    plot(x(1,:),abs(v_out(sx_size/2,:)))
    axis tight
    title(sprintf('Field Amplitude x-cut z=%.3f\\lambda',z));
    subplot(2,2,3)
    mesh(Kx,Ky,abs(V_out))
    view(0,90)
    axis tight
    colorbar
    title(sprintf('Transform Amplitude xy z=%.3f\\lambda',z));
    subplot(2,2,4)
    mesh(Kx,Ky,angle(V_out))
    view(0,90)
    axis tight
    title(sprintf('Transform Phase z=%.3f\\lambda',z));
    drawnow
    
    figure(2);
    subplot(1,2,1)
    mesh(x,y,real(P_f))
    view(0,90)
    axis tight
    colorbar
    title(sprintf('Real propagator exponent xy z=%.3f\\lambda',z));
    subplot(1,2,2)
    mesh(x,y,imag(P_f))
    view(0,90)
    colorbar
    axis tight
    title(sprintf('Imag propagator exponent x-cut z=%.3f\\lambda',z));
    drawnow
%GIF-storing
% % %    figure(1);
% % %    filename = 'diff0p1l_10l_presCP_p5.gif';
% % %    frame = getframe(1);
% % %    im = frame2im(frame);
% % %    [imind,cm] = rgb2ind(im,256);
% % %    if z == 0.1*lambda
% % %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
% % %    else
% % %          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
% % %    end
% % %    figure(2);
% % %    filename = 'diff0p1l_10l_presCP_p6.gif';
% % %    frame = getframe(2);
% % %    im = frame2im(frame);
% % %    [imind,cm] = rgb2ind(im,256);
% % %    if z == 0.1*lambda
% % %         imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1);
% % %    else
% % %          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1);
% % %    end
% % %  end

%Part three!
%Continuous field propagation along z, ZX-plane, definitions given above

% % for z = linspace(0.01*lambda,end_vector*lambda,z_size)
% %     [v_out,V_out,V_in,Kx,Ky,P_f,p_f]=propagation(v_in,x,y,z,lambda,'biDFT');
% %     if z==0.01*lambda
% %         vx_out = [v_out(sx_size/2,:)];
% %     else
% %         vx_out = [vx_out;v_out(sx_size/2,:)];
% %     end
% % end
% %     figure(3);
% %     mesh(X,Z,abs(vx_out));
% %     view(0,90)
% %     axis tight
% %     title(sprintf('Field amplitude along z-axis',z));
% %     xlabel('propagation z-coordinates [\mum]') % y-axis label
% %     ylabel('x [\mum]') % x-axis label