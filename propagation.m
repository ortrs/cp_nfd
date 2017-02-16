function [v_out,V_out,V_in,K_x,K_y,Prop,prop]=propagation(v_in,x,y,z,lambda,sw)    
% Propagation over a predefined distance    
% Solution to the Helmholtz equation    
% All spatial coordinates in the unit of the wavelength    
% function call [v_out,V_out,V_in,Kx,Ky]=propagation(v_in,x,y,z);    
% v_in: input field (2D-matrix)    
% x,y: Spatial coordinates (2D-matrices)    
% z: propagation distance    
% v_out: output field in the spatial domain (2D-matrix)    
% V_out: output field in the Fourier domain (2D-matrix)    
% V_in:  input field in Fourier domain (2D-matrix)    
% Kx,Ky: Coordinates of the Fourier space (2D-matrix)
% sw: Use one or another method available for computation:
    % "manDFT" = Manual DFT. Slow and probably would crash
    % "redDFT" = reduced DFT. Using vectorization
    % "biDFT" = built-in DFT. 
    
    %wave-vector defintion
    k = 2*pi/lambda;
    %Pixel size dimensions
    [Nx,Ny] = size(v_in);
    %delta_x,y definition
    deltax = max(max(x))/(Nx-1); %-1 factor for odd-numbered pixels
    deltay = max(max(y))/(Ny-1);
    const = (1./(Nx*Ny)); % Normalization constant
    %Initial definition for field size
    Kx = linspace(Nx/4/min(min(x)),Ny/4/max(max(x)),Nx);
    Ky = linspace(Ny/4/min(min(y)),Ny/4/max(max(y)),Ny);
    %Preferred definition of field size
    K_x= linspace(-pi/deltax,pi/deltax,Nx);
    K_y= linspace(-pi/deltay,pi/deltay,Nx);
    % Mesh grid definitions
    [Kx,Ky] = meshgrid(Kx,Ky);
    [K_x,K_y] = meshgrid(K_x,K_y);

    if (strcmp(sw,'lens'))
        
        n_lens=1.2; %refractive index of glass
        R_1=+50; %radius of lens  in(in um)
        R_2=-50; %radius of lens out(in um)
        delta_0=0;
        V_in =  fftshift(fft2(v_in));
        %focal point calculation
        f = 1/((n_lens-1).*(1/R_1 - 1/R_2))       
        deltaxy = -R_1.*(1 -sqrt(1-((x.^2 +y.^2)./(R_1^2))) ) ...
                  +R_2.*(1 -sqrt(1-((x.^2 +y.^2)./(R_2^2))));
        prop =  exp(-1i*k*n_lens*delta_0) .* ...
                exp(-1i*k*(n_lens-1) .* deltaxy);
        %propagator definition
        Prop =  fftshift(fft2(prop));
        %v_out = v_in .*prop;
        v_out = V_in .*prop;
        Prop=fftshift(fft2(prop));
        V_out=V_in.*prop;
        v_out=ifft2(ifftshift(v_out));   
    end
    
    if (strcmp(sw,'biDFT'))
       %FFT of input field
       V_in =  fftshift(fft2(v_in));
       %Prop = exp ((1i).*sqrt(k^2 - Kx.^2 -Ky.^2).* z);
       %Propagator definition in reciprocal space
       Prop = (1i.*sqrt(k^2 - K_x.^2 -K_y.^2).*z);
       %Prop = sqrt(k^2 - Kx.^2 -Ky.^2).* z;
       %Euclidean (i.e. Actual) space propagator (visualization only)
       prop = ifft2(fftshift(Prop));
       %Reciprocal space output field 
       V_out = V_in .*  exp (1i.*sqrt(k^2 - K_x.^2 -K_y.^2).*z);
       %V_out = V_in .*  exp (1i.*sqrt(k^2 - Kx.^2 -Ky.^2).*z);
       %V_out = fftshift(V_out);
       
       %Euclidean space output field 
       v_out = ifft2(ifftshift(V_out));
    end
       
    if (strcmp(sw,'redDFT'))
        v_in = reshape (v_in,1,Nx*Ny);
        V_in = zeros (1,Nx.*Ny);
        m=0:Nx*Ny-1;
        m=repmat(m,Nx*Ny,1);
        m=m.*m';
        m = exp(-1i.*(2*pi/(Nx*Ny).*(m)));
        mi = exp(+1i.*(2*pi/(Nx*Ny).*(m)));
        V_out = (m*v_in')';
        V_out = reshape(V_out,Nx,Ny);
        %fft-shift re-shaping
        V_out = [circshift(V_out(1:Ny/2,1:end),[0 Nx/2-1]);
                 circshift(V_out(Ny/2+1:end,1:end),[0 Nx/2])];
        V_out = circshift(V_out,[Nx/2 0]);
        %V_out =  [V_out(Nx/2+1:end,Ny/2+1:end), ...
        %          V_out(Nx/2+1:end,1:Ny/2);
        %          V_out(1:Nx/2,Ny/2+1:end), ...
        %          V_out(1:Nx/2,1:Ny/2)];  
        %wave propagation
        V_out = V_out .* ...
             exp (1i .*sqrt(k^2 - Kx.^2 -Ky.^2).* z);
        trans = exp (1i .*sqrt(k^2 - Kx.^2 -Ky.^2).* z);
        %fft-shift
        %V_cout = circshift(V_out,[Ny/2 Nx/2]);
        V_cout = [circshift(V_out(1:Ny/2,1:end),[0 Nx/2-1]);
                 circshift(V_out(Ny/2+1:end,1:end),[0 Nx/2])];
        V_cout = circshift(V_out,[Nx/2 0]);
        %V_cout =  [V_out(Nx/2+1:end,Ny/2+1:end), ...
        %          V_out(Nx/2+1:end,1:Ny/2);
        %          V_out(1:Nx/2,Ny/2+1:end), ...
        %          V_out(1:Nx/2,1:Ny/2)];
        v_out=reshape (V_cout,1,Nx*Ny);
        v_out = (1/(Nx*Ny)).*(m*v_out')';
        v_out = reshape(v_out,Nx,Ny);
        %mesh(x,y,real(v_out));
        %figure;
        %mesh(x,y,imag(v_out));
    end
 
    if (strcmp(sw,'manDFT'))

        for m=1:Nx*Ny
            for n=1:Nx*Ny
                V_in(m) =  v_in(n) .* exp(-1i .* (2*pi*const).*m.*n); 
            end
        end
    end
    %V_in = 1/const .* reshape(Nx,Ny);
    %tac = 1;
end

