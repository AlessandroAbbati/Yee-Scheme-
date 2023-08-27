%==========================================================================
% EN.625.719 Numerical Diff. EQ.
% Alessandro Abbati
% Implement Yee Numerical Scheme
%==========================================================================
clear;
close all;
clc;


alpha = 8;
c = 3*(10^8);
Z = 376.7;
dtau = alpha/16; 
dx = alpha/8; 
dy = alpha/8; 

t = 0;
x = linspace(((50*alpha) - (c*t)),((58*alpha) - (c*t)),2*alpha);
EzSource = sin(((x' - (50*alpha) + (c*t))*pi)/(8*alpha));
EzSource = repmat(EzSource,1,8*alpha);

Nx = 81;
Ny = 98; 
nSteps = 95; 

putObstacle = 1; 
showPlots = 1; 

waveStart = 50;
waveStop = waveStart + (2*alpha-1); 

Ez(:,:,1) = zeros(Nx,Ny);
Ez(waveStart:waveStop,19:82,1) = EzSource;

Hy(:,:,1) = (1/Z)*Ez(:,:,1);
Hx(:,:,1) = zeros(Nx,Ny);

M=moviein(nSteps);
for t = 1:nSteps

    %boundary conditions
    Hx(:, 98, t) = 0; 
    Hx(:, 1, t) = 0; 
    Hx(81, :, t) = 0; 
    Hx(1, :, t) = 0; 
    Ez(81, :, t) = 0; 
    Ez(1, :, t) = 0; 

    %obstacle boundary conditions
    if putObstacle == 1
        %top and bottom boundary conditions of the square
        Ez(17:49, 65, t) = 0;
        Hy(17:49, 65, t) = 0;
        Ez(17:49, 33, t) = 0;
        Hy(17:49, 33, t) = 0;

        %left and right boundary conditions of the square
        Ez(49, 33:65, t) = 0;
        Hx(49, 33:65, t) = 0;
        Ez(17, 33:65, t) = 0;
        Hx(17, 33:65, t) = 0;
    end


    % Update the Magnetic Field
    for i = 2:Nx-1
        for j = 2:Ny-1
            Hy(i,j,t+1) = Hy(i,j,t) + ((1/Z)*(dtau/dx)*(Ez(i+1,j,t)-Ez(i,j,t)));
            Hx(i,j,t+1) = Hx(i,j,t) - ((1/Z)*(dtau/dy)*(Ez(i,j+1,t)-Ez(i,j,t)));
        end
    end 


    % Update the Electric Field
    for i = 2:Nx-1
        for j = 2:Ny-1
            Ez(i,j,t+1) = Ez(i,j,t) + ((Z*dtau/dx)*(Hy(i,j,t+1)-Hy(i-1,j,t+1))) - ((Z*dtau/dy)*(Hx(i,j,t+1)-Hx(i,j-1,t+1)));
        end 
    end 


    surf(transpose(Ez(:, :, t+1))); 
    title_str = strcat("Propagating Half Sine Wave, t = ", num2str(t+1));
    title(title_str);
    xlabel("i")
    ylabel("j")
    M(:,t) = getframe ;


end 
movie(M,1);


if showPlots == 1

    t = 5;
    j = 30; 
    figure(2)
    y = transpose(Ez(:, j, t+1));
    plot(y)
    title_str = strcat("Ez by means of (14a)-(14c). (j, t) = (", num2str(j), ", ", num2str(t), " )");
    title(title_str)
    xlabel("i")
    ylabel("Amplitude")

    t = 35;
    j = 30; 
    figure(3)
    y = transpose(Ez(:, j, t+1));
    plot(y)
    title_str = strcat("Ez by means of (14a)-(14c). (j, t) = (", num2str(j), ", ", num2str(t), " )");
    title(title_str)
    xlabel("i")
    ylabel("Amplitude")


    t = 65;
    j = 30;     
    figure(4)
    y = transpose(Ez(:, j, t+1));
    plot(y)
    title_str = strcat("Ez by means of (14a)-(14c). (j, t) = (", num2str(j), ", ", num2str(t), " )");
    title(title_str)
    xlabel("i")
    ylabel("Amplitude")


    t = 95;
    j = 30;     
    figure(5)
    y = transpose(Ez(:, j, t+1));
    plot(y)
    title_str = strcat("Ez by means of (14a)-(14c). (j, t) = (", num2str(j), ", ", num2str(t), " )");
    title(title_str)
    xlabel("i")
    ylabel("Amplitude")   

    % t = 5; surf(transpose(Ez(:, :, t))); 
    % title_str = strcat("Propagating Half Sine Wave, t = ", num2str(t));
    % title(title_str);
    % xlabel("i")
    % ylabel("j")
    
end




