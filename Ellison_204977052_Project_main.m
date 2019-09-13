% Ellison_204977052_Project_main.m
% Samuel Ellison
% 204977052
% FINAL PROJECT: This script simulates the Spinodal Decomposition over time of an
% Initial Distribution composed of two pure phases. The script solves the
% Cahn-Hillard Equation using different Runge-Kutta methods and Laplacian
% stencil approximations. 4 different simulations are run for different RK,
% stencil, and time-step combinations. Plots and videos of the simulations
% are produced. A function that approximates the Laplacian is used in this
% script. 

% clear cache
clc; close all; clear all; rng('shuffle');

% Constants
a = 1;
b = 1;
gamma = 1;
D = 3;
h = 1;
tFinal = 20;

% Choose Simulation/Video --> Defines stencil, method, and time-step
number = input('Please Enter Simulation:\n\n(1) FIVE POINT STENCILS\n(2) NINE POINT STENCILS\n(3) SECOND ORDER RUNGE KUTTA\n(4) ANIMATION OF TIME EVOLUTION\n(5) EXTRA CREDIT\n');
switch(number)
    case 1
        stencil = 5;
        method = 1;
        dt = 0.001;
    case 2
        stencil = 9;
        method = 1;
        dt = 0.0001;
    case 3
        stencil = 5;
        method = 2;
        dt = 0.005;
    case 4       
        vid = input('Enter video:\n\n(1) FIVE POINT STENCIL || RK = 1\n(2) NINE POINT STENCIL || RK = 1\n(3) FIVE POINT STENCIL || RK = 2\n(4) NINE POINT STENCIL || RK = 4\n');
        
        switch(vid)
            case 1
                stencil = 5;
                method = 1;
                dt = 0.001;
            case 2
                stencil = 9;
                method = 1;
                dt = 0.0001;
            case 3
                stencil = 5;
                method = 2;
                dt = 0.005;
            case 4
                stencil = 9;
                method = 4;
                dt = 0.0001;
            otherwise
                error('Invalid Video Number. Try again.')
        end
        
    case 5
        stencil = 9;
        method = 4;
        dt = 0.0001;
    otherwise
        error('Invalid Simulation Number. Try again.')
end

% Phi matrix dimensions
rows = 150;
cols = 100;

% Fill Phi matrix with 1 or -1
% for k = 1:1:rows
%     for j = 1:1:cols
%         r = rand;
%         if r <= 0.5
%             phi(k,j) = 1;
%         else
%             phi(k,j) = -1;
%         end
%     end
% end

% Load the previoulsy created initial distribution
load mydata.mat phi

% Normal Plots without video
if number ~= 4
        
    % Index and position for plot
    j = 1;
    position = {'northwest', 'northeast', 'southwest', 'southeast', 'center'};

    % Iterate
    for k = 0:dt:tFinal

        % Plot every 5 seconds
        if k == 0 || k == 5 || k == 10 || k == 15 || k == 20

            figure(j)
            movegui(position{j})
            imagesc(phi)
            drawnow
            axis([0 cols 0 rows])
            xticks(0:20:cols)
            yticks(0:30:rows)
            if k == 0
                title(sprintf('Initial Field Distribution [RK %i, %i pt Stencil]', method, stencil))
            else
                title(sprintf('Field Distribution [RK %i, %i pt Stencil] (t = %is)', method, stencil, k))
            end
            colorbar
            caxis([-1 1])

            % Update figure index
            j = j + 1;

        end
        
        % Update the phi matrix depending on method
        % Use Cahn-Hillard
        if (method == 1)

            % RK 1 %
            y = phi;
            C = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f = (D)*(Laplacian_2D(C, h, stencil));
            phi = phi + dt*f;

        elseif (method == 2)

            % RK 2 %
            y = phi;
            C1 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f1 = (D)*(Laplacian_2D(C1, h, stencil));
            c1 = dt*f1;

            y = phi + 0.5*c1;
            C2 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f2 = (D)*(Laplacian_2D(C2, h, stencil));
            c2 = dt*f2;

            phi = phi + c2;

        elseif (method == 4)

            % RK 4 %
            y = phi;
            C1 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f1 = (D)*(Laplacian_2D(C1, h, stencil));
            c1 = dt*f1;

            y = phi + 0.5*c1;
            C2 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f2 = (D)*(Laplacian_2D(C2, h, stencil));
            c2 = dt*f2;

            y = phi + 0.5*c2;
            C3 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f3 = (D)*(Laplacian_2D(C3, h, stencil));
            c3 = dt*f3;

            y = phi + c3;
            C4 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f4 = (D)*(Laplacian_2D(C4, h, stencil));
            c4 = dt*f4;

            phi = phi + (c1/6) + (c2/3) + (c3/3) + (c4/6);

        end

    end

    
% Videos
else
    
    % Create video with approproiate name
    myVideo = VideoWriter(sprintf('%iPS_RK%i.mp4', stencil, method),'MPEG-4');
    
    % set frame rate and open
    myVideo.FrameRate = 30;
    open(myVideo)
    
    % Set indices to track
    count = 0;
    t_out = 0;
    
    
    % Iterate through 10s
    for k = 0:dt:10
            
        % Add a frame to the video
        if k >= t_out

            % "plot" figure
            figure('visible', 'off')
            imagesc(phi)
            axis([0 cols 0 rows])
            xticks(0:20:cols)
            yticks(0:30:rows)
            title (sprintf('Field Distribution [RK %i, %i pt Stencil]', method, stencil))
            colorbar
            caxis([-1 1])

            % Add it as a frame to video
            Fr(count+1) = getframe(gcf);           
            writeVideo(myVideo,Fr(count+1));

            % update counters
            count = count + 1;
            t_out = t_out + (1/30);
            

        end
                     
        % Update the phi matrix depending on method
        % Use Cahn-Hillard
        if method == 1
            
            % RK 1 %
            y = phi;
            C = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f = (D)*(Laplacian_2D(C, h, stencil));
            phi = phi + dt*f;
            
            
        elseif method == 2
            
            % RK 2 %
            y = phi;
            C1 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f1 = (D)*(Laplacian_2D(C1, h, stencil));
            c1 = dt*f1;

            y = phi + 0.5*c1;
            C2 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f2 = (D)*(Laplacian_2D(C2, h, stencil));
            c2 = dt*f2;

            phi = phi + c2;
            
            
        elseif method == 4
            
            % RK 4 %
            y = phi;
            C1 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f1 = (D)*(Laplacian_2D(C1, h, stencil));
            c1 = dt*f1;

            y = phi + 0.5*c1;
            C2 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f2 = (D)*(Laplacian_2D(C2, h, stencil));
            c2 = dt*f2;

            y = phi + 0.5*c2;
            C3 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f3 = (D)*(Laplacian_2D(C3, h, stencil));
            c3 = dt*f3;

            y = phi + c3;
            C4 = (b^4)*(y.^3) - (a)*(b^2)*(y) - (gamma)*(Laplacian_2D(y, h, stencil));
            f4 = (D)*(Laplacian_2D(C4, h, stencil));
            c4 = dt*f4;

            phi = phi + (c1/6) + (c2/3) + (c3/3) + (c4/6);
            
            
        end
            
    end
    
    % close the video
    close(myVideo)
    
end



