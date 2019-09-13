% Ellison_204977052_HW_01_main.m
% Samuel Ellison
% 204977052
% This homework has 2 problems. The first calculates and prints an exact and
% approximate surface area for an oblate spheroid with inputted radii. The
% second calculates and prints eight different approximations for the perimeter of an
% ellipse with inputted major and minor axes.

% clear cache
clc; close all; clear all;

% choose what problem to call
x = input('Enter the problem you would like to call:\n');
switch (x)
    case 1
        %% Problem 1
        % Input proper values for equatorial and polar radius
        r1 = input('Enter equatorial radius:\n');
        r2 = input('Enter polar radius:\n');
        if (r1 <= 0) || (isreal(r1) == 0)
            disp('Error: The equatorial radius must be a positive real number. Rerun script.');
        elseif (r2 <= 0) || (isreal(r2) == 0)
            disp('Error: The polar radius must be a positive real number. Rerun script.');
        elseif (r1 <= r2)
            disp('Error: The equatorial radius must be larger than the polar radius. Rerun script.');
        else
        % Define gamma, Exact Area, Approximated Area, and Percent Error
        gamma = acos(r2/r1);
        A_exact = 2*pi*(r1^2+((r2^2/sin(gamma))*log((cos(gamma))/(1-sin(gamma)))));
        A_approx = 4*pi*((r1+r2)/2)^2;
        error = (abs(A_approx - A_exact)/A_exact)*100;
        % Display Areas
        fprintf('Exact Surface Area: %10.10f\n', A_exact);
        fprintf('Approximate Surface Area: %10.10f\n', A_approx);
        fprintf('Percent Error: %10.10f', error);
        disp('%');
        end
    case 2
        %% Problem 2
        % Input proper values for major and minor axes a and b
        a = input('Enter Major Axis:\n');
        b = input('Enter Minor Axis:\n');
        if (a <= 0) || (isreal(a) == 0)
            disp('Error: The Major Axis must be a positive real number. Rerun script.');
        elseif (b <= 0) || (isreal(b) == 0)
            disp('Error: The Minor Axis must be a positive real number. Rerun script.');
        else
            % Define the departure from circle-hood, h
            h = ((a-b)/(a+b))^2;
            % Define the Perimeter Approximations for an ellipse
            P1 = pi*(a+b);
            P2 = pi*sqrt(2*(a^2+b^2));
            P3 = pi*sqrt((2*(a^2+b^2))-(0.5*(a-b)^2));
            P4 = pi*(a+b)*(1+0.125*h)^2;
            P5 = pi*(a+b)*(1+(3*h)/(10+sqrt(4-3*h)));
            P6 = pi*(a+b)*(64-3*h^2)/(64-16*h);
            P7 = pi*(a+b)*(256-48*h-21*h^2)/(256-112*h+3*h^2);
            P8 = pi*(a+b)*(3-sqrt(1-h))/2;
            % Display Perimeters and h
            fprintf('An ellipse with an h value of %1.10f', h);
            disp(' gives the following Perimeter Approximations:');
            fprintf('\n\tMethod 1: P1 = %1.10f\n', P1);
            fprintf('\tMethod 2: P2 = %1.10f\n', P2);
            fprintf('\tMethod 3: P3 = %1.10f\n', P3);
            fprintf('\tMethod 4: P4 = %1.10f\n', P4);
            fprintf('\tMethod 5: P5 = %1.10f\n', P5);
            fprintf('\tMethod 6: P6 = %1.10f\n', P6);
            fprintf('\tMethod 7: P7 = %1.10f\n', P7);
            fprintf('\tMethod 8: P8 = %1.10f\n', P8);
        end
    otherwise
        disp('Please enter a valid problem number, 1 or 2. Rerun script.');
end
