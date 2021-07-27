% ============================================= smoothspline(x,y,mode,type)
% smoothspline.m version 6.0
%
% Created:  Brenden Epps, August  9, 2008
% Modified: Brenden Epps, April   5, 2012
%
% --------------------------------------------- Copyright 2009 Brenden Epps
% This program is free software.  You can redistribute it and/or modify it
% under the terms of the GNU General Public License version 2, as published
% by the Free Software Foundation.  This program is distributed in the hope
% that it will be useful, but WITHOUT ANY WARRANTY; without even the
% implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% -------------------------------------------------------------- Reference:
% If you find this software useful, please cite the following
% articles.  Theoretical background information is given therein:
%
%   B.P. Epps, T.T. Truscott, and A.H. Techet, "Evaluating derivatives of
%   experimental data using smoothing splines," Mathematical Methods in 
%   Engineering International Symposium, Coimbra, Portugal, October 2010.
%
%   Epps, B.P. (2010) "An impulse framework for hydrodynamic force 
%   analysis: fish propulsion, water entry of spheres, and marine 
%   propellers," Ph.D. thesis, MIT, chapter 5.
%
% See also:
%   Truscott, T.T.; Epps, B.P.; and Techet, A.H.; "Unsteady forces on 
%   spheres during free-surface water entry," Journal of Fluid Mechanics.
%
%
% ------------------------------------------------------------ Description:
% smoothspline(...) determines the best smoothing spline fit to a set of    
% presumably noisy (x,y) data. It uses the Matlab function spaps(x,y,tol) 
% to create a smoothing spline for a given error tolerance, tol, and 
% smoothspline(...) determines the "best" tol to use to capture the 
% roughness of the function that the data represents, without capturing
% the roughness due to the noise in the data.  Derivatives of this spline
% can be computed analytically, with no amplification of measurement noise. 
%
% smoothspline(...) is intended for use with high-resolution (large N), 
% high-precision (small error) data.  The user can investigate the effects
% of resolution and precision using the analytic example case given below.
%
% Inputs:
%   x     = [N,1] vector, independent variable
%   y     = [N,1] vector, dependent variable
%   mode  = 'auto'   -- automatically select best curve fit
%         = 'manual' -- manually      select best curve fit
%   type  = 2 for (default) cubic  smoothing spline (d^2f/dx^2 == 0 at endpoints) 
%         = 3 for quintic          smoothing spline (d^3f/dx^3 == 0 at endpoints) 
%
% Outputs:
%   SP    = smoothing spline data structure 
%   YY    = [N,1] vector, y(x)      spline values 
%   YYp   = [N,1] vector, dy/dx     spline values
%   YYpp  = [N,1] vector, d^2y/dx^2 spline values
%   YYppp = [N,1] vector, d^3y/dx^3 spline values 
%   ESP   = error     of output spline
%   RSP   = roughness of output spline
%   Cbest = index of selected "best" spline (1-5 in manual mode, 0 in auto mode)
% 
% ------------------------------------------------------------ Usage Notes:
%
% This function requires the MATLAB function "spaps", which is given in
% the "splines" toolbox.
%
% This function prompts the user with the following and then takes commands:
%
%         Input an action: A,H,V,D, or S
%         A ll new curve(s)
%         H ide    curve(s)
%         V iew    curve(s)
%         D ata      removed
%         S elect  curve
%
%   Description (by way of example): 
%
%   >> A
%
%   brings the user to the Roughness versus Error tolerance figure, and 
%   gives the user a cursor.  The user can then select various points on
%   this curve -- which represent different splines -- to compute.
%  
%
%   >> H125
%
%   hides splines 1, 2, and 5 (where the user can input any number spline 
%   that he/she wishes to hide).
%
%   
%   >> V12
%
%   displays splines 1 and 2.  In other workds Vxx is the opposite of Hxx, 
%   where xx can be any of the current spline(s).
%
%
%   >> D
%
%   removes the data from the Y vs X figure.
%
%
%   >> S3
%
%   selects and outputs spline 3 (again, the user can select any spline).
%   This terminates the program.
%
% ---------------------------------------------------------------- Example:
%     % Define a function to analyze
%     x  = 0:0.01:10;   
%     yt = exp(-x).*sin(x);                     % true  y(x)
%     y  = yt +  0.001*randn(size(x));          % noisy y(x)
% 
%     type  = 3;  % quintic spline fit
%
%     % Run smoothspline.m
%     [SP, YY, YYp, YYpp, YYppp] = smoothspline(x,y,'auto',type);
%     figure, plot(x,yt,'k--', x,y,'.b', x,YY,'-r')
%
%       OR
%
%     [SP, YY, YYp, YYpp, YYppp] = smoothspline(x,y,'manual',type);
%
%     % Error and roughness of output spline
%     E = trapz(x, ( y - fnval(SP,x) ).^2 )    % error     of SP
%     R = trapz(x,fnval(fnder(SP,type),x).^2)  % roughness of SP
%
%
%     % Compare output with analytic derivatives:
%     yp    =   -exp(-x).*sin(x) +   exp(-x).*cos(x);
%     ypp   = -2*exp(-x).*cos(x);
%     yppp  =  2*exp(-x).*cos(x) + 2*exp(-x).*sin(x);
%     R2    = trapz(x, ypp.^2)
%     R3    = trapz(x,yppp.^2)
% 
% =========================================================================

function [SP, YY, YYp, YYpp, YYppp, ESP, RSP, Cbest, figs] = smoothspline(x,y,mode,type)

% close all,
 
warning('off'),

if nargin == 2     % nargin is "Number of Input Arguments".  When
   mode = 'auto';  % smoothspline(x,y) is invoked, then mode = 'auto'
   type  = 2;      % use cubic smoothing spline
elseif nargin == 3
   type  = 2;      % use cubic smoothing spline
end  

% --------------------------- Perform error checking on the input arguments
if size(x) ~= size(y)
    disp('ERROR: x and y must be the same size')
    disp(' ')
    return,    
    
elseif size(x,1) > 1  && size(x,2) > 1
    disp('ERROR: x and y must be vectors')
    disp(' ')
    return,  
    
elseif size(x,1) == 1 && size(x,2) > 1
    x = x';
    y = y';
        
elseif x(end) < x(1)
    disp('ERROR: x must be strictly increasing.')
    disp(' ')
    disp('Set: x = x(end:-1:1); y = y(end:-1:1); and retry.')  
    return,
    
elseif ~(strcmp(mode,'auto') | strcmp(mode,'manual'))
    strcmp(mode,'manual')
    disp('ERROR: mode must be "auto" or "manual"')
    disp(' ')
    return,  
    
elseif type ~= 2 & type ~= 3
    disp('ERROR: type must be 2 for cubic spline or 3 for quintic spline')
    disp(' ')
    return, 
end

% ------------------------------------------------- Remove NaN and inf data
x0 = x;                    % original data sites
y  = y(find(isfinite(y))); %   finite data
x  = x(find(isfinite(y))); %   finite data sites

N = length(x); % number of data sites

%%
if strcmp(mode,'auto')   
    % ====================================================== AUTO mode code
    Cbest = 0;  % not used in this mode

    % ---------------------------- Find upper bound for E using:
    % type == 3  ---->  type-1 == 2 is a quadratic fit
    % type == 2  ---->  type-1 == 1 is a linear    fit
    Emax = trapz(x, ( y - polyval(polyfit(x,y,type-1),x) ).^2 );

    % --------- Find upper bound for R using interpolating spline (E == 0)
    Rmax = trapz(x, fnval(fnder(spaps(x,y,0,type),type),x).^2);

    % ----- Estimate Ecr
    % Ecr   ~ N e^2 dt
    % R2max = 36 N e^2 dt^-3 ~ 36 Ecr dx^-4   (cubic   splines)
    % R3max = 31 N e^2 dt^-5 ~ 31 Ecr dx^-6   (quintic splines)
    dx  = mean(diff(x));

    if type == 2     % (cubic   splines)
        Ecr = Rmax / (36 * dx^-4);
    elseif type == 3 % (quintic splines)
        Ecr = Rmax / (31 * dx^-6);
    end

    % --------------------- Define initial double-bisection stencil, E(1:5)
    disp('Defining initial double-bisection stencil...')
    E  = zeros(5,1); % error tolerance
    R  = zeros(5,1); % roughness

    E(1) = 1e-14;            % note: machine zero is 2.2204e-16
    E(2) = sqrt(Ecr*1e-14);  % log(E2) = 1/2 (log(Ecr)+log(1e-14))
    E(3) = Ecr;
    E(4) = sqrt(Ecr*Emax);
    E(5) = Emax;

    % ---------------------------------- Find roughness at initial stencil, 
    %                                    and augment stencil if necessary
    % Find roughness of point (1)
    disp('Finding roughness at point (1)...')
    R(1) = trapz(x,fnval(fnder(spaps(x,y,E(1),type),type),x).^2); 

    if R(1) == 0
        disp('ABORT: Roughness is zero with error tolerance 1e-14')
        return
    end

    % Make sure that roughness is nonzero at E(3)
    Rflag = 0;

    while Rflag == 0
        disp('Finding point (3) with non-zero roughness...')

        R(3)  = trapz(x,fnval(fnder(spaps(x,y,E(3),type),type),x).^2);  

        if R(3) == 0
            E(3) = sqrt(E(3)*1e-14);  % make log(E(3)) half as much
        else
            Rflag = 1;
        end
    end

    % Find roughness of point (2)
    disp('Finding roughness at point (2)...')
    E(2) = sqrt(E(3)*1e-14);  % log(E(2)) = 1/2 (log(E(3))+log(1e-14))
    R(2) = trapz(x,fnval(fnder(spaps(x,y,E(2),type),type),x).^2);  


    % Find roughness of point (5), and set to near-zero if necessary
    disp('Finding roughness at point (5)...')
    R(5) = trapz(x,fnval(fnder(spaps(x,y,E(5),type),type),x).^2);  

    if log10(R(5)) == -Inf
        R(5) = 10^-100;  
    end

    % Make sure that roughness is nonzero at point (4)
    E(4) = sqrt(E(3)*E(5));  % log(E(4)) = 1/2 (log(E(3))+log(E(5)))

    Rflag = 0;

    while Rflag == 0
        disp('Finding point (4) with non-zero roughness...')

        R(4)  = trapz(x,fnval(fnder(spaps(x,y,E(4),type),type),x).^2);  

        if R(4) == 0
            E(4) = sqrt(E(3)*E(4));
        else
            Rflag = 1;
        end
    end

    % --- Ensure that the initial stencil has a point of positive curvature
    Rpp  = [NaN 0 0 0 NaN]';  % R-E curvature:  d^2(log(R))/d(log(E))^2
    
    for i = 2:4
       Rpp(i) = ((log(R(i+1))-log(R(i)))/(log(E(i+1))-log(E(i))) - (log(R(i))-log(R(i-1)))/(log(E(i))-log(E(i-1))))/(0.5*(log(E(i+1))-log(E(i-1))));    
    end

    % If no point with positive curvature, adjust the E(3) and try again
    if isempty(find(Rpp>0))
        disp('Adjusting point (3)...'),
        E(3) = sqrt(E(3)*E(4));
        R(3) = trapz(x,fnval(fnder(spaps(x,y,E(3),type),type),x).^2); 
    
        for i = 2:4
            Rpp(i) = ((log(R(i+1))-log(R(i)))/(log(E(i+1))-log(E(i))) - (log(R(i))-log(R(i-1)))/(log(E(i))-log(E(i-1))))/(0.5*(log(E(i+1))-log(E(i-1))));
        end
    end
    
    % If no point with positive curvature, adjust the E(2) and E(4) and try again
    if isempty(find(Rpp>0))
        disp('Adjusting points (2) and (4)...'),
        E(2) = sqrt(E(3)*E(2));
        R(2) = trapz(x,fnval(fnder(spaps(x,y,E(2),type),type),x).^2); 
        
        E(4) = sqrt(E(3)*E(4));
        R(4) = trapz(x,fnval(fnder(spaps(x,y,E(4),type),type),x).^2); 
        
        for i = 2:4
            Rpp(i) = ((log(R(i+1))-log(R(i)))/(log(E(i+1))-log(E(i))) - (log(R(i))-log(R(i-1)))/(log(E(i))-log(E(i-1))))/(0.5*(log(E(i+1))-log(E(i-1))));
        end
    end
    
    % If no point with positive curvature, then suggest manual mode
    if isempty(find(Rpp>0))
        disp('Warning: there may be no critical point.  Try manual mode.')
        return,
    end    
    % ---------------------------------------------------------------------

    % --------------- Perform double-bisection iterations until convergence
    done  = 0;
    count = 0;

    while done == 0 & count < 100
        count = count + 1;
        disp(['Performing double-bisection iteration ',num2str(count)]),

        % Find curvature of interior stencil points
        for i = 2:4
           Rpp(i) = ((log(R(i+1))-log(R(i)))/(log(E(i+1))-log(E(i))) - (log(R(i))-log(R(i-1)))/(log(E(i))-log(E(i-1))))/(0.5*(log(E(i+1))-log(E(i-1))));    
        end
        
        % Find stencil point with max curvature
        p = find(Rpp == max(Rpp));

            %         figure(1),
            %             plot(E,R,'.r')
            %             plot(E(p),R(p),'.g')
            %             axis([1e-14 Emax 1 Rmax])
            %             pause,

        % Redraw stencil and evaluate roughness of two new splines
        Eold = E;
        E(1) = Eold(p-1);
        E(2) = sqrt(Eold(p-1)*Eold(p));
        E(3) = Eold(p);
        E(4) = sqrt(Eold(p)*Eold(p+1));
        E(5) = Eold(p+1);

        SP2  = spaps(x,y,E(2),type);        
        SP4  = spaps(x,y,E(4),type);        
        
        Rold = R;
        R(1) = Rold(p-1);
        R(2) = trapz(x,fnval(fnder(SP2,type),x).^2);
        R(3) = Rold(p);
        R(4) = trapz(x,fnval(fnder(SP4,type),x).^2);
        R(5) = Rold(p+1);


        % Check for convergence (if 1% difference between stencil points)
        if abs((log(E(4)) - log(E(3)))/log(E(4))) < 0.01 
            done = 1;
        end
    end  % double-bisection while loop
    
    % ------------------------------------ Output selected smoothing spline
       SP = SP4;
      ESP = trapz(x, ( y - fnval(SP4,x) ).^2 );    % error     of SP
      RSP = trapz(x,fnval(fnder(SP4,type),x).^2);  % roughness of SP       
       YY = fnval(fnxtr(      SP   ),x0);   
      YYp = fnval(fnxtr(fnder(SP,1)),x0);  
     YYpp = fnval(fnxtr(fnder(SP,2)),x0);  
    YYppp = fnval(fnxtr(fnder(SP,3)),x0);     

    % ================================================== END AUTO mode code     
else
    % ==================================================== MANUAL mode code
%% 

    % ---------------------------- Find upper bound for E using:
    % type == 3  ---->  type-1 == 2 is a quadratic fit
    % type == 2  ---->  type-1 == 1 is a linear    fit
    Emax = trapz(x, ( y - polyval(polyfit(x,y,type-1),x) ).^2 );

    % ---------- Find upper bound for R using interpolating spline (E == 0)
    Rmax = trapz(x,fnval(fnder(spaps(x,y,0,type),type),x).^2);  

    % ------------------------- Generate R-E frontier for smoothing splines 
    disp(' '),
    disp('Generating (E,R) relation for smoothing splines ...'),

    E = 10.^[floor(10*log10(Emax))/10:-0.1:-14];    % error tolerance

    R = zeros(size(E));                             % roughness

    for i = 1:length(E)
        disp(['Evaluating spline ',num2str(i),' of ',num2str(length(E))]),

        sp    = spaps(x,y,E(i),type);                 % spline fit

        R(i)  = trapz(x,fnval(fnder(sp,type),x).^2);  

        
        % Stop computing splines if R(i) is nearly Rmax
        if abs((Rmax - R(i))/Rmax) < 0.1
            R(i+1:end) = NaN;

            break,
        end
    end

    % -------------------------- Find (E,R) point of max positive curvature
    Rpp = 0*R;

    for i = 2:length(E)-1
       Rpp(i) = ((log(R(i+1))-log(R(i)))/(log(E(i+1))-log(E(i))) - (log(R(i))-log(R(i-1)))/(log(E(i))-log(E(i-1))))/(0.5*(log(E(i+1))-log(E(i-1)))); 
    end

    Ibest = find(Rpp == max(Rpp));  % index of point with max positive curvature

    C     = 5;                     % number of choices of spline fit
        
    if Ibest < C
        Ibest = C;  
    elseif length(Ibest) > 1
        Ibest = find(R == max(R));
    end
    
     
    % -------------------------- Generate C choices of possible spline fits
    disp(' '),
    disp('Generating possible smoothing spline choices...'),
     
    DC       = ones(1,C);             % Display Curves    flag
    DD       = 1;                     % Display Data      flag

    TOL      = E(Ibest:-1:Ibest-C+1); % E of each spline fit choice
    Rsfc     = R(Ibest:-1:Ibest-C+1); % R of each spline fit choice
%     Rsfc     = zeros(C,1);            % R of each spline fit choice
    yy       = zeros(C,N);            % y(x)
    yyp      = zeros(C,N);            % dy/dx
    yypp     = zeros(C,N);            % d^2y/dx^2
    yyppp    = zeros(C,N);            % d^3y/dx^3

    for i = 1:C
           sp(i)   = spaps(x,y,TOL(i),type);         

           yy(i,:) = fnval(      sp(i)   ,x);   
          yyp(i,:) = fnval(fnder(sp(i),1),x);  
         yypp(i,:) = fnval(fnder(sp(i),2),x);  
        yyppp(i,:) = fnval(fnder(sp(i),3),x);  

%         Rsfc(i)    = trapz(x,fnval(fnder(sp(i),type),x).^2);  
    end
    % ---------------------------------------------------------------------

    % -------------------------------------------------------- Plot results
    disp(' '),
    disp('Plotting results...'),
    
    set(0,'DefaultAxesFontSize',10,'DefaultAxesFontName','Times')
    
    % ------------------------------------ Allocate memory for plot handles
    HHyy    = zeros(1,C);
    HHyyp   = zeros(1,C);
    HHyypp  = zeros(1,C);
    HHyyppp = zeros(1,C);
    HHerr   = zeros(1,C);
    HHint   = zeros(1,C);

    % --------------------------------------------------------- Color matrix
    CLR = [     1       0       0;      ... % (1) Red
                0       0.9     0;      ... % (2) Green
                0       0       1;      ... % (3) Blue
                0.75    0       0.75;   ... % (4) Purple
                1       0.5     0;      ... % (5) Orange
                0       1       1;      ... % (6) Cyan
                1       0       1;      ... % (7) Magenta
                0.75    0.5     0.25;   ... % (8) Brown
                0.25    0.25    0.75];      % (9) Navy blue

    CLRlegend_all = {'(1)','(2)','(3)','(4)','(5)','(6)','(7)','(8)','(9)'};           

    CLRlegend = CLRlegend_all(1:C);
    DATlabel = {'data'};  
%%   
    % --------------------------------------------------------------- (X,Y)
    Hy = figure;
        set(Hy,'Position',[130 680 560 420])
        set(Hy,'Position',[-25 450 560 420]) 
        hold off,
        plot(x,0*x,'k--'), hold on,
        HHdata = plot(x,y,'ko');
        for i = 1:C
            HHyy(i) = plot(x,yy(i,:),'color',CLR(i,:),'linewidth',1);    
        end
        set(gca,'XLim',[x(1) x(end)])
        xlabel('X','Fontsize',12,'Fontname','Times')
        ylabel('Y','Fontsize',12,'Fontname','Times')
        legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')
    % ---------------------------------------------------------------------     

    % -------------------------------------------------------------- (X,Yp)
    Hyp = figure;
        set(Hyp,'Position',[130 185 560 420])
        set(Hyp,'Position',[-25 0 560 420])
        hold off,
        plot(x,0*x,'k--'), hold on,
        for i = 1:C    
            HHyyp(i) = plot(x,yyp(i,:),'color',CLR(i,:),'linewidth',1);
        end
        set(gca,'XLim',[x(1) x(end)])
        xlabel('X' ,'Fontsize',12,'Fontname','Times')
        ylabel('dY / dX','Fontsize',12,'Fontname','Times')
        legend(HHyyp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
    % --------------------------------------------------------------------- 

    % ------------------------------------------------------------- (X,Ypp)
    Hypp = figure;
        set(Hypp,'Position',[695 185 560 420])
        set(Hypp,'Position',[460 0 560 420])
        hold off,
        plot(x,0*x,'k--'), hold on,
        for i = 1:C
            HHyypp(i) = plot(x,yypp(i,:),'color',CLR(i,:),'linewidth',1);
        end
        set(gca,'XLim',[x(1) x(end)])
        xlabel('X'  ,'Fontsize',12,'Fontname','Times')
        ylabel('d^2Y / dX^2','Fontsize',12,'Fontname','Times')
        legend(HHyypp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')  
    % ---------------------------------------------------------------------   

    % ------------------------------------------------------------ (X,Yppp)
    Hyppp = figure;
        set(Hyppp,'Position',[1260 185 560 420])
        set(Hyppp,'Position',[945 0 560 420])
        hold off,
        plot(x,0*x,'k--'), hold on,
        for i = 1:C
            HHyyppp(i) = plot(x,yyppp(i,:),'color',CLR(i,:),'linewidth',1);
        end
        set(gca,'XLim',[x(1) x(end)])
        xlabel('X'   ,'Fontsize',12,'Fontname','Times')
        ylabel('d^3Y / dX^3','Fontsize',12,'Fontname','Times')
        legend(HHyyppp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
    % -------------------------------------------------------------------------   

   % --------------------------------------------------------- (tol,int(Yppp))
    Hint = figure;
        pp = find(isnan(R),1);
        
        if length(pp) == 0
            pp = length(R);
        end
        
        set(Hint,'Position',[1260 680 560 420])
        set(Hint,'Position',[945 450 560 420])
        hold on, 
        HHRvsTOL = loglog(E(1:pp),R(1:pp),'k.');
        loglog([E(1),E(pp)],[Rmax,Rmax],'k--');
        for i = 1:C
            HHint(i) = loglog(TOL(i),Rsfc(i),'.','color',CLR(i,:));
        end
  
        xlabel('Error tolerance, E','Fontsize',12,'Fontname','Times')
        ylabel('Roughness, R'      ,'Fontsize',12,'Fontname','Times')
        set(gca,'XScale','log','YScale','log')
        if any(DC == 1)
            legend(HHint(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
        end
        
    % Store figures in array to return.
    figs = [Hy, Hyp, Hypp, Hyppp, Hint];
    
    % -------------------------------------------------------------------------
    % -------------------------------------------------------------------------

    
    % ------------------------------- Manually select best smoothing spline
    done = 0;

    while done == 0
        
        disp(' ')
        disp('Input an action: A,H,V,D, or S')
        disp('A ll new curve(s)')
        disp('H ide    curve(s)')
        disp('V iew    curve(s)') 
        if DD == 1    % display estimates flag is currently on
            disp('D ata      removed')
        else          % display estimates flag is currently off
            disp('D ata      displayed')
        end
        disp('S elect  curve   ')
        disp(' ')

        pause(0.0001),
        action = input('  >> ','s');
        
        
        if any(['H','h'] == action(1))
            if length(action) == 1
                RMcurve = input('Hide curve(s): ','s');
            else
                RMcurve = action(2:end);
            end

            for i = 1:length(RMcurve)
                DC(str2num(RMcurve(i))) = 0;
            end

            % ----------------------- Update which curve fits are displayed
            figure(Hy),
                set(HHyy          ,'visible','off')
                set(HHyy(find(DC)),'visible','on')
                legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')

            figure(Hyp),
            	set(HHyyp          ,'visible','off')  
            	set(HHyyp(find(DC)),'visible','on')    
            	legend(HHyyp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')

            figure(Hypp),
            	set(HHyypp          ,'visible','off')  
            	set(HHyypp(find(DC)),'visible','on')    
            	legend(HHyypp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')

            figure(Hyppp),
            	set(HHyyppp          ,'visible','off')  
            	set(HHyyppp(find(DC)),'visible','on')    
            	legend(HHyyppp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
                
            figure(Hint),
            % -------------------------------------------------------------

        elseif any(['V','v'] == action(1))
            if length(action) == 1
                ADDcurve = input('View curve(s): ','s');
            else
                ADDcurve = action(2:end);
            end    

            for i = 1:length(ADDcurve)
                DC(str2num(ADDcurve(i))) = 1;
            end    

            % ----------------------- Update which curve fits are displayed
            figure(Hy),
                set(HHyy          ,'visible','off')
                set(HHyy(find(DC)),'visible','on')
                legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')

            figure(Hyp),
            	set(HHyyp          ,'visible','off')  
            	set(HHyyp(find(DC)),'visible','on')    
            	legend(HHyyp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')

            figure(Hypp),
            	set(HHyypp          ,'visible','off')  
            	set(HHyypp(find(DC)),'visible','on')    
            	legend(HHyypp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')

            figure(Hyppp),
            	set(HHyyppp          ,'visible','off')  
            	set(HHyyppp(find(DC)),'visible','on')    
            	legend(HHyyppp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')  

            figure(Hint),
            % -------------------------------------------------------------

        elseif any(['D','d'] == action(1))
            if DD == 1
                DD = 0;
                
                % ------- Update whether original data is displayed
                figure(Hy),
                        set(HHdata  ,'visible','off')  
                    legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')
                
                figure(Hyp),
                figure(Hypp),
                figure(Hyppp),
                figure(Hint),    

            else
                DD = 1;
                
                % ------- Update whether original data is displayed
                figure(Hy),
                        set(HHdata  ,'visible','on')  
                    legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')
                
                figure(Hyp),
                figure(Hypp),
                figure(Hyppp),
                figure(Hint),
            end
            % -------------------------------------------------------------            

        elseif any(['S','s'] == action(1))
            if length(action) == 1
                Cbest = input('Select curve: ');
            else
                Cbest = str2num(action(2));
            end    

            DC        = zeros(1,C);
            DC(Cbest) = 1;

            done = 1;

%             % ----------------------- Update which curve fits are displayed
%             figure(Hy),
%                 set(HHyy          ,'visible','off')
%                 set(HHyy(find(DC)),'visible','on')
%                 legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')
% 
%             figure(Hyp),
%             	set(HHyyp          ,'visible','off')  
%             	set(HHyyp(find(DC)),'visible','on')    
%             	legend(HHyyp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
% 
%             figure(Hypp),
%             	set(HHyypp          ,'visible','off')  
%             	set(HHyypp(find(DC)),'visible','on')    
%             	legend(HHyypp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
% 
%             figure(Hyppp),
%             	set(HHyyppp          ,'visible','off')  
%             	set(HHyyppp(find(DC)),'visible','on')    
%             	legend(HHyyppp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')  
% 
%             figure(Hint),
%             % -------------------------------------------------------------    
 
        elseif any(['A','a'] == action(1))
            disp(' ')
            disp('Click on (E,R) points to generate curve fits.')
            disp('Press "return" when done.')
            
            % -------------------- All new selection of (Ippp,tol) point(s)
            [tolDES,junk] = getpts(Hint);

            % ----------------------- Round tolDES to the nearest tol value
            for i = 1:length(tolDES)
                [min_difference, array_position] = min(abs(E - tolDES(i)));

                tolDES(i) = E(array_position);
            end

            % --------------------------------------- Update curves to plot
            TOL = unique(tolDES');
            C   = length(TOL); 
            
            if C > 9
                disp(' ')
                disp('Note: max display of 9 curves.')
                C   = 9;
                TOL = TOL(1:9);
            end
            
            HHyy     = zeros(1,C);
            HHyyp    = zeros(1,C);
            HHyypp   = zeros(1,C);
            HHyyppp  = zeros(1,C);
            HHerr    = zeros(1,C);
            HHint    = zeros(1,C);
            yy       = zeros(C,N);            % y(x)
            yyp      = zeros(C,N);            % dy/dx
            yypp     = zeros(C,N);            % d^2y/dx^2
            yyppp    = zeros(C,N);            % d^3y/dx^3
            Rsfc     = zeros(C,1);            % roughness of spline fit choices 
            
            
            for i = 1:C
                disp(['Computing new spline fit choice ',num2str(i),' of ',num2str(C)])
            
                   sp(i)   = spaps(x,y,TOL(i),type);        
                      
                   yy(i,:) = fnval(      sp(i)   ,x);   
                  yyp(i,:) = fnval(fnder(sp(i),1),x);  
                 yypp(i,:) = fnval(fnder(sp(i),2),x);  
                yyppp(i,:) = fnval(fnder(sp(i),3),x);  

                Rsfc(i)    = trapz(x,fnval(fnder(sp(i),type),x).^2);  
            end            
            
            % ----------------------------------------- Plot all new curves
            CLRlegend = CLRlegend_all(1:C); % Legend
            DC        = ones(1,C);          % Display Curves flag

            figure(Hy)
                hold off,
                plot(x,0*x,'k--'), hold on,
                HHdata    = plot(x,y,'ko');
                for i = 1:C
                    HHyy(i) = plot(x,yy(i,:),'color',CLR(i,:),'linewidth',1);    
                end
                xlabel('X','Fontsize',12,'Fontname','Times')
                ylabel('Y','Fontsize',12,'Fontname','Times')
                legend([HHdata(find(DD)),HHyy(find(DC))],[DATlabel(find(DD)),CLRlegend(find(DC))],'Fontsize',12,'Fontname','Times')
                
            figure(Hyp)
                hold off,
                plot(x,0*x,'k--'), hold on,
                for i = 1:C    
                    HHyyp(i) = plot(x,yyp(i,:),'color',CLR(i,:),'linewidth',1);
                end
                xlabel('X'      ,'Fontsize',12,'Fontname','Times')
                ylabel('dY / dX','Fontsize',12,'Fontname','Times')
            	legend(HHyyp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
                
            figure(Hypp)
                hold off,
                plot(x,0*x,'k--'), hold on,
                for i = 1:C
                    HHyypp(i) = plot(x,yypp(i,:),'color',CLR(i,:),'linewidth',1);
                end
                xlabel('X'          ,'Fontsize',12,'Fontname','Times')
                ylabel('d^2Y / dX^2','Fontsize',12,'Fontname','Times')
            	legend(HHyypp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
                
            figure(Hyppp)    
                hold off,
                plot(x,0*x,'k--'), hold on,
                for i = 1:C
                    HHyyppp(i) = plot(x,yyppp(i,:),'color',CLR(i,:),'linewidth',1);
                end
                xlabel('X'          ,'Fontsize',12,'Fontname','Times')
                ylabel('d^3Y / dX^3','Fontsize',12,'Fontname','Times')
            	legend(HHyyppp(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')

            figure(Hint)
                hold off,              
                HHRvsTOL = loglog(E(1:pp),R(1:pp),'k.'); hold on, 
                loglog([E(1),E(pp)],[Rmax,Rmax],'k--');
                for i = 1:C
                    HHint(i) = loglog(TOL(i),Rsfc(i),'.','color',CLR(i,:));
                end
                xlabel('Error tolerance, E','Fontsize',12,'Fontname','Times')
                ylabel('Roughness, R'      ,'Fontsize',12,'Fontname','Times')
                set(gca,'XScale','log','YScale','log')
                legend(HHint(find(DC)),CLRlegend(find(DC)),'Fontsize',12,'Fontname','Times')
            % -------------------------------------------------------------
       end
    end

    % ---------------------------------------- Output selected smoothing spline
       SP =   sp(Cbest);
      ESP =  TOL(Cbest);
      RSP = Rsfc(Cbest);
       YY = fnval(fnxtr(      SP   ),x0);   
      YYp = fnval(fnxtr(fnder(SP,1)),x0);  
     YYpp = fnval(fnxtr(fnder(SP,2)),x0);  
    YYppp = fnval(fnxtr(fnder(SP,3)),x0); 
    % ------------------------------------------------ END manual mode code
%%    
end
% ============================================================= END program
% =========================================================================

