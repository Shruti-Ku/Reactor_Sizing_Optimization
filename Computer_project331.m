% Given Input Points
% Inlet conc.
Ca0 = [2, 5, 6, 6, 11, 14, 16, 24]; 
% Outlet conc.
Ca = [0.5, 3, 1, 2, 6, 10, 8, 4];
% Tau
t = [30, 1, 50, 8, 4, 20, 20, 4]; 
%inverse of rate(1/(-ra))
l_rate = t./(Ca0-Ca);  

%fitting curve through Input Points points
curve = interp1(Ca, l_rate,'spline','pp');

% range of x values
Ca_range = linspace(min(Ca), max(Ca), 100); 

% inverse rate function values
l_rate_func = ppval(curve, Ca_range); 

% Plotting given Input Points and the inverse rate function curve
plot(Ca, l_rate, '*', 'DisplayName', 'Given Input Points');
hold on;

% plot fitted curve
Ca_range = linspace(min(Ca), max(Ca), 100);
plot(Ca_range, l_rate_func, 'DisplayName', 'Fitted inverse rate function curve');

% set plot labels, title, and legend
xlabel('C_a (mmol/m^{3})');
ylabel('-1/r_a');
title('Inverse Rate Function Curve');
legend('show');


%given conversion, conc. and volume flow rate
Xa=0.9;
Cin=10;
Cout=Cin*(1-Xa);
vo=0.1;



%Providing options to user
Method = menu("Select a method",'A single PFR', 'A single CSTR', 'Two stirred tanks of any size', 'Combination of a PFR and a MFR', 'A PFR with recycle');

if Method==1
    %For PFR

    fitted_curve = @(x) ppval(curve, x);
    Ca_pfr=linspace(Cout,Cin,100); 
    l_ra_pfr = ppval(curve,Ca_pfr);
    % integral function to compute the area under the curve
    Area = integral(fitted_curve, Cout, Cin);

    Volume_of_pfr=vo*Area;
    fprintf('Minimum required volume of PFR: %.5f m^3\n', Volume_of_pfr);
    % Plotting given Input Points and the inverse rate function curve
    figure(2)
    plot(Ca, l_rate,'*');
    hold on
    plot(Ca_range, l_rate_func);
    xlabel('C_a (mmol/m^{3})');
    ylabel('-1/r_a');
    hold on
    area(Ca_pfr , l_ra_pfr, 'FaceColor', '#F4C2C2', 'EdgeColor', 'r','linestyle',':','linewidth',2);
        legend('Data', 'Fitted curve','PFR Area');
        title('PFR');

end 
if Method==2
    
   % For CSTR
   % range of x values
    Ca_pfr = linspace(Cout,Cin, 100);
    % fitted y values
    l_ra_pfr = ppval(curve, Ca_pfr); 
    Area = (Cin-Cout)*ppval(curve,Cout);
    Volume_of_cstr=vo*Area;
    fprintf('Minimum required volume of CSTR: %.5f m^3\n', Volume_of_cstr);
     ploty=ppval(curve,Cout)+zeros(1,100);
    % Plotting given Input Points and the inverse rate function curve
    figure(3)

area(Ca_pfr, ploty, 'FaceColor', '#F4C2C2', 'EdgeColor', 'r','linestyle',':','linewidth',2);
  hold on
  plot(Ca, l_rate,'*');
hold on
plot(Ca_range, l_rate_func);
xlabel('C_a (mmol/m^{3})');
ylabel('-1/r_a');
    legend('CSTR Area','Input Points', 'Curve');
    title('CSTR');
end


if Method==3
    % Derivative of the fitted inverse rate function at the point x
    deriv = @(x) ppval(fnder(curve), x);

    % Find the root of the equation (by the method of maximization of
    % rectangles)
    Maximization_of_rect= @(x) ((10 - ppval(curve, x)) / (x - Cin)) - deriv(x);
    xnot = fzero(Maximization_of_rect, [min(Ca), max(Ca)-0.001]);
    %For 1st CSTR 
    Area_CSTR1=(Cin-xnot)*ppval(curve, xnot);
    V_CSTR1=Area_CSTR1*vo;
    %For 2nd CSTR
    Area_CSTR2=(xnot-Cout)*ppval(curve, Cout);
    V_CSTR2=Area_CSTR2*vo;
    fprintf('Minimum required volume of CSTR1: %.4f m^3\n', V_CSTR1);
    fprintf('Minimum required volume of CSTR2: %.4f m^3\n', V_CSTR2);

    figure(4);
    xlabel('C_a (mmol/m^{3})');
    ylabel('-1/r_a')
    
    hold on;
    xval1 = linspace(Cin,xnot, 100);
   plot1=ppval(curve,xnot)+zeros(1,100);
    xval2 = linspace(xnot,Cout, 100);
    plot2=ppval(curve,Cout)+zeros(1,100);
    hold on;
    area(xval1, plot1, 'FaceColor', '#F4C2C2','linestyle',':', 'linewidth',2,'EdgeColor', 'r');
    area(xval2, plot2, 'FaceColor', '#3EB489','linestyle',':','linewidth',2, 'EdgeColor', 'r');
     plot(Ca, l_rate,'*');
     hold on;
plot(Ca_range, l_rate_func);
    legend('CSTR1', 'CSTR2','Input Points', 'Curve','Location', 'northeast');
    title('Two CSTR in Series');
    

end
if Method ==4
    
    Ca_range_PFR_MFR = linspace(min(Ca), max(Ca), 100);

% Determining point of minima of the fitted inverse rate function
    

% Search for the minimum value of the spline function within the range of x-values
    [x0, Y0] = fminbnd(@(x) ppval(curve, x), min(Ca_range_PFR_MFR), max(Ca_range_PFR_MFR));
    
    fitted_curve = @(x) ppval(curve, x);

    % integral function to compute the area under the curve
    Area = integral(fitted_curve, Cout, x0);
    Volume_PFR=Area*vo;
    Ca_range_PFR=linspace(Cout,x0,100);
    l_ra_PFR=ppval(curve,Ca_range_PFR);
    Area_CSTR =(Cin-x0)*ppval(curve,x0);
    Volume_CSTR=Area_CSTR*vo;
    
    fprintf('Minimum required volume of PFR: %.4f m^3\n', Volume_PFR);
    fprintf('Minimum required volume of CSTR: %.4f m^3\n', Volume_CSTR);

    figure(5);
    
    plot1=ppval(curve,x0)+zeros(1,100);
    xval1 = linspace(Cin,x0, 100);
    area(xval1, plot1, 'FaceColor', 'm','linestyle',':', 'linewidth',2, 'EdgeColor', 'r');
    hold on
     area(Ca_range_PFR, l_ra_PFR, 'FaceColor', '#F4C2C2','linestyle',':', 'linewidth',2,'EdgeColor', 'r');
     hold on
  plot(Ca, l_rate,'*');
hold on
plot(Ca_range, l_rate_func);
xlabel('C_a (mmol/m^{3})');
ylabel('-1/r_a');
   
    legend('CSTR', 'PFR','Input Points', 'Curve');
    title('PFR and MFR in Series');
   
end 

if Method==5
   Ca_recycle=linspace(1, 10, 1000);
    l_ra_recycle=ppval(curve, Ca_recycle);
    [x0, Y0] = fminbnd(@(x) ppval(curve, x), min(Ca_recycle), max(Ca_recycle));
    for i = 1:length(Ca_recycle)
        if (Ca_recycle(i) < x0)
          
            Cr1=fzero(@(x) ppval(curve, x) - ppval(curve, Ca_recycle(i)),10);
            if Cr1>10
                continue;
            end 
            h=ppval(curve,Cr1);
            Cspan=linspace(Cout,Ca_recycle(i));
            y = ppval(curve, Cspan) - h;
            area_under = trapz(Cspan, y);
            Cspan2 = linspace(Ca_recycle(i),Cr1);
            y = h - ppval(curve, Cspan2);
            area_above=trapz(Cspan2, y);
            if abs(area_above-area_under)<=0.1
            Cr2=Ca_recycle(i);
            break
        end 
        

    else 
      
        Cr1=fzero( @(x) ppval(curve, x) - ppval(curve, Ca_recycle(i)),0);
        if Cr1<1
            continue;
        end 
        h=ppval(curve,Cr1);
        Cspan = linspace(Cout,Cr1); % range of x values
        % Compute the area between the two curves
        y = ppval(curve, Cspan) - h; % function to integrate
        area_above = trapz(Cspan, y); % integrate over x_range
        Cspan2 = linspace(Cr1,Ca_recycle(i)); % range of x values
      
        % Compute the area between the two curves
        y = h- ppval(curve, Cspan2); % function to integrate
        area_under = trapz(Cspan2, y); % integrate over x_range
        if abs(area_above-area_under)<=0.1
            Cr2=Ca_recycle(i);
            break
        end 
    end   
    end
    area_recycle=h*(Cin-Cout);
    R_ratio=(Cin-Cr1)/(Cr1-Cout);
    Volume_PFR_Recycle = area_recycle*vo;
    fprintf('PFR volume with recycle: %.5f m^3\n', Volume_PFR_Recycle);
     figure(6)
     xr=linspace(Cin,Cout, 100);
     ploty=ppval(curve,Cr2)-1+ones(1,100);
      area(xr, ploty, 'FaceColor', '#F4C2C2','linestyle',':','linewidth',2, 'EdgeColor', 'r');
      hold on;
      plot(Ca, l_rate,'*');
 hold on;
 plot(Ca_range, l_rate_func);
 xlabel('C_a (mmol/m^{3})');
 ylabel('-1/r_a');
      
     legend('PFR with Recycle','Input Points', 'Curve');
     title('PFR with Recycle');
end 