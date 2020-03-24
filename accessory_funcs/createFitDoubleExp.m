function [fitresult, gof] = createFitDoubleExp(tToFit, yToFit)
%CREATEFIT1(TTOFIT,YTOFIT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : tToFit
%      Y Output: yToFit
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 23-Mar-2020 18:02:27


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( tToFit, yToFit );

% Set up fittype and options.
ft = fittype( 'a*exp(-b*x) + c*exp(-d*x) + e', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Robust = 'Bisquare';
opts.StartPoint = [0.171186687811562 0.706046088019609 0.757740130578333 0.743132468124916 0.421761282626275];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'yToFit vs. tToFit', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'tToFit', 'Interpreter', 'none' );
% ylabel( 'yToFit', 'Interpreter', 'none' );
% grid on

