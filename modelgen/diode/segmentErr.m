% ===============================================================
% == Compute error and optionally plot the fit                 ==
% ==                                                           ==
% ==   Supplement of the work:                                 ==
% ==                                                           ==
% ==     From Cardiac Cells to Genetic Regulatory Networks     ==
% ==                                                           ==
% ==   Description:                                            ==
% ==     Compute the linear interpolation for the input        ==
% ==     curve x,y. Output the coefficients error.             ==                                          ==
% ==                                                           ==
% ==   Input:                                                  ==
% ==         x,y:   Curves given as an x-points vector and     ==
% ==         a y-points vectors                                ==
% ==         p:     0 = do not plot,   1 = plot                ==
% ==                                                           ==
% ==   Output:                                                 ==
% ==         e:     error                                      ==
% ==         a,b:   segment coefficients                       ==
% ==                                                           ==
% ==   Authors:                                                ==
% ==                                                           ==
% ==     R. Grosu, E. Bartocci                                 ==
% ==                                                           == 
% ==   Date:  11/05/10                                         ==
% ==                                                           ==
% ==   Free distribution with authors permission               ==
% ==                                                           ==
% ==   SUNY Stony Brook, Stony Brook, NY                       ==
% ==                                                           == 
% ===============================================================  
function [e,a,b] = segmentErr(x,y,p)
  % Find out size of x,y  
  n = size(x,2);
  % Compute line segment coefficients   
  a = (y(n)-y(1))/(x(n)-x(1));
  b = (y(1)*x(n)-y(n)*x(1))/(x(n)-x(1));
  % Compute error for above line segment
  e=0; 
  for (k=1:n) e = e + (y(k)-a*x(k)-b)^2/(a^2+1); end; %LS perpend
  % Optionally plot segment
  if (p) hold on, plot(x,a*x+b); end;
end % error  
