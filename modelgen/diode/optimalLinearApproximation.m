% ===============================================================
% == Global Linear Approximation by Dynamic Programming        ==
% ==                                                           ==
% ==   Supplement of the work:                                 ==
% ==                                                           ==
% ==     From Cardiac Cells to Genetic Regulatory Networks     ==
% ==                                                           ==
% ==   Description:                                            ==
% ==     An adapted version (it takes several curves in input) ==
% ==     of the algorithm by Prez and Vidal, Pattern           ==
% ==     Recognition Letters 15, 743-750.                      ==
% ==                                                           ==
% ==   Input:                                                  ==
% ==         x,y:   Curves given as an x-points vector and     ==
% ==         a y-points vectors                                ==
% ==         S:     number >= 2 of desired segments            ==
% ==         p:     0 = do not plot,   1 = plot                ==
% ==                                                           ==
% ==   Output:                                                 ==
% ==         er:    error vector                               ==
% ==         a,b:   segment coefficients                       ==
% ==         xb:    x-coordinate at breaking point             ==
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
function [er,a,b,xb] = optimalLinearApproximation(x,y,S,p)
  %% Initialize constants and matrices
  
  P = size(x,2);       % number of points of the digitized curves
  nCurves = size(y,1); % number of curves
  se = zeros(1,4);
  % Cost tables
  cost = ones(P,S)*inf;  % cost(30,4)= min cost to point 30 with 4 segments
  hCst = ones(P,P)*inf;  % hCst(i, n)= cached error of segment i to n  
  cost(2,1) = 0;         % zero segment cost from point 1 to point 2
  
  % Predecessor table
  father = ones(P,S)*inf; % father(30,4) no of starting pt of 4th segment

  %% Compute optimal segmentation
  % Initialize 1 segment cost and father from point 1 to all other points
  for n=2:P      
       for i=1:nCurves
           se(i) = segmentErr(x(1:n),y(i,1:n),0); 
       end 
       cost(n,1) = max(se); %segmentError(x(1:n),y(i,1:n),r,0) 
       father(n,1) = 1;
  end;
  

  % Compute 2 <= n <= S segments cost from point 1 to all other points
  for m=2:S                                         % segments loop
     for n=3:P                                      % points loop
          minErr=cost(n-1,m-1); minIndex = n-1;
          for i=m:n-2                               % intermediate points
              if (hCst(i,n)==Inf)
                    for k=1:nCurves
                        se(k) = segmentErr(x(i:n),y(k,i:n),0); 
                    end
                  hCst(i,n) = max(se); %segmentError(x(i:n),y(i:n),r,0);
              end;
              currErr = cost(i,m-1) + hCst(i,n);
              if (currErr < minErr)
                  minErr = currErr; minIndex = i;
              end;
          end;
          cost(n,m) = minErr;
          father(n,m) = minIndex;
      end;
  end;
  %totalErr = cost(n,m);

  %% Extract answer and optionally plot it

  % Intialize answer matrices
  ib = zeros(S,S+1); xb = zeros(S,S+1); yb = zeros(nCurves, S,S+1);
   a = zeros(nCurves, S,S);    b = zeros(nCurves, S,S);   er = zeros(nCurves, S,S);
 
  % Extract fit for all segmentations up to S
  for (j=S:-1:1)
    ib(j,j+1) = P; 
    xb(j,j+1) = x(ib(j,j+1)); 
    for k=1:nCurves
        yb(k,j,j+1) = y(k,ib(j,j+1));
    end
    
    if (p) 
        figure; 
        for k=1:nCurves
            hold on;
            plot(x,y(k,:),xb(j,j+1),yb(k,j,j+1),'or'); 
        end
    end; 
    
    % Compute other points and optionally plot
    for (i=j:-1:1)
        ib(j,i) = father(ib(j,i+1),i);
        xb(j,i) = x(ib(j,i)); 
        for k=1:nCurves
            yb(k,j,i) = y(k,ib(j,i));
        end
        if (p) 
            for k=1:nCurves
                hold on, 
                plot(x(ib(j,i)),y(k,ib(j,i)),'or'); 
            end
        end;
        for k=1:nCurves
            [er(k,j,i),a(k,j,i),b(k,j,i)] = segmentErr(x(ib(j,i):ib(j,i+1)),y(k,ib(j,i):ib(j,i+1)),p);
        end    
    end;
  end;
      
end % optLinAppDynProg
     
