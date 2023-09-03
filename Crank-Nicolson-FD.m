classdef CN
    %   The class contains analytical BS solution of an European and American option,
    %   Crank Nicolson finite difference method and calculation of relative
    %   error
    
    %   The Crank Nicolson function implements the Thomas algorithm
    %   (Thomas.m)
    
    properties (Constant)  
    end
    
    methods (Static)
        
        % Compute analytical price of an European option
        function [BSPrice, BSPricevec] = BSmodel(S0,S, T, K, r, sigma,xSteps, Type)
            
            %Preallocation
            BSPricevec = zeros(xSteps+1,1);
            
            % Black Scholes model
            std = sqrt(T) * sigma;
            d1 = (log(S/K) + r*T + T*sigma*sigma*0.5) / std; 
            d2 = d1 - std;
            
        switch Type
            case 'CALL'
                BSPricevec = exp(-r*T)*K * normcdf(-d2) - S.*normcdf(-d1);
                BSPrice = interp1(S,BSPricevec,S0);
                
            case 'PUT'
               BSPricevec = normcdf(-d2)*exp(-r*T)*K - S.*normcdf(-d1);
               BSPrice = interp1(S,BSPricevec,S0);

        end
        end
        
        % Calculate relative error 
        function relativeError = relativeError(CNPrice, BSPrice)
            relativeError = abs(BSPrice-CNPrice)/BSPrice;
        end
            
        function BSPrice = truesoln(S0, tgrid, K, r, sigma, Type)
            
        std = sqrt(tgrid) * sigma;
        d1 = (log(S0/K) + r.*tgrid + sigma*sigma*0.5.*tgrid) / std; 
        d2 = d1 - std;
        
        switch Type
            case 'CALL'
            BSPrice = exp(-r*T)*K * normcdf(-d2) - S0.*normcdf(-d1);
                
            case 'PUT'
            BSPrice = normcdf(-d2).*exp(-r.*tgrid)*K - S0.*normcdf(-d1);

        end
        end
        
        % CN pricing function of European and American options
        function [Price, Pricevec] = CN_Method(K,S0,r,T,sigma,xmax,xmin,S,t,theta,Type)
        
        % Get the number of grid points
        xSteps = length(S)-1;
        tSteps = length(t)-1;

        % Get the grid sizes (assuming equi-spaced points)
        dtau = T/tSteps;
        dx = (xmax-xmin)/xSteps;

        % Preallocation of w matrix to collect prices
        w(xSteps+1,tSteps+1) = nan;

        % Specify the boundary conditions
        switch Type
            case 'EURCALL'
                % Time boundary condition
                w(:,end) = max(S-K,0);
                w(end,:) = CN.truesoln(xmax, t(end:-1:1), K, r, sigma, 'CALL');
                w(1,:) = CN.truesoln(xmin, t(end:-1:1), K, r, sigma, 'CALL');

            case 'EURPUT'
                % Specify the expiry time boundary condition
                w(:,end) = max(K-S,0);
                w(end,:) = 0;
                w(1,:) = (K-xmin)*exp(-r*t(end:-1:1));
                
            case 'AMRPUT'
                % Specify the expiry time boundary condition
                w(:,end) = max(K-S,0);
                w(end,:) = 0;
                w(1,:) = (K-xmin)*exp(-r*t(end:-1:1));
        end

        % LHS
        i = 0:xSteps;
        sigma2 = sigma*sigma;
        a = (1-theta)/(2*dx)*(r.*(xmin+i*dx)-sigma2*(xmin+i*dx).^2/(dx));
        b = 1/dtau+(1-theta)*(r+sigma2*(xmin+i*dx).^2/(dx^2));
        c = (1-theta)/(2*dx)*(-r*(xmin+i*dx) - sigma2*(xmin+i*dx).^2/(dx));
    
        % RHS
        alpha = -theta/(2*dx)*(r*(xmin+i*dx)-sigma*sigma*(xmin+i*dx).^2/(dx));
        beta = 1/dtau-theta*(r+sigma2*(xmin+i*dx).^2/(dx^2));
        gamma = theta/(2*dx)*(r*(xmin+i*dx) + sigma2*(xmin+i*dx).^2/(dx));
         
        A = diag(a(3:xSteps),-1) + diag(b(2:xSteps)) + diag(c(2:xSteps-1),1);
        B = diag(alpha(3:xSteps),-1) + diag(beta(2:xSteps)) + diag(gamma(2:xSteps-1),1);

        % Solve equations
        Cvector = zeros(size(B,2),1);
        
       for j = tSteps:-1:1
                Cvector(1) = alpha(2)*w(1,j+1)-a(2)*w(1,j);
                Cvector(end) = gamma(end)*w(end,j+1)-c(end)*w(end,j);
                w(2:xSteps,j) = Thomas(a(3:xSteps),b(2:xSteps),c(2:xSteps),B*w(2:xSteps,j+1) + Cvector);
       
       switch Type
           case 'AMRPUT'
               for i = 2:xSteps
                   w(i,j) = max(w(i,j),K-S(i));
               end
       end
       end

        % Calculate the option price
        Price = interp1(S,w(:,1),S0);
        Pricevec = w(:,1);
        end
    end
end
  