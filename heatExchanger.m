
%inputs mean: 
% x1_ is a term controlling A, steepness of peak.
% x2_ is a term controlling B, hight of peak.
% nx is the number of peaks 
generateHeatExchanger(0,0,0.1,-.01,0,100,10);

function [M,T] = generateHeatExchanger(x11,x12,x13,x21,x22,x23,nx)
    
    x = linspace(-100,100,201);
    y = zeros(201);
    
    %add all the peaks
    for i = 1:nx
        dx = i*201/(nx+1);
        X = -100+dx;
        
        a = getA(x11,x12,x13,X);
        b = getB(x21,x22,x23,X);
    
        yP = getPeak(a,b,X,x);
        y = plus(y(1,:),yP(1,:));
    end
    
    %add an offset
    y = y+5;
    y = cast(y,"int64");
    %plot(x,y);
    %axis([-100,100,0,100])

    %turn the function into a computational domain
    f = min(y,100); %clip the function to simplify life 
    G = functionToGrid(f);

    %send the computational domain to a heat equation solver 
    [M,T] = heatEquationSolver(G,f)

end 

function [mass,steadyStateTemp] = heatEquationSolver(domain,f)
    %solve the steady state heat eqation. Easy... right...please..
    E = zeros(101,203,6); %ae aw an as sp su 
    Tau = 10;
    Conv = .0001;
    Tinf = 30; %make sure the source temp is like 100 and the initial guess is high 

    %iteration 1: set up the equation for each point.
    %fancy loop boundaries ensure i+-1, j+1 always exist 
    for i = 2:202
        for j = 2:f(i-1)

            %north boundary
            if (domain(j+1,i) == 0) 
                E(j,i,5) = E(j,i,5) - (Conv/Tau);
                E(j,i,6) = E(j,i,6) + (Conv/Tau) * Tinf;
            else 
                E(j,i-1,3) = Tau; %add to an
            end 
            
            %east boundary
            if (domain(j,i+1) == 0) 
                E(j,i,5) = E(j,i,5) - (Conv/Tau);
                E(j,i,6) = E(j,i,6) + (Conv/Tau) * Tinf;
            else 
                E(j,i,1) = Tau; %add to ae
            end  
            
            % west boundary 
            if (domain(j,i) == 0) 
                E(j,i,5) = E(j,i,5) - (Conv/Tau);
                E(j,i,6) = E(j,i,6) + (Conv/Tau) * Tinf;
            else 
                E(j,i,2) = Tau; %add to aw
            end 

            E(j,i,4) = Tau; %always add to south 

        end 
    end 

    %imagesc(equations(:,:,6));

    %at this point the equations can actually be solved 
    S = ones(101,203)*50;
    S(1,2:202) = 10000; %this adds the source
    
    for z=1:2000

        %jacobi method goes in here!
        for i = 2:202
            for j = 2:f(i-1)
                ap = (E(j,i,1)+E(j,i,2)+E(j,i,3)+E(j,i,4)-E(j,i,5));
                S(j,i) = (E(j,i,1)*S(j,i+1) +E(j,i,2)*S(j,i-1) + E(j,i,3)*S(j+1,i) + E(j,i,4)*S(j-1,i) + E(j,i,6)) / ap;
            end 
        end 
        
    end
    
    %imagesc(S)
    
    mass = sum(f)
    steadyStateTemp = mean(S(4,:))

end 

function[G] = functionToGrid(f)
    grid = zeros(101,203);
    
    for i = 2:202
        for j = 1:f(i-1)
            grid(j,i) = 1;
        end
    end
    
    %imagesc(grid);
    G = grid; 
end 

%peak creation function
function [yPeak] = getPeak(a,b,c,x)
    yPeak = b*2.^(-a*(x-c).^2);
end 

function [a] = getA(x1,x2,x3,X)
    a = max([(x1*X^2)+(x2*X)+x3,0]);
end 

function [b] = getB(x1,x2,x3,X)
    b = max([(x1*X^2)+(x2*X)+x3,0]);
end 

