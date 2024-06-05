
%A is the new basis, the coefficients are mulitples for legendre
%polynomials 
%Note we are only using the even legendre polynomials, because odd means we
%can get an asymetric result which will look like shit. Higher order = more
%squiggly. Each coefficient is just a number, no real restrictions apply.
%you can add as many numbers to A as you want. 
A = [19.4770   38.9742   36.2418   39.5011   37.4831];
GHE(A)

%it is absolutely possible for these to be too stable, leading to nothing.
%should start with high dx and decrease later on 
%size(df(A,1))
%h(A,1)

%newtonMethod(A)
%steepestDescent(A)
%GHE([13.0170   13.0170   13.0170   13.0170   13.0170   13.0170    0.1000    0.1000])

%steepest descent is a lot better than newtons method for this problem
%because its computationally dirt cheap and we struggle with second
%derriviatives because our nonlinearity is high and resolution low 
function [Af, Tf] = steepestDescent(A0)

    fprintf("Inital temperature:")
    GHE(A0)
    fprintf("Initial A0")
    A0
    
    %REPLACE THIS WITH A WHILE LOOP THAT WILL END 
    for i = 1:5
        DF = transpose(df(A0,0.1));
        
        %now have to solve minimization probelm for step size 
        An = A0 - 0.1*DF;

        while GHE(An) <= GHE(A0)
            A0 = An;
            An = An - 0.1*DF; %replace this to modify the optimization method to sommething less naive 
        end

        GHE(A0)
    end 

    fprintf("Final temperature:")
    Tf =  GHE(A0)
    fprintf("Final A0")
    Af = A0
end 

%newtons method does not work great for the problem to be honest
function [AF,Tf] = newtonMethod(A0)
    
    fprintf("Inital temperature:")
    GHE(A0)
    fprintf("Initial A0")
    A0
        
    dx = 0.1;
    
    %large dx
    for i = 1:5 
        H = inv(h(A0,dx));
        F =  df(A0,dx);
        update = (H*F);
        A0 = A0(1,:)- transpose(update(:,1));
    end 
    
    fprintf("Final temperature:")
    GHE(A0)
    fprintf("Final A0")
    A0

end 

%cacluates the hamiltonian at A with custom dx 
function [H] = h(A,dx)
    derriviatives = zeros(length(A),length(A));

    %a better man would use the symetric proerty of the hamiltionian
    for i = 1:length(A)
        for j = 1:i
            if i ~= j 
                A1 = A;
                A2 = A;
                A3 = A;
                A4 = A;
    
                A1(i) = A(i)+dx;
                A1(j) = A(j)+dx;
    
                A2(i) = A(i)+dx;
                A2(j) = A(j)-dx;
    
                A3(i) = A(i)-dx;
                A3(j) = A(j)+dx;
    
                A4(i) = A(i)-dx;
                A4(j) = A(j)-dx;

                s = (GHE(A1)-GHE(A2)-GHE(A3)+GHE(A4))/(4*dx*dx);
                derriviatives(i,j) = s;
                derriviatives(j,i) = s;
                
            else 
                Am = A;
                Ap = A;
                Am(i) = Am(i)-dx;
                Ap(i) = Ap(i)+dx;
                derriviatives(i,i) = (GHE(Ap) - 2*GHE(A) + GHE(Am))/(dx*dx);
            end 
        end 
    end 

    H = derriviatives; 

end 

%calculates derrivaitive vector at 0 with custom dx 
%GOOD IDEA TO STAY ABOVE DX = 0.1 OR RISK DERRIVIATIVE FAILIND DUE TO
%INSUFFICIENT RESOLUTION OF HEAT SOLVER
function [DF] =df(A,dx)
    derriviatives = zeros(1,length(A));
    for i = 1:length(A)
        Am = A;
        Ap = A;
        Am(i) = Am(i)-dx;
        Ap(i) = Ap(i)+dx;
        derriviatives(i) = (GHE(Ap) - GHE(Am))/(2*dx);
    end
    DF = transpose(derriviatives);
end 

%modified to return exclusively temperature 
function [MT] = GHE(A)
    
    x = linspace(-100,100,201);
    y = zeros(201);
    
    %add all the peaks
    for i = 1:length(A)
        yP = A(i)*legendreP(2*(i-1),x/100);
        y = plus(y(1,:),yP(1,:));
    end
    
    y = max(0,y);
    %add an offset
    y = y+5;
    y = cast(y,"int64");
    %plot(x,y);
    %axis([-100,100,0,100])

    %turn the function into a computational domain
    f = min(y,100); %clip the function to simplify life 
    G = functionToGrid(f);

    %send the computational domain to a heat equation solver 
    [MT] = heatEquationSolver(G,f);

end 

function [MT] = heatEquationSolver(domain,f)
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
    
    imagesc(S)
    
    mass = sum(f);
    steadyStateTemp = mean(S(4,:));
    MT = steadyStateTemp;
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

