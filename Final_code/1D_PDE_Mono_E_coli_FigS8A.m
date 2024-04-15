function pde1
    m=0;
    options=odeset('RelTol',1e-13,'AbsTol',1e-13);
    x = linspace(0,5,100);
    %[0 0.05 0.1 0.15 0.20 0.25 0.30 0.35 0.40 0.45 0.50 0.55 0.60 0.65 0.70 0.75 0.80 0.95 1.00 1.05 1.10 1.15 1.20
    %1.25 1.30 1.35 1.40 1.45
    % 1.50 1.55 1.60 1.65 1.70 1.75 1.80 1.95 2.00];
    t=linspace(0,40,400);
    sol = pdepe(m,@pdefun,@pdeic,@pdebc,x,t,options);
    meth = sol(:,:,1);
    e=sol(:,:,2);
    lactose=sol(:,:,3);

    figure (1)
    surf(x,t,e)
    title('E(x,t)')
    xlabel('Distance x')
    ylabel('Time t')

    figure (2)
    surf(x,t,meth)
    title('Methionine(x,t)')
    xlabel('Distance x')
    ylabel('Time t')

    figure (3)
    surf(x,t,lactose)
    title('Lactose(x,t)')
    xlabel('Distance x')
    ylabel('Time t')
    
    fileID = fopen('z_e.txt','w');
    fprintf(fileID,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n',e);
    fclose(fileID);

    fileID = fopen('z_methionine.txt','w');
    fprintf(fileID,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n',meth);
    fclose(fileID);

    fileID = fopen('z_lactose.txt','w');
    fprintf(fileID,'%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t\n',lactose);
    fclose(fileID);

end

function [c,f,s] = pdefun(x,t,u,dudx) % Equation to solve
    
    if u(1)<1e-13
        u(1)=0;
    end

    if u(2)<1e-13
        u(2)=0;
    end

    if u(3)<1e-13
        u(3)=0;
    end
    
    %D_M=0.01138; %Diffusion rate of nutrient methionine (centimeter per hour).
    %D_su=0.05138; %Diffusion rate of lactose and acetate (centimerter per hour).
    %p_met=0.000106593; %Methionine production rate from S. enterica.
    %p_ac=6.5488E-05; %Acetate production rate from E. coli. 
    %c_E_met=1.45601E-05; %Methionine consumption rate from E.
    %c_E_lac=0.000646081; %Lactose consumption rate from E.
    %c_S=0.002094978; %Acetate consumption rate from S.
    %r_E=0.5452575; %E_g growth rate.
    %r_S=0.2319854; %S_g growth rate.
    %kappa_M=5e-9; %Natural decay rate of methionine.
    %kappa_A=5e-9; %Natural decay rate of acetate.
    %kappa_E=5e-9; %Natural death rate of E_g.
    %kappa_S=5e-9; %Natural death rate of S_g.
    %kappa_L=5e-9; %Nature decay rate of lactose.
    %K_M=1.67546E-07;%Half-saturation methionine concentration for E. coli growth.
    %K_A=1.50E-05; %Half-saturation acetate concentration for S. enterica growth.
    %K_Lac=3.44534E-06; %Half-saturation lactose concentration for E. coli growth.
   
    D_M=0.01; %Diffusion rate of nutrient methionine (centimeter per hour).
    D_su=0.05; %Diffusion rate of lactose and acetate (centimerter per hour).
    p_met=1.56; %Methionine production rate from S. enterica.
    p_ac=1.01; %Acetate production rate from E. coli. 
    c_E_met=0.1; %Methionine consumption rate from E.
    c_E_lac=1; %Lactose consumption rate from E.
    c_S=1; %Acetate consumption rate from S.
    r_E=1; %E_g growth rate.
    r_S=0.5; %S_g growth rate.
    kappa_M=5e-9; %Natural decay rate of methionine.
    kappa_A=5e-9; %Natural decay rate of acetate.
    kappa_E=5e-9; %Natural death rate of E_g.
    kappa_S=5e-9; %Natural death rate of S_g.
    kappa_L=5e-9; %Nature decay rate of lactose.
    K_M=1;%Half-saturation methionine concentration for E. coli growth.
    K_A=1; %Half-saturation acetate concentration for S. enterica growth.
    K_Lac=1; %Half-saturation lactose concentration for E. coli growth.
   
    M=-c_E_met*u(1)/(u(1)+K_M)*(u(3)/(u(3)+K_Lac))*u(2)-kappa_M*u(1);
    E=r_E*u(2)*(u(1)/(u(1)+K_M))*(u(3)/(u(3)+K_Lac))-kappa_E*u(2);
    L=-c_E_lac*u(3)/(u(3)+K_Lac)*u(1)/(u(1)+K_M)*u(2)-kappa_L*u(3);
    
    c = [1; 1; 1];
    %f = [0; 0; 0] .* dudx;
    f = [D_M; 0; D_su] .* dudx;
    s = [M; E; L]; 
    
end
% ---------------------------------------------
function u0 = pdeic(x) % Initial Conditions
    
    u0 = [1000;0;1000];
    
    %for e_where=[2.09  2.50  2.34  1.50  0.49  0.09  3.46  1.80  0.050  1.54  2.08  3.41  0.52  1.63  1.08  1.34  4.93  2.35  3.44  3.20  4.17  0.56  3.68  2.16  4.90] %1st simulation
    %for e_where=[4.82  0.30  1.33  3.45  3.12  0.06  4.57  3.13  4.53  0.33  0.48  0.07  3.22  3.94  2.30  1.29  2.77  0.09  0.97  2.08  0.64  1.92  0.85  0.08  1.47] %2nd simulation
    for e_where=[3.61  0.86  3.30  0.94  2.08  3.55  3.02  3.75  0.43  2.15  0.05  0.32  3.03  1.93  0.19  3.49  2.67  0.84  0.51  0.40  1.69  4.96]%3rd simulation
    %for e_where=[0.75  0.77  0.46  0.69  0.60  0.89  0.62  0.67  0.30  0.36  0.88  0.79  0.84  0.29  0.62  0.35  0.10  0.95  0.34  0.71]
        if x>=e_where && x<=e_where+0.05
            u0(2) = 100;
        end
    end
end
% ---------------------------------------------
function [pl,ql,pr,qr] = pdebc(xl,ul,xr,ur,t) % Boundary Conditions
    pl = [0; 0; 0];
    ql = [1; 1; 1];
    pr = [0; 0; 0];
    qr = [1; 1; 1];
end


