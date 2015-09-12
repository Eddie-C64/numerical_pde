function diffusion1(N,steps)
%diffusion1(N,steps) solves the diffusion (heat) equation u_t = u_xx 
%    with N+1 gridpoints 1,...,N+1 to time t = steps*dt
%    with Dirichlet BCs u(-1,t)=0=u(1,t) & ICs u(x,0)=1-x^2 
%    Uses TRBDF2 method with fixed timestep

tic
% Note there are N-1 interior points
GAMMA = 2-sqrt(2);
u = zeros(N-1,1); % set u to be an N-1 dimensional column vector
umid = zeros(N-1,1);
x = zeros(N-1,1);
h = 2/N;
%dt = h^2/2; % forward Euler stability limit
dt = N*h^2/2;
e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1);
j = (1:N-1)';
x = -1+h*j;
u = 1-x.^2;

CONST = GAMMA/(2*h^2);
CONST1 = (1-GAMMA)/((2-GAMMA)*h^2);
CONST2 = 1/(GAMMA*(2-GAMMA));
CONST3 = (1-GAMMA)^2/(GAMMA*(2-GAMMA));
for n = 1:steps % timestep loop
    umid = (speye(N-1)-CONST*dt*D2)\((speye(N-1)+CONST*dt*D2)*u); % TR
    u = (speye(N-1)-CONST1*dt*D2)\(CONST2*umid-CONST3*u); % BDF2
    plot([-1; x; 1],[0; u; 0],'r-')
    xlabel('x','FontSize',16)
    ylabel('u','FontSize',16)
    axis([-1 1 0 1.1]);
    M(n) = getframe;
end
toc
movie(M,1,10)
end
