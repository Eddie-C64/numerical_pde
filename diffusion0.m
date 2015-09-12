function diffusion0(N,steps)
%diffusion0(N,steps) solves the diffusion (heat) equation u_t = u_xx 
%    with N+1 gridpoints 1,...,N+1 to time t = steps*dt
%    with Dirichlet BCs u(-1,t)=0=u(1,t) & ICs u(x,0)=1-x^2 
%    Uses backward Euler method
%		It also has an option for run time.
%    Try: diffusion0(100,50)
%		
tic
% Note there are N-1 interior points
u = zeros(N-1,1); % set u to be an N-1 dimensional column vector
x = zeros(N-1,1);
h = 2/N;
%dt = h^2/2; % forward Euler stability limit
dt = N*h^2/2;
e = ones(N-1,1);
D2 = spdiags([e -2*e e],[-1 0 1],N-1,N-1);
j = (1:N-1)';
x = -1+h*j;
u = 1-x.^2;
for n = 1:steps % timestep loop
    % u = u+dt*D2*u/h^2; % forward Euler
    u = (speye(N-1)-dt*D2/h^2)\u; % backward Euler
    plot([-1; x; 1],[0; u; 0],'r-')
    xlabel('x','FontSize',16)
    ylabel('u','FontSize',16)
    axis([-1 1 0 1.1]);
    M(n) = getframe;
end
toc
movie(M,1,5)
end
