% Application of simple constraints
1; % see https://docs.octave.org/interpreter/Script-Files.html
function s=simplebounds(s,Lb,Ub)
  % Apply the lower bound
  ns_tmp=s;
  I=ns_tmp<Lb;
  ns_tmp(I)=Lb(I);

  % Apply the upper bounds
  J=ns_tmp>Ub;
  ns_tmp(J)=Ub(J);
  % Update this new move
  s=ns_tmp;
endfunction

% Test if inequalities hold
% Index function H(g) for inequalities
function H=getH(g)
  if g<=0
      H=0;
  else
      H=1;
  endif
endfunction

% Index function for equalities
function H=getHeq(geq)
  if geq==0
     H=0;
  else
     H=1;
  endif
endfunction

function Z=apply_pressure_vessel_constraints(u)
  Z=0;
  % Penalty constant
  lam=10^15;

  % Inequality constraints
  g(1)= -u(1)+0.0193*u(3);
  g(2) = -u(2) + 0.00954*u(3);
  g(3) = -pi*u(3)^2*u(4)-(4/3)*pi*u(3)^3+1296000;
  g(4) = u(4)-240;
  g

  % No equality constraint in this problem, so empty;
  geq=[];

  % Apply inequality constraints
  for k=1:length(g)
      Z=Z+ lam*g(k)^2*getH(g(k));
  endfor
  % Apply equality constraints
  for k=1:length(geq)
     Z=Z+lam*geq(k)^2*getHeq(geq(k));
  endfor
endfunction

function Z=apply_compression_spring_constraints(u)
  Z=0;
  % Penalty constant
  lam=10^15;

  % Inequality constraints
  g(1) = 1 - (u(3)*u(2)^3)/(71785*u(1)^4);
  g(2) = (4*u(2)^2 - u(1)*u(2))/(12566*(u(2)*u(1)^3 - u(1)^4)) + 1/(5108*u(1)^2) - 1;
  g(3) = 1 - 140.45*u(1)/(u(3)*u(2)^2);
  g(4) = (u(2) + u(1))/1.5 - 1;
  g

  % No equality constraint in this problem, so empty;
  geq=[];

  % Apply inequality constraints
  for k=1:length(g)
      Z=Z+ lam*g(k)^2*getH(g(k));
  endfor
  % Apply equality constraints
  for k=1:length(geq)
     Z=Z+lam*geq(k)^2*getHeq(geq(k));
  endfor
endfunction
