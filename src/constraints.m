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

function Z=apply_constraints(u)
  Z=0;
  % Penalty constant
  lam=10^15;

  % Inequality constraints
  g(1)= -u(1)+0.0193*u(3);
  g(2) = -u(2) + 0.00954*u(3);
  g(3) = -pi*u(3)^2*u(4)-(4/3)*pi*u(3)^3+1296000;
  g(4) = u(4)-240;

  % No equality constraint in this problem, so empty;
  geq=[];

  % Apply inequality constraints
  for k=1:length(g),
      Z=Z+ lam*g(k)^2*getH(g(k));
  endfor
  % Apply equality constraints
  for k=1:length(geq),
     Z=Z+lam*geq(k)^2*getHeq(geq(k));
  endfor
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
