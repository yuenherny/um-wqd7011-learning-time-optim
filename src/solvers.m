% collection of solvers
1; % see https://docs.octave.org/interpreter/Script-Files.html
function c_next = steepest_descent_univariate(c_curr, df, alpha)

  p_next = -df(c_curr);
  c_next = c_curr + alpha * p_next;

endfunction


function c_next = quasi_newton_sr1_univariate(c_curr, df, alpha, B_curr)

  p_curr = -inv(B_curr) * df(c_curr);
  c_next = c_curr + alpha * p_curr

 endfunction

function c_next = steepest_descent_bivariate(c_curr, df, alpha)

  p_next = -df(c_curr(1), c_curr(2));
  c_next = c_curr + alpha * p_next;

endfunction


function c_next = quasi_newton_sr1_bivariate(c_curr, df, alpha, B_curr)

  p_curr = -inv(B_curr) * df(c_curr(1), c_curr(2));
  c_next = c_curr + alpha * p_curr

 endfunction

