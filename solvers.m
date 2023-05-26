% collection of solvers

function c_next = steepest_descent_univariate(c_curr, df, alpha)

  % calculate direction at first point
  p_next = -df(c_curr);
  c_next = c_curr + alpha * p_next;

endfunction


function c_next = quasi_newton_sr1_univariate(c_curr, df, ddf, alpha, B_curr)

  p_curr = -inv(B_curr) * df(c_curr);
  c_next = c_curr + alpha * p_curr

 endfunction


