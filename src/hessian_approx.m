% hessian approximation
1; % https://docs.octave.org/interpreter/Script-Files.html
function B_next = sr1_univariate(B_curr, c_next, c_curr, df)

  % calculate B at second point, B1 using BFGS
  s_curr = c_next - c_curr;
  y_curr = df(c_next) - df(c_curr);
  B_next = B_curr + ((y_curr - B_curr*s_curr) * transpose(y_curr - B_curr*s_curr)) / (transpose(y_curr - B_curr*s_curr) * s_curr);

endfunction


function B_next = sr1_bivariate(B_curr, c_next, c_curr, df)

  % calculate B at second point, B1 using BFGS
  s_curr = c_next - c_curr;
  y_curr = df(c_next(1),c_next(2)) - df(c_curr(1),c_curr(2));
  B_next = B_curr + ((y_curr - B_curr*s_curr) * transpose(y_curr - B_curr*s_curr)) / (transpose(y_curr - B_curr*s_curr) * s_curr);

endfunction

function B_next = sr1_trivariate(B_curr, c_next, c_curr, df)

  % calculate B at second point, B1 using BFGS
  s_curr = c_next - c_curr;
  y_curr = df(c_next(1),c_next(2), c_next(3)) - df(c_curr(1),c_curr(2), c_next(3));
  B_next = B_curr + ((y_curr - B_curr*s_curr) * transpose(y_curr - B_curr*s_curr)) / (transpose(y_curr - B_curr*s_curr) * s_curr);

endfunction

function B_next = sr1_quadravariate(B_curr, c_next, c_curr, df)

  % calculate B at second point, B1 using BFGS
  s_curr = c_next - c_curr;
  y_curr = df(c_next(1),c_next(2), c_next(3), c_next(4)) - df(c_curr(1),c_curr(2), c_next(3), c_next(4));
  B_next = B_curr + ((y_curr - B_curr*s_curr) * transpose(y_curr - B_curr*s_curr)) / (transpose(y_curr - B_curr*s_curr) * s_curr);

endfunction


function B_next = bfgs_bivariate(B_curr, c_next, c_curr, df)

  % calculate B at second point, B1 using BFGS
  s_curr = c_next - c_curr;
  y_curr = df(c_next(1),c_next(2)) - df(c_curr(1),c_curr(2));
  B_next = B_curr - (B_curr*s_curr*transpose(s_curr)*B_curr / (transpose(s_curr)*B_curr*s_curr)) + (y_curr*transpose(y_curr) / (transpose(y_curr)*s_curr));

endfunction

function B_next = bfgs_quadravariate(B_curr, c_next, c_curr, df)

  % calculate B at second point, B1 using BFGS
  s_curr = c_next - c_curr;
  y_curr = df(c_next(1),c_next(2), c_next(3), c_next(4)) - df(c_curr(1),c_curr(2), c_next(3), c_next(4));
  B_next = B_curr - (B_curr*s_curr*transpose(s_curr)*B_curr / (transpose(s_curr)*B_curr*s_curr)) + (y_curr*transpose(y_curr) / (transpose(y_curr)*s_curr));

endfunction

