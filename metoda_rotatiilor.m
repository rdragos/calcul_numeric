%preventing Octave recognizing as a function file
1;

function [retval] = my_norm(A)
  num_rows = size(A)(1);
  num_columns = size(A)(2);

  retval = 0;

  for i = 1:num_rows
    for j = 1:num_columns
      if i ~= j
        retval += A(i,j) * A(i, j);
      end
    end
  end
  retval = sqrt(retval);
endfunction

function [] = rotating_solve(A, B)
  %no idea why this works
	how_many = 0;

  num_rows = num_columns = size(A)(1);
  toy_matrix = A;
	while(my_norm(toy_matrix) > 1e-10)
    how_many += 1;
    mx_matrix = 0;

    px = 0; py = 0;

    for i = 1:num_rows
      for j = i+1:num_columns
        if mx_matrix < abs(toy_matrix(i,j))
          mx_matrix = abs(toy_matrix(i,j));
          px = i;
          py = j;
        end
      end
    end
    theta = pi / 4;
    if toy_matrix(px, px) ~= toy_matrix(py, py)
      %cica scoatem arctangenta
      theta = 0.5 * atan(2 * toy_matrix(px, py) / (toy_matrix(px, px) - toy_matrix(py, py)) );
    end
    cos_theta = cos(theta);
    sin_theta = sin(theta);

    partial_copy = toy_matrix;

    for j = 1: num_columns
      if (j ~= px)
        if (j ~= py)
          partial_copy(px, j) = partial_copy(j, px) = cos_theta * toy_matrix(px, j) + sin_theta * toy_matrix(py, j);
          partial_copy(py, j) = partial_copy(j, py) = -sin_theta * toy_matrix(px, j) + cos_theta * toy_matrix(py, j);
        end
      end
    end
    % such is eigen values
    partial_copy(px, px) = (cos_theta ** 2) * toy_matrix(px, px) + 2 * cos_theta * sin_theta * toy_matrix(px, py) + (sin_theta ** 2) * toy_matrix(py, py);
    partial_copy(py, py) = (sin_theta ** 2) * toy_matrix(px, px) - 2 * cos_theta * sin_theta * toy_matrix(px, py) + (cos_theta ** 2) * toy_matrix(py, py);

    partial_copy(px, py) = partial_copy(py, px) = 0;

    toy_matrix = partial_copy;
  endwhile

  %toy_matrix
  [V, lambda] = eig(A);
  printf("Octave gives me; then custom method gives me; \n");
  lambda = sort(diag(lambda));
  fake_lambda = sort(diag(toy_matrix));
  for i = 1:size(lambda)(1)
    printf("%f ----------------- %f\n", lambda(i), fake_lambda(i));
  end
  A * V;

endfunction

