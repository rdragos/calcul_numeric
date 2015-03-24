
function metoda_rotatii(A, B)
  %no idea why this works
	how_many = 0;

  num_rows = num_columns = size(A)(1);
  toy_matrix = A;
  printf("IV Matrix");

	while(norm(toy_matrix, "fro") > 1e-10)
    norm(toy_matrix, "fro");
    how_many += 1;
    mx_matrix = 0;

    px = 0; py = 0;

    for i = 1:num_rows
      for j = 1:num_columns
        if mx_matrix < toy_matrix(i,j)
          mx_matrix = toy_matrix(i,j);
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
          partial_copy(px, j) = partial_copy(j, px) = cos_theta * toy_matrix(px, j) + sin_theta * toy_matrix(py, py);
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
  V;
  lambda;
  toy_matrix;
  A * V;
  norm(toy_matrix, "fro")

endfunction



test = [10];
for idx = 1:size(test)(1)
	%equation matrix
  dimension = test(idx);
	A = eye(test(idx)) * 2;
  for i = 1:(dimension - 1)
    A(i, i + 1) = 1;
    A(i + 1, i) = 1;
  end
	%column vector
	B = ones(test(idx))(:,1);
  X = A \ B
  A * X - B
	metoda_rotatii(A, B);
endfor
