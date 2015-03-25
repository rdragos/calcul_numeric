
function [] = main()
  ALL_TESTS = [10; 20; 50; 100];
  source metoda_rotatiilor.m
  for idx = 1:size(ALL_TESTS)(1)
    %equation matrix
    printf("Running on test %d \n", ALL_TESTS(idx));
    dimension = ALL_TESTS(idx);
    A = eye(ALL_TESTS(idx)) * 2;
    for i = 1:(dimension - 1)
      A(i, i + 1) = 1;
      A(i + 1, i) = 1;
    end
    %column vector
    B = ones(ALL_TESTS(idx))(:,1);
    rotating_solve(A, B);
  endfor
endfunction


main()
