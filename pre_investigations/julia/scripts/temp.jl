A = [-173.913  -434.783        0.0        0.0       0.0          0.0;
    100000.0       0.0    -100000.0        0.0       0.0          0.0;
     0.0     434.783    -6260.87       0.0       0.0      -6086.96;
     0.0       0.0          0.0     -173.913  -434.783        0.0;
     0.0       0.0          0.0   100000.0       0.0    -100000.0;
     0.0       0.0      -6086.96       0.0     434.783    -6260.87]

B = [434.783    0.0;
    0.0      0.0;
    0.0      0.0;
    0.0    434.783;
    0.0      0.0;
    0.0      0.0]

function f!(du, u, p, t)
  du[:] = A * u + B * p
end