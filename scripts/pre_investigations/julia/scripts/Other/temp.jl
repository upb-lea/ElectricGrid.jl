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


Base.@kwdef mutable struct hello
  a = 0
  b = a > 2 ? 200 : 0
end





function aaaa()
  a = rand(30_000, 30_000)

  b = maximum(a)
  return b
end

@time aaaa()

@time begin
  aaaa()
end





using CUDA

@time begin
  CUDA.@sync begin
    aaaa()
  end
end