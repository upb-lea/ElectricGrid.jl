using DrWatson
@quickactivate "JEG"

include(srcdir("classical_control.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

x = [-5.; -1; 0; 2]
y = [-2.; 6; 1; 3]

coef = DividedDiff(x, y)

xnew = -5:0.1:2
ynew = NewtonInterpolation(coef, x, x[4])


print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")