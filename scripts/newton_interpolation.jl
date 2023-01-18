using DrWatson
@quickactivate "dare"

include(srcdir("Classical_Control.jl"))

print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

x = [-5.; -1; 0; 2]
y = [-2.; 6; 1; 3]

coef = Divided_Diff(x, y)

xnew = -5:0.1:2
ynew = Newton_Interpolation(coef, x, x[4])


print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")