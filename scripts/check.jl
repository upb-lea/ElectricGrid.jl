print("\n...........o0o----ooo0o0ooo~~~  START  ~~~ooo0o0ooo----o0o...........\n\n")

L1 = 0.002
R1 = 0.04
L2 = 0.002
R2 = 0.05
R_C = 0.09
C = 0.0003

f = 50
ω = 2π*f

v1 = 1.02*230*sqrt(2)
v2 = 335.881*exp(-1.125706im*π/180)

x1 = ω*L1
x2 = ω*L2
xc = -1/(ω*C)

z1 = R1+ 1im*x1
z2 = R2+ 1im*x2
zc = R_C + 1im*xc

z1c = z1 + zc
z2c = z2 + zc

a = v2/zc
b = -z2c/zc
c = v1/zc
d = -z1c/zc

I1 =(a + b*c)/(1 - b*d)

I2 = (v1 - I1*z1c)/zc

i1 = (v2 - I2*z2c)/zc

Ic = I1 + I2

-v1 + I1*z1 + Ic*zc

-v2 + I2*z2 + Ic*zc

println("|I_inv| = ", abs(I1), " angle = ", angle(I1)*180/pi)
println("|I_poc| = ", abs(I2), " angle = ", 180 + angle(I2)*180/pi)
println("|v_cap| = ", abs(Ic*zc), " angle = ", angle(Ic*zc)*180/pi)

print("\n...........o0o----ooo0o0ooo~~~  END  ~~~ooo0o0ooo----o0o...........\n")


a = sqrt(100000)
b = 101 
c = -300

if a^2 < b^2 + c^2

    k = 3 # always an integer (we have infinite solutions)
    δ₁ = 2*atan((b + (- a^2 + b^2 + c^2)^(1/2))/(a + c)) + 2*pi*k # set of solutions 1
    δ₂ = 2*atan((b - (- a^2 + b^2 + c^2)^(1/2))/(a + c)) + 2*pi*k # set of solutions 2

    eq1 = b*sin(δ₁) + c*cos(δ₁) - a
    eq2 = b*sin(δ₂) + c*cos(δ₂) - a

    if round(eq1, digits = 3) != 0.0 || round(eq2, digits = 3) != 0.0

        println("No viable solution found")
    end

else
    println("No viable solution found")
end

println("δ₁ = ", round((δ₁*180/pi)%(2π), digits = 3), " degrees")
println("δ₂ = ", round((δ₂*180/pi)%(2π), digits = 3), " degrees")