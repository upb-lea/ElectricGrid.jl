using DrWatson
@quickactivate "MicroGridSimWithRL"

using PyCall
@pyinclude(srcdir("nodeconstructor.py"))

p = Dict()

p["R_source"] = 0.4
p["L_source"] = 2.3e-3
p["C_source"] = 10e-6
p["L_cabel"] = 2.3e-3
p["R_cabel"] = 0.4
p["R_load"] = 14
p["V_dc"] = 300

n = py"NodeConstructor"(2, 2, p)

A, B, C, D = n.get_sys()