(pwd() != @__DIR__) && cd(@__DIR__) # allow starting app from bin/ dir

using DAREApp
push!(Base.modules_warned_for, Base.PkgId(DAREApp))
DAREApp.main()
