#using SpecialFunctions
#using LinearAlgebra
#using Plots
#using PyPlot
#using StatsBase
#using LinearAlgebra
#using FFTW
#using Expokit #for exponentials of sparse matrices
#using LightGraphs
#using GraphRecipes
#using GraphMakie
using DataStructures: CircularBuffer
using DynamicalSystems
using OrdinaryDiffEq
#using CairoMakie

using GLMakie#, Makie

# %% 0. Learn about `Observable`s
# An `Observable` is a mutable container of an object of type `T`.
# `T` can be any type. The value of an `Observable` can then be
# changed on the spot, just like updating any mutable container.
# (This is similar to the base structure `Ref`, if you're familiar)
# The important part here is that `Observable`s can be "listened" to.
# What does this mean...?
o = Observable(1)

l1 = on(o) do val # Create a listener `l1` of observable.
    println("Observable now has value $val")
end

# `l1` is triggered each time the value of `o` is updated.
# (demo in REPL, set `o[] = 2`.)

# We care about `Observable`s because Makie.jl is hooked up
# to this "listener" system. If any plotted element is
# initialized as an observable, then Makie.jl understands this.
# Updating the observable would update the plot values.

# For example:
ox = 1:4
oy = Observable(rand(4))
lw = Observable(2)

fig, ax = lines(ox, oy; linewidth = lw)
#ylims!(ax, 0, 1)

lw[] = 50
oy[] = rand(4)

#display(fig)

# %% 1. Initialize simulation in a stepping manner
# (this can also be done with a pre-run simulation)
# the goal is to create a "step" function which
# once called it progresses the data for one animation frame
const L1 = 1.0
const L2 = 0.9
M = 2
u0 = [π/3, 0, 3π/4, -2]
dp = Systems.double_pendulum(u0; L1, L2)

# Solve diffeq with constant step for smoother curves
diffeq = (alg = Tsit5(), adaptive = false, dt = 0.005)

integ = integrator(dp, u0; diffeq)

function xycoords(state)
    θ1 = state[1]
    θ2 = state[3]
    x1 = L1 * sin(θ1)
    y1 = -L1 * cos(θ1)
    x2 = x1 + L2 * sin(θ2)
    y2 = y1 - L2 * cos(θ2)
    return x1,x2,y1,y2
end

# `integ` is an integrator. `step!(integ)` progresses the integrator
# for one step. `integ.u` is the system state at current step. 
# Then `xycoords` converts the states of the integrator
# to their plottable format. So we can imagine something like 
function progress_for_one_step!(integ)
    step!(integ)
    u = integ.u
    return xycoords(u)
end
# to be our stepping function that returns the new data.

# If we had finite data instead of a forever-running animation, 
# then the "stepping function" would simply be to progress the index `i`
# of the existing data one step forwards...


# %% 2. Initialize the `Observable`s of the animation 
# You need to think of this in advance: what things will to be 
# animated, and what other plotting elements will be static? 
# Animated elements will need to become `Observable`s.

# Here the animated elements will be: balls and rods making the
# double pendulum, and the tail (trajectory) of the pendulum.
x1,x2,y1,y2 = xycoords(u0)
rod   = Observable([Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)])
balls = Observable([Point2f(x1, y1), Point2f(x2, y2)])
# (Remember: the most optimal way to plot 2D things in Makie.jl is to
# give it a vector of `Point2f`, the coordinates for the plot)

# Here we have initialized two _different_ observables, because
# rods and balls will be plotted in a different manner (lines/scatter)

# Next is the observable for the tail
tail = 300 # length of plotted trajectory, in units of `dt`
# The circular buffer datastructure makes making stepping-based
# animations very intuitive
traj = CircularBuffer{Point2f}(tail)
fill!(traj, Point2f(x2, y2)) # add correct values to the circular buffer
traj = Observable(traj) # make it an observable


# %% 3. Plot the `Observable`s and any other static elements
# Before plotting we need to initialie a figure
fig = Figure(); 
display(fig)
# in my experience it leads to cleaner code if we first initialize 
# an axis and populate it accordingly.
ax = Axis(fig[1,1])

# Now we plot the observables _directly_! First the pendulum
lines!(ax, rod; linewidth = 4, color = :purple)
scatter!(ax, balls; marker = :circle, strokewidth = 2, 
    strokecolor = :purple,
    color = :black, markersize = [8, 12]
)

# then its trajectory, with a nice fadeout color
c = to_color(:purple)
tailcol = [RGBAf(c.r, c.g, c.b, (i/tail)^2) for i in 1:tail]
lines!(ax, traj; linewidth = 3, color = tailcol)

# We can also plot now any other static elements
ax.title = "double pendulum"
ax.aspect = DataAspect()
l = 1.05(L1+L2)
xlims!(ax, -l, l)
ylims!(ax, -l, 0.5l)

# %% 4. Create the "animation stepping function"
# Using the functions of step 1, we now define a function
# that updates the observables. Makie.jl understands observable
# updates and directly reflects this on the plotted elements.
function animstep!(integ, rod, balls, traj)
    x1,x2,y1,y2 = progress_for_one_step!(integ)
    rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    balls[] = [Point2f(x1, y1), Point2f(x2, y2)]
    push!(traj[], Point2f(x2, y2))
    traj[] = traj[] # <- important! Updating in-place the value of an
                    # `Observable` does not trigger an update!
end

# %% 5. Test it
for i in 1:1000
    animstep!(integ, rod, balls, traj)
    sleep(0.001)
end

# cool it works. Let's wrap up the creation of the observables
# and plots in a function (just to re-initialie everything)
function makefig(u0)
    dp = Systems.double_pendulum(u0; L1, L2)
    integ = integrator(dp, u0; diffeq...)
    x1,x2,y1,y2 = xycoords(u0)
    rod   = Observable([Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)])
    balls = Observable([Point2f(x1, y1), Point2f(x2, y2)])
    traj = CircularBuffer{Point2f}(tail)
    fill!(traj, Point2f(x2, y2)) # add correct values to the circular buffer
    traj = Observable(traj) # make it an observable
    fig = Figure(); display(fig)
    ax = Axis(fig[1,1])
    lines!(ax, rod; linewidth = 4, color = :purple)
    scatter!(ax, balls; marker = :circle, strokewidth = 2, 
        strokecolor = :purple,
        color = :black, markersize = [8, 12]
    )
    lines!(ax, traj; linewidth = 3, color = tailcol)
    ax.title = "double pendulum"
    ax.aspect = DataAspect()
    l = 1.05(L1+L2)
    xlims!(ax, -l, l)
    ylims!(ax, -l, 0.5l)
    # also return the figure object, we'll ues it!
    return fig, integ, rod, balls, traj
end
    
# %% 6. Save animations to videos
fig, integ, rod, balls, traj = makefig(u0)
frames = 1:200
record(fig, "video.mp4", frames; framerate = 60) do i # i = frame number
    for j in 1:5 # step 5 times per frame
        animstep!(integ, rod, balls, traj)
    end
    # any other manipulation of the figure here...
end # for each step of this loop, a frame is recorded

# %% 7. Interactive application
# Makie.jl has tremendously strong capabilities for real-time
# interactivity. To learn all of this takes time of course,
# and you'll need to consult the online documentation.
# Here we will do two interactions: 1) a play/stop button
# 2) clicking on the screen and getting a new initial condition!

fig, integ, rod, balls, traj = makefig(u0)
# The run button is actually pretty simple, we'll add it below the plot
run = Button(fig[2,1]; label = "run", tellwidth = false)
# This button will start/stop an animation. It's actually surprisingly
# simple to do this. The magic code is:
isrunning = Observable(false)
on(run.clicks) do clicks; isrunning[] = !isrunning[]; end
on(run.clicks) do clicks
    @async while isrunning[]
        isopen(fig.scene) || break # ensures computations stop if closed window
        animstep!(integ, rod, balls, traj)
        sleep(0.001) # or `yield()` instead
    end
end

# `on` an important Observables function when one starts
# doing advanced stuff. It triggers a piece of code once an observable
# triggers its update.

# We'll add one more interactive feature which will trigger once
# we click on the axis. Notice that by default makie performs a zoom
# once one clicks on the axis, so we'll disable this
ax = content(fig[1,1])
Makie.deactivate_interaction!(ax, :rectanglezoom)
# and we'll add a new trigger using the `select_point` function:
spoint = select_point(ax.scene)

# Now we see that we can click on the screen and the `spoint` updates!

# Okay, let's do something useful when it triggers
function θωcoords(x, y)
    θ = atan(y,x) + π/2
    return SVector(θ,0,0,0)
end

on(spoint) do z
    x, y = z
    u = θωcoords(x, y)
    reinit!(integ, u)
    # Reset tail and balls to new coordinates
    x1,x2,y1,y2 = xycoords(u)
    traj[] .= fill(Point2f(x2, y2), length(traj[]))
    traj[] = traj[]
    rod[] = [Point2f(0, 0), Point2f(x1, y1), Point2f(x2, y2)]
    balls[] = [Point2f(x1, y1), Point2f(x2, y2)]
end

# And that's the end of the tutorial!