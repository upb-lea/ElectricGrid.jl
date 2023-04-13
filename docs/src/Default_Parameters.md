# __Default Parameters__

There are a number of default values automatically loaded if they are not specified by the user. These may be altered by provided a dictionary to the "parameters" dictionary containing the following variables.

_______________________________________________________________________________
## Default Grid Values

| Property | Default Value | Unit | Description |
| --- | --- | --- | --- |
| phase | 3 |   | number of phases* |
| v_rms | 230.0 | V | global nominal voltage root mean square reference |
| f_grid | 50.0 | Hz | global nominal frequency reference |
| ramp_end | 2/f_grid | s | The time interval for the sources to ramp up their set-points |
| process_start | 2/f_grid | s | The point in time where changes to reference real and imaginary powers take effect, i.e. when sources in PQ or VSG mode accept "p_set" and "q_set" |
| Δfmax | 0.5 | % | The drop (increase) in source frequency that causes a 100% increase (decrease) in real power (from rated) - determines frequency droop coefficients |
| ΔEmax | 5 | % | The drop (increase) in voltage over the capacitor nearest to a source that causes a 100% increase (decrease) in imaginary power (from rated)  - determines voltage droop coefficients|

#### Note
*Only 3 phase systems are currently supported for sources that are classically controlled.


_______________________________________________________________________________
## Default Source Values

| Property | Default Value | Unit | Description |
| --- | --- | --- | --- |
| pwr |  | VAr | rated apparent power |
| control_type | "classic" |  | The inner control structure. Either PI control ("classic") or Reinforment Learning ("RL") |
| mode | "Synchronverter" | | The higher level control functions.  |
| p_set | 0.0 | W | real power set point |
| q_set | 0.0 | VAi | imaginary power set point |
| v_pu_set | 1.0 | per unit | nominal voltage set point |
| fltr | "LCL" | "L", "LC", or "LCL" | The topology of the inverter's filter |
| v_rip | 0.01537 |  |  | # TODO is the default unit in percentage or a ratio?
| i_rip | 0.15 |  |  |
| vdc | 800.0 | V |  |
| pf | 0.8 |  | power factor |
| Dp | * | N⋅m⋅s/rad | frequency droop coefficient |
| Dq | * | VAi/V | voltage droop coefficient |
| I_kp | * | V/A | current proportional gain |
| I_ki | * | V/A⋅s | current integral gain |
| V_kp | * | A/V | voltage proportional gain |
| V_ki | * | A/V⋅s | voltage integral gain |
| v_δ_set | 0.0 | degrees | reference angle for "Swing" mode |
| τv | 0.002 | s | time constant of the voltage loop |
| τf | 0.002 | s | time constant of the frequency loop |
| Observer | false | true, or false | activate the discrete Luenberger observer (i.e. estimator) for sources with LCL filters |
| L1 | * | H |  |
| R1 | * | Ω  |  |
| L2 | * | F |  |
| R2 | * | Ω |  |
| C | * | F |  |
| R_C | * | Ω |  |
| γ | p_set | W | stochastic process asymptotic mean |
| std_asy | pwr / 8 | W | stochastic process asymptotic standard deviation |
| X₀ | p_set | W | stochastic process initial starting values |
| σ | 0.0 | W | brownian motion scale i.e. ∝ diffusion, volatility parameter |
| Δt | 1/f_grid | s | stochastic process time step |
| k | 0 | 0,1, or 2 | stochastic process interpolation degree |

## Note
*These values are automatically calculated from others.



```julia

```
