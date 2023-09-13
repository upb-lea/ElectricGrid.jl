Base.@kwdef mutable struct battery_module
    R = [0.010342 0.0012244
        0.0067316 0.0011396
        0.0051156 0.0012661
        0.0043447 0.0012265
        0.0038826 0.0011163
        0.0034226 0.0009968
        0.003346 0.0011458
        0.0033222 0.001345
        0.0033201 0.0013091
        0.0032886 0.0010986
        0.0028114 0.0010309]     # Resistance for n-th RC pair
    C = [2287.7 11897
        6122 24515
        18460 42098
        20975 44453
        15254 33098
        10440 24492
        13903 32975
        16694 40007
        15784 35937
        12165 26430
        9118 24795]   # Capacitance for n-th RC pair
end;

Base.@kwdef mutable struct battery_block
    battery_module::battery_module
    T_0 = 273.15    # temperature offset
    # n = 100                      # number of RC pairs
    R_0 = 0.016                  # series resistance
    V_0 = [3.5042 3.5136
        3.5573 3.5646
        3.6009 3.6153
        3.6393 3.6565
        3.6742 3.6889
        3.7121 3.7214
        3.7937 3.8078
        3.8753 3.8945
        3.97 3.9859
        4.0764 4.0821
        4.1924 4.193]          # Battery open-circuit voltage
    Q = 0.0                   # Current Ladung in A*s
    SOC = 1.                    # State of Charge (SOC)
    SOC_BP = vec([0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1])
    T_BP = vec([293.15 313.15])
    tau = 1e-4
    mode = nothing

    LT_V0 = LinearInterpolation((SOC_BP, T_BP), V_0)
    LT_R = LinearInterpolation((SOC_BP, T_BP), battery_module.R)
    i_bat_limit = 10000
    P_1h = nothing
    parallel = 3
    serial = 12
    Q_0 = parallel * serial * 26 * 3600         # Nominal Ladung
    v_dc = 0
    v_next = serial * V_0[end, end]
    i_dc = 0.0
end;

function get_V(self::battery_block, I_batt, T)

    self.v_dc = self.v_next
    self.i_dc = I_batt
    dT = self.T_0 + T

    self.Q = clamp(self.Q - I_batt * self.tau, 0, self.Q_0) # Clamp

    self.SOC = clamp(self.Q / self.Q_0, 0, 1)

    if sign(I_batt) == -1
        self.mode = "charge"
    else
        self.mode = "discharge"
    end

    self.P_1h = self.Q / 3600 * self.serial * self.LT_V0(self.SOC, dT)

    I_batt = clamp(I_batt, -self.i_bat_limit, self.i_bat_limit) # Current limitation

    self.v_next = self.serial * (self.LT_V0(self.SOC, dT) - self.LT_R(self.SOC, dT) * I_batt - self.R_0 * I_batt)

    return self.v_dc
end
