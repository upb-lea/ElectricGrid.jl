

"""
Basic example for a diffrent source
"""
Base.@kwdef mutable struct SolarModule
    I_0 = 2.0381e-10           # Diode seturation
    ni = 1.2                   # Diode ideality factor for Si-mono
    k = 1.3806e-23             # Boltzman constant
    q = 1.6022e-19             # Electron charge
    T_0 = 273                  # Temperature offset
    N_cell = 36                # Number of cells

    G_ref = 1000               # Irradiance at STC
    mu_sc = 1.3e-3             # Coefficient tempreture of short circuit current
    I_ph_ref = 3.11            # Phototcurrent at STC

end;

Base.@kwdef mutable struct SolarArray
    ModuleParameters::SolarModule           # module parameters
    serial = 10
    parallel = 4

end;



function get_I(SolarArr::SolarArray, V, G, T)

    self = SolarArr.ModuleParameters
    function I_photo(self::SolarModule, G, T)
        dT = self.T_0 + T
        I_ph = G / self.G_ref * (self.I_ph_ref + self.mu_sc * dT)
        return I_ph
    end

    function I_diode(self::SolarModule, V, G, T)
        dT = self.T_0 + T
        V_T = self.k * dT / self.q
        I_d = self.I_0 * (exp(V / (self.ni * self.N_cell * SolarArr.serial * V_T)) - 1)
        return I_d
    end

    I = (I_photo(self, G, T) - I_diode(self, V, G, T)) * SolarArr.parallel
    return I
end

function get_V(SolarArr::SolarArray, I, G, T)

    self = SolarArr.ModuleParameters
    function I_photo(self::SolarModule, G, T)
        dT = self.T_0 + T
        I_ph = G / self.G_ref * (self.I_ph_ref + self.mu_sc * dT)
        return I_ph
    end

    I = maximum([0.0, I[1]])
    res = (I - I_photo(self, G, T)) / (self.I_0 * SolarArr.parallel)

    dT = self.T_0 + T
    V_T = self.k * dT / self.q

    if res >= 1
        V = 0
    else
        V = self.ni * self.N_cell * SolarArr.serial * V_T * log(1 - res)
    end

    return V
end
