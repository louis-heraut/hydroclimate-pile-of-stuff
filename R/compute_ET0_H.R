# Copyright 2026 Louis Héraut (louis.heraut@inrae.fr)*1,
#                Louise Mimeau (louise.mimeau@inrae.fr)*1
#
# *1   INRAE, UR RiverLy, Villeurbanne, France
#
# This file is part of hydroclimate-pile-of-stuff R package.
#
# hydroclimate-pile-of-stuff R package is free software: you can
# redistribute it and/or modify it under the terms of the GNU General
# Public License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# hydroclimate-pile-of-stuff R package is distributed in the hope that
# it will be useful, but WITHOUT ANY WARRANTY; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE. See the GNU xGeneral Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with hydroclimate-pile-of-stuff R package.
# If not, see <https://www.gnu.org/licenses/>.


compute_Pa = function (elevation) {
    Pa = 101.3 * ((293 - 0.0065 * elevation)/293)^5.26
    return (Pa)
}

get_gamma = function (elevation) {
    gamma = 0.665 * 10^(-3) * compute_Pa(elevation)
    return (gamma) 
}

compute_ff2m = function (ff) {
    ff2m = ff * 4.87/log(672.58)
    return (ff2m)
}

compute_es = function (T) {
    es = 0.6108 * exp((17.27*T)/(T+273.3))
    return (es)
}

get_es = function (Tmin, Tmax) {
    es = (compute_es(Tmin) + compute_es(Tmax))/2
    return (es)
}

# get_ea = function (Tmin, Tmax, HRmin, HRmax) {
#     ea = (compute_es(Tmin)*HRmax/100 + compute_es(Tmax)*HRmin/100)/2
#     return (ea)
# }

get_ea_mean = function(Tmin, Tmax, HRmean) {
    Tmean = (Tmin + Tmax)/2
    ea = (HRmean/100) * compute_es(Tmean)
    # ea = HRmean/100*(compute_es(Tmin) + compute_es(Tmax))/2
    return (ea)
}

compute_Tmean = function (Tmin, Tmax) {
    Tmean = (Tmin + Tmax) / 2
    return (Tmean)
}

get_Delta = function (Tmin, Tmax) {
    Tmean = compute_Tmean(Tmin, Tmax)
    Delta = 4098 * compute_es(Tmean)/(Tmean + 273.3)^2
    return (Delta)
}

# get_nb_days_in_year = function (date) {
#     NBD = ifelse(lubridate::leap_year(year(date)), 366, 365)
#     return (NBD)
# }

compute_Ra = function (latitude, yearday) {
    Gsc = 0.0820
    phi = latitude * pi / 180
    dr = 1 + 0.033 * cos(2 * pi * yearday / 365)
    delta = 0.409 * sin(2 * pi * yearday / 365 - 1.39)
    omega_s = acos(-tan(phi) * tan(delta))
    Ra = (24*60/pi) * Gsc * dr * (omega_s*sin(phi)*sin(delta) +
                                  cos(phi)*cos(delta)*sin(omega_s))
    return (Ra)
}

compute_Rso = function (elevation, latitude, yearday) {
    Ra = compute_Ra(latitude, yearday)
    Rso = (0.75 + 2*10^(-5)*elevation)*Ra
    return (Rso)
}

# compute_Rn = function (Rg, elevation, latitude, yearday, Tmin, Tmax, HRmin, HRmax) {
#     Rso = compute_Rso(elevation, latitude, yearday)
#     ea = get_ea(Tmin, Tmax, HRmin, HRmax)
#     Rn = 0.77*Rg -
#         (0.34 - 0.14*sqrt(ea)) *
#         4.903e-9 * 
#         ((Tmin + 273.15)^4 + (Tmax + 273.15)^4)/2 *
#         (1.35 * pmin((Rg/Rso), 1) - 0.35)
#     return (Rn)
# }
    
compute_Rs = function (latitude, yearday, Tmin, Tmax, Krs=0.175) {
    Ra = compute_Ra(latitude, yearday)
    # Rs = Krs * Ra * (Tmax - Tmin)**0.5
    Rs = Krs * Ra * sqrt(pmax(Tmax - Tmin, 0))
    return (Rs)
}

compute_Rn_H = function (elevation, latitude, yearday, Tmin, Tmax, HRmean) {
    Rs = compute_Rs(latitude, yearday, Tmin, Tmax)
    Rso = compute_Rso(elevation, latitude, yearday)
    ea = get_ea_mean(Tmin, Tmax, HRmean)
    Rn = 0.77*Rs -
        (0.34 - 0.14*sqrt(ea)) *
        4.903e-9 *
        ((Tmin + 273.15)^4 + (Tmax + 273.15)^4)/2 *
        (1.35 * pmin((Rs/Rso), 1) - 0.35)
    return(Rn)
}

compute_ET0_H = function (elevation, latitude, yearday, Tmin, Tmax, HRmean, ff) {
    Delta = get_Delta(Tmin, Tmax)
    Tmean = compute_Tmean(Tmin, Tmax)
    es = get_es(Tmin, Tmax)
    ea = get_ea_mean(Tmin, Tmax, HRmean)
    gamma = get_gamma(elevation)
    ff2m = compute_ff2m(ff)
    Rn = compute_Rn_H(elevation, latitude, yearday, Tmin, Tmax, HRmean)
    ET0_H = (0.408*Delta*Rn +
             gamma*(900/(Tmean+273))*ff2m *
             pmax(es - ea, 0)) /
        (Delta + gamma*(1 + 0.34*ff2m))
    return (ET0_H)    
}

