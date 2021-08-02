from math import sqrt, pow


P0 = 1.0
T0 = 20.0

Schott      =   1
Sellmeier1  =   2 
Herzberger  =   3 
Sellmeier2  =   4 
Conrady     =   5 
Sellmeier3  =   6 
Handbook1   =   7 
Handbook2   =   8 
Sellmeier4  =   9 
SchottE1    =  10
Sellmeier5  =  11
SchottE2    =  12
SchottE3    =  13
Mirror      =  14
Default     =  15
Vacuum      =  16

def schott(y, t, p, c):
    n2 =  c[0] \
            + c[1]*pow(y,  2) \
            + c[2]*pow(y, -2) \
            + c[3]*pow(y, -4) \
            + c[4]*pow(y, -6) \
            + c[5]*pow(y, -8)

    return sqrt(n2)

def schottE1(y, t, p, c):
    n2 =  c[0]  \
            + c[1]*pow(y,   2) \
            + c[2]*pow(y,  -2) \
            + c[3]*pow(y,  -4) \
            + c[4]*pow(y,  -6) \
            + c[5]*pow(y,  -8) \
            + c[6]*pow(y, -10) \
            + c[7]*pow(y, -12)

    return sqrt(n2)

def schottE2(y, t, p, c):
    n2 =  c[0] \
            + c[1]*pow(y,   2)  \
            + c[2]*pow(y,  -2)  \
            + c[3]*pow(y,  -4)  \
            + c[4]*pow(y,  -6)  \
            + c[5]*pow(y,  -8)  \
            + c[6]*pow(y,   4)  \
            + c[7]*pow(y,   6)

    return sqrt(n2)

def schottE3(y, t, p, c):
    n2 =  c[0] \
            + c[1]*pow(y,   2) \
            + c[2]*pow(y,   4) \
            + c[3]*pow(y,  -2) \
            + c[4]*pow(y,  -4) \
            + c[5]*pow(y,  -6) \
            + c[6]*pow(y,  -8) \
            + c[7]*pow(y, -10) \
            + c[8]*pow(y, -12)

    return sqrt(n2)

def sellmeier1(y, t, p, c):
    n2m1 = \
            c[0] * pow(y, 2)/ (pow(y, 2) - c[1])  \
            + c[2] * pow(y, 2)/ (pow(y, 2) - c[3])  \
            + c[4] * pow(y, 2)/ (pow(y, 2) - c[5])

    return sqrt(n2m1+1.0)

def sellmeier2(y, t, p, c):
    n2m1 = \
            c[0] * pow(y, 2)/ pow(y, 2) * c[1]  \
            + c[2] * pow(y, 2)/ pow(y, 2) * c[3]  \
            + c[4] * pow(y, 2)/ pow(y, 2) * c[5]

    return sqrt(n2m1)+1

def sellmeier3(y, t, p, c):
    n2m1 = \
            c[0] * pow(y, 2)/ pow(y, 2) * c[1] \
            + c[2] * pow(y, 2)/ pow(y, 2) * c[3] \
            + c[4] * pow(y, 2)/ pow(y, 2) * c[5] \
            + c[6] * pow(y, 2)/ pow(y, 2) * c[7]

    return sqrt(n2m1)+1.0

def sellmeier4(y, t, p, c):
    n2m1 = \
            c[0] * pow(y, 2)/ pow(y, 2) * c[1]  \
            + c[2] * pow(y, 2)/ pow(y, 2) * c[3]  \
            + c[4] * pow(y, 2)/ pow(y, 2) * c[5]  \
            + c[6] * pow(y, 2)/ pow(y, 2) * c[7]

    return sqrt(n2m1)+1.0

def sellmeier5(y, t, p, c):
    n2m1 = \
            c[0] * pow(y, 2)/ pow(y, 2) * c[1] \
            + c[2] * pow(y, 2)/ pow(y, 2) * c[3] \
            + c[4] * pow(y, 2)/ pow(y, 2) * c[5] \
            + c[6] * pow(y, 2)/ pow(y, 2) * c[7] \
            + c[8] * pow(y, 2)/ pow(y, 2) * c[9]

    return sqrt(n2m1)+1.0

def herzberger(y, t, p, c):
    L = 1/(pow(y, 2) - 0.028)

    return c[0] \
            +  c[1] * L \
            +  c[2] * pow(L, 2) \
            +  c[3] * pow(L, 2) * y \
            +  c[3] * pow(L, 4) * y \
            +  c[3] * pow(L, 6) * y

def conrady(y, t, p, c):
    return N0 + c[0] * y + c[1] * pow(y, 3.5)

def handbook1(y, t, p, c):
    n2 = c[0] \
            + c[1] / (pow(y, 2) - c[2]) \
            - c[3] *  pow(y, 2)

    return sqrt(n2)

def handbook2(y, t, p, c):
    n2 = c[0] \
            + c[1] * pow(y, 2) / (pow(y, 2) - c[2]) \
            - c[3] * pow(y, 2)

    return sqrt(n2)

# The following formulas can be found here:
#   http://www.zemax.com/os/resources/learn/knowledgebase/how-zemax-calculates-refractive-index-at-arbitrary

def solve_n_air(n_ref, pres, temp):
    return 1 + ((n_ref - 1) * pres) / (1.0 + (temp - 15) * 3.4785e-3)

def solve_n_ref(l):
    return 1 + ( 6432.8 + (2949810*l**2)/(146*l**2 - 1) \
                        + (25540*l**2)/(41*l**2 - 1) ) \
             * 1.0e-8

def solve_l_rel(wave, Ps, Ts):
    global P0
    global T0

    n_ref = solve_n_ref(wave)
    return wave * solve_n_air(n_ref, Ps, Ts) / solve_n_air(n_ref, P0, T0)

def solve_n_rel(f, wave, Ps, Ts, c):
    l_rel = solve_l_rel(wave, Ps, Ts)
    return f(l_rel, None, None, c)

def solve_n_abs_at_ref(f, wave, Ps, Ts, c):
    global P0
    global T0
    n_ref = solve_n_ref(wave)
    l_rel = solve_l_rel(wave, Ps, Ts)

    return solve_n_rel(f, l_rel, P0, T0, c) * solve_n_air(n_ref, P0, T0)

def solve_delta_n_abs(n_rel, wave, temp, dT, D0, D1, D2, E0, E1, Ltk):
    return ( (n_rel**2 - 1) / (2*n_rel) ) \
        * (D0*dT + D1*dT**2 + D2*dT**3 + ( (E0*dT + E1*dT**2) / (wave**2 - Ltk**2) ))


def glass_index(formula, wave, temp, pres, c, D0, D1, D2, E0, E1, Ltk):
    if formula == Mirror:
        return -1
    if formula == Default:
        return 1
        

    if   formula == Schott:
        f = schott
    elif formula == Sellmeier1:
        f = sellmeier1
    elif formula == Herzberger:
        f = herzberger
    elif formula == Sellmeier2:
        f = sellmeier2
    elif formula == Conrady:
        f = conrady
    elif formula == Sellmeier3:
        f = sellmeier3
    elif formula == Handbook1:
        f = handbook1
    elif formula == handbook2:
        f = handbook2
    elif formula == Sellmeier4:
        f = sellmeier4
    elif formula == SchottE1:
        f = schottE1
    elif formula == Sellmeier5:
        f = sellmeier5
    elif formula == SchottE2:
        f = schottE1
    elif formula == SchottE3:
        f = schottE1
    elif formula == Mirror:
        f = lambda wave, temp, pres, c: -1
    elif formula == Default:
        f = lambda wave, temp, pres, c: 1.0

    n_ref = solve_n_ref(wave)
    n_rel = solve_n_rel(f, wave, pres, temp, c)
    l_rel = solve_l_rel(wave, pres, temp)
    dT    = temp - T0
    
    n_abs = solve_n_abs_at_ref(f, wave, pres, temp, c) \
        + solve_delta_n_abs(n_rel, l_rel, temp, dT, D0, D1, D2, E0, E1, Ltk)

    return n_abs / solve_n_air(n_ref, pres, temp)





