from math import sqrt, pow


Schott	    =   1
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

# The following two formulas can be found here:
#   http://www.zemax.com/os/resources/learn/knowledgebase/how-zemax-calculates-refractive-index-at-arbitrary

def solve_n_air(n_ref, temp, pres):
    return 1 + ((n_ref - 1) * pres) / (1.0 + (temp - 15) * 3.4785e-3)

def solve_n_ref(l):
    return 1 + ( 6432.8 + (2949810*l**2)/(146*l**2 - 1) \
                        + (25540*l**2)/(41*l**2 - 1) ) \
             * 1.0e-8

def glass_index(formula, wave, temp, pres, c):
    if   formula == Schott:
        n_glass = schott(wave, temp, pres, c)
    elif formula == Sellmeier1:
        n_glass = sellmeier1(wave, temp,pres, c)
    elif formula == Herzberger:
        n_glass = herzberger(wave, temp, pres, c)
    elif formula == Sellmeier2:
        n_glass = sellmeier2(wave, temp, pres, c)
    elif formula == Conrady:
        n_glass = conrady(wave, temp, pres, c)
    elif formula == Sellmeier3:
        n_glass = sellmeier3(wave, temp, pres, c)
    elif formula == Handbook1:
        n_glass = handbook1(wave, temp, pres, c)
    elif formula == handbook2:
        n_glass = handbook2(wave, temp, pres, c)
    elif formula == Sellmeier4:
        n_glass = sellmeier4(wave, temp, pres, c)
    elif formula == SchottE1:
        n_glass = schottE1(wave, temp, pres, c)
    elif formula == Sellmeier5:
        n_glass = sellmeier5(wave, temp, pres, c)
    elif formula == SchottE2:
        n_glass = schottE1(wave, temp, pres, c)
    elif formula == SchottE3:
        n_glass = schottE1(wave, temp, pres, c)
    elif formula == Mirror:
        n_glass = -1
    elif formula == Default:
        n_glass = 1.0

    return n_glass / n_air
