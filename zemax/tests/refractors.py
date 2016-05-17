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


def glass_index(formula, wave, temp, pres, c):
    if   formula == Schott:
        return schott(wave, temp, pres, c)
    elif formula == Sellmeier1:
        return sellmeier1(wave, temp,pres, c)
    elif formula == Herzberger:
        return herzberger(wave, temp, pres, c)
    elif formula == Sellmeier2:
        return sellmeier2(wave, temp, pres, c)
    elif formula == Conrady:
        return conrady(wave, temp, pres, c)
    elif formula == Sellmeier3:
        return sellmeier3(wave, temp, pres, c)
    elif formula == Handbook1:
        return handbook1(wave, temp, pres, c)
    elif formula == handbook2:
        return handbook2(wave, temp, pres, c)
    elif formula == Sellmeier4:
        return sellmeier4(wave, temp, pres, c)
    elif formula == SchottE1:
        return schottE1(wave, temp, pres, c)
    elif formula == Sellmeier5:
        return sellmeier5(wave, temp, pres, c)
    elif formula == SchottE2:
        return schottE1(wave, temp, pres, c)
    elif formula == SchottE3:
        return schottE1(wave, temp, pres, c)
    elif formula == Mirror:
        return -1
    elif formula == Default:
        return 1

    return 0.0
