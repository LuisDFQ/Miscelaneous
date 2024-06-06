import math
from math import *
from sympy import diff, Symbol, sympify

# This script is for facilitating the work in Physics laboratory
# mainly the referent to Error analysis


def mean(Q: list):  # Mean
    if Q==[]:
        raise ValueError('Q must be a non-empty list.')
    return sum(Q)/len(Q)


def variance(n, Q: list):  # Variance
    if len(Q) != n:
        raise ValueError('You must insert n data points.')
    T = []
    for i in range(n):
        T.append((Q[i]-mean(Q))**2)
    return 1/(n-1)*sum(T)


def variance_1m(f, IV: dict, U:dict):  # Variance for 1 measurement (quadrature sum)
    T = []
    obj = {}
    for ob in IV:
        obj[ob] = Symbol(ob)
    for s in obj:
        T.append((diff(f, obj[s]).doit().subs(
            {Symbol(o): IV[o][0] for o in IV}))**2*U[str(s)][0]**2)
    return sum(T).doit()


def variance_m(n, Q):  # Mean's variance
    return 1/n*variance(n, Q)


def e_covariance(n, X1, X2):  # Estimated covariance associated to X1 and X2 (according Cuban Rules)
    T = []
    for i in range(n):
        T.append((X1[i]-mean(X1))*(X2[i]-mean(X2)))
    return 1/n/(n-1)*sum(T)


def uncertainty(variance):  # Standard deviation (or Standard Error of the Mean) or typical deviation (or mean's typical deviation) (depends on type of variance)
    return (variance)**0.5


def variance_c_nocor(f, IV):  # Combined variance for independent (non-correlated) variables
    T = []
    obj = {}
    for ob in IV:
        obj[ob] = Symbol(ob)
    for s in obj:
        if len(IV[str(s)]) == 1:
            T.append((diff(f, obj[s]).doit().subs(
                {Symbol(o): mean(IV[o]) for o in IV}))**2*Unc1[str(s)][0]**2)
        else:
            print("\nFor variable "+str(s)+", it is necessary to know how was measured ...")
            typeB = thumb_rules()
            T.append((diff(f, obj[s]).doit().subs(
                {Symbol(o): mean(IV[o]) for o in IV}))**2*(variance_m(n, IV[str(s)])+typeB**2))
    return sum(T).doit()


def variance_c_cor(f, IV):  # Combined variance for correlated variables (it generalises combined variance for independent variables)
    obj = {}
    Corr = []
    for ob in IV:
        obj[ob] = Symbol(ob)
    T1 = []
    Variab = [var for var in IV.keys()]
    for i in Variab[:-1]:
        T2 = []
        if len(IV[str(i)]) == 1:
            pass
        else:
            for j in Variab[1:]:
                if len(IV[str(j)]) == 1:
                    pass
                else:
                    T2.append(diff(f, obj[j]).doit().subs(
                        {Symbol(o): mean(IV[o]) for o in IV})*e_covariance(n, IV[str(i)], IV[str(j)]))
                    if e_covariance(n, IV[str(i)], IV[str(j)]) != 0 and str(i) != str(j):
                        if [str(i), str(j)] not in Corr and [str(j), str(i)] not in Corr:
                            print('\nVariables ' + str(i) + ' and ' +
                                str(j) + ' are correlated.')
                            Corr.append([str(i), str(j)])
            T1.append(diff(f, obj[i]).doit().subs(
                    {Symbol(o): mean(IV[o]) for o in IV})*sum(T2))
    return variance_c_nocor(f, IV) + 2*sum(T1).doit()


def thumb_rules():
    que = input("\nDo you know what instrument was used for obtaining the measurements? (Yes or no)\n")
    if que.upper() == 'Y' or que.upper() == 'YES':
        kt = input(
            "\nDo you know manufacturer's specifications about absolute tolerance? (Yes or no)\n")
        if kt.upper() == 'Y' or kt.upper() == 'YES':
            try:
                at = float(input('Insert the instrumental absolute tolerance: '))
            except:
                raise TypeError('You must insert a number.')
            unc = at/sqrt(3)
        elif kt.upper() == 'N' or kt.upper() == 'NO':
            print('\n** ATTENTION **: Rules of thumb will be applied.')
            try:
                ti = int(input(
                    '\nWhat kind of instrument did you use for measurement? [1] Digital or [2] Analog\n'))
                if ti != 1 and ti != 2:
                    raise ValueError('Type an existing option.')
            except:
                raise TypeError('You must insert an integer.')
            if ti == 1:
                try:
                    lsd = float(
                        input('\nInsert the value of the least significant digit place displayed: '))
                except:
                    raise TypeError('You must insert a number.')
                unc = 5*lsd
            else:
                try:
                    lsd = float(input('\nInsert the smallest division: '))
                except:
                    raise TypeError('You must insert a number.')
                unc = lsd/2
        else:
            raise ValueError("It's a Yes or No question.")
    elif que.upper() == 'N' or que.upper() == 'NO':
        unc = 0
    else:
        raise ValueError("It's a Yes or No question.")
    return unc


def FunDef(iv):
    F = input(
        '\nType the function that relates measured variables y = f(x,z,...)\nNOTE: If there is constants, type the numeric value.\nPlease write the equation: ')
    F0 = F
    if '^' in F:
        F = F.replace('^', '**')
        F0 = F
    if '!' in F:
        if F[F.index('!')-1] != ')':
            F = F.replace(F[F.index('!')-1]+'!',
                          'factorial('+F[F.index('!')-1]+')')
            F0 = F
        else:
            ind = F.index('!')
            ind0 = ind
            while True:
                ind -= 1
                if F[ind] == '(':
                    break
            F = F.replace(F[ind:ind0]+'!',
                          'factorial('+F[ind+1:ind0-1]+')')
            F0 = F
    for elem in dir(math):
        if elem in F and elem != 'e' and elem != 'tau':
            F = F.replace(elem, '')
    if '(' in F or ')' in F:
        F = F.replace('(', '').replace(')', '')
    if '*' in F:
        F = F.replace('*', ',')
    if '-' in F:
        F = F.replace('-', ',')
    if '+' in F:
        F = F.replace('+', ',')
    if '/' in F:
        F = F.replace('/', ',')
    if '.' in F:
        F = F.replace('.', ',')
    if ',,' in F:
        F = F.replace(',,', ',')
    varL = F.split(',')  # List of measurable variables

    for i in range(2):  # Sweep two times for detect possible additional variables
        for t in varL:
            if t.isdigit():
                varL.remove(t)
            if t == '':
                varL.remove(t)
            while varL.count(t) > 1:
                varL.remove(t)

    ivL = iv.split(',')
    varL.sort()
    ivL.sort()

    if varL != ivL:
        raise ValueError(
            "You've inserted additional variables or have omitted them in equation in comparison to reported ones.")
    return F0


ChOp = input(
    '\nDo you want to find a variable by indirected measurement? (Yes or no)\n')
if ChOp.upper() == 'Y' or ChOp.upper() == 'YES':
    print('** ATTENTION **: Uncertainty Propagation Law will be applied.')
    UPL = True
elif ChOp.upper() == 'N' or ChOp.upper() == 'NO':
    UPL = False
    try:
        ChAn = int(input(
            "\nWhat do you want to calculate for a given variable?\n[1] Mean\n[2] Variance\n[3] Mean's variance\n[4] Uncertainty\n[5] Mean and ucertainty\n\nType an option: "))
    except:
        raise TypeError('You must insert an integer.')
    for op in range(1, 6):
        if ChAn != op:
            if op == 5:
                raise ValueError('Type an existing option.')
            pass
        else:
            break
    Dat = input('\nInsert data points for variable (Ex.: 1,3.0,4,...): ')
    for p in Dat.split(','):
        try:
            float(p)
        except:
            raise TypeError('You must insert numbers.')
        Ld = [float(i) for i in Dat.split(',')]
        def Excep(med):
            if len(med) == 1:
                raise ValueError('You cannot calculate this selection')
    if ChAn == 1:
        print('\nMean: '+str(mean(Ld)))
    elif ChAn == 2:
        Excep(Ld)
        print('\nVariance: '+str(variance(len(Ld), Ld)))
    elif ChAn == 3:
        Excep(Ld)
        print("\nMean's variance: "+str(variance_m(len(Ld), Ld)))
    elif ChAn == 4:
        if len(Ld) == 1:
            print('\nUncertainty: ' + str(thumb_rules()))
        else:
            TR = thumb_rules()
            print('\nUncertainty: ' +
              str(sqrt(variance(len(Ld), Ld)+TR**2)))
            print('\nUncertainty (best estimation): ' +
              str(sqrt(variance_m(len(Ld), Ld)+TR**2)))
    else:
        if len(Ld) == 1:
            print('\nResult: '+str(mean(Ld)) +' +/- '+ str(thumb_rules()))
        else:
            print('Result (best estimated): '+str(mean(Ld)) +
              ' +/- '+str(sqrt(variance_m(len(Ld), Ld)+thumb_rules()**2)))
else:
    raise ValueError("It's a Yes or No question.")
M1 = False
if UPL:
    IV = {}
    iv = input('\nInsert symbols of each measured variable (Ex.: x,y,z,...): ')
    n = int(input('\nHow many measurements did you realize? '))
    if n < 0:
        raise ValueError('You must insert a positive integer.')
    elif n == 1:
        print('\nYou are gonna approximate a variable taking one measurement by measurable ones.')
        Om = input(
            '\n** ATENTION **: Instrumental Uncertainty Method will be applied.\n')
        M1 = True
    if M1:
        # This procedure is according to Smith,W.F.(ed.)(2020). Experimental Physics. Principles and Practice for the Laboratory. CRC Press.
        U = {}
        for a in iv.split(','):
            try:
                meas = float(input('\nInsert the measurement of '+a+': '))
            except:
                raise TypeError('You must insert a number.')
            IV[a] = [meas]
            U[a] = [thumb_rules()]
        f = sympify(FunDef(iv))
        print('Result: '+str(f.doit().subs({Symbol(v): IV[v][0] for v in IV}))+' +/- '+str(
            uncertainty(variance_1m(f, IV, U))))
    else:
        if n == 0 or n == 1:
            raise ValueError(
                'The number of measurements must be bigger than 1.')
        Unc1 = {}
        for a in iv.split(','):
            D = input('\nInsert data points for variable ' +
                      a + ' (Ex.: 1,3.0,4,...): ')
            for u in D.split(','):
                try:
                    float(u)
                except:
                    raise TypeError('You must insert numbers.')
            L = [float(i) for i in D.split(',')]
            IV[a] = L
            if len(L) != n:
                Ans = input(
                    'Was this variable measured one time? (Y)es or (N)o\n')
                if Ans.upper() == 'Y' or Ans.upper() == 'YES':
                    Unc1[a] = [thumb_rules()]
                elif Ans.upper() == 'N' or Ans.upper() == 'NO':
                    raise ValueError("You've not inserted "+str(n) +
                                     ' measurements as you declared before.')
                else:
                    raise ValueError("It's a Yes or No question.")

        f = sympify(FunDef(iv))

        CV = input(
            '\nDo you want to know if correlated variables exist in your data? (Yes or no)\n')
        if CV.upper() == 'Y' or CV.upper() == 'YES':
            if variance_c_cor(f, IV) == variance_c_nocor(f, IV):
                print('\nSummary: No correlated variables found.')
            else:
                print('\nSummary: Correlated variables found.')
        elif CV.upper() == 'N' or CV.upper() == 'NO':
            pass
        else:
            raise ValueError("It's a Yes or No question.")
        try:
            ChC = int(input(
                "\nWhat do you want to calculate?\n[1] Combined Variance (no correlated)\n[2] Combined Variance (correlated)\n[3] Uncertainty\n\nType an option: "))
        except:
            raise TypeError('You must insert an integer.')
        for opt in range(1, 4):
            if ChC != n:
                if n == 3:
                    raise ValueError('Type an existing option.')
                pass
            else:
                break
        if ChC == 1:
            print('Result: '+variance_c_nocor(f, IV))
        elif ChC == 2:
            print('Result: '+variance_c_cor(f, IV))
        else:
            print(40*'-')
            print('\nResult: '+str(f.doit().subs({Symbol(v): mean(IV[v]) for v in IV}))+' +/- '+str(uncertainty(variance_c_nocor(f, IV))))
            print(40*'-')
            print('\nResult (best estimated with correlations): '+str(f.doit().subs({Symbol(v): mean(IV[v]) for v in IV}))+' +/- '+str(
                uncertainty(variance_c_cor(f, IV))))
