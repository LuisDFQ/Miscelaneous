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
    return sum(T)


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
    return variance_c_nocor(f, IV) + 2*sum(T1)


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

def Exper_round(a: str, b: str): # Uncertainty report according to Significative-digits's Rule (a> value, b> uncertainty)
    B = [bi for bi in b]
    A = [ai for ai in a]
    count = 0
    pos: bool
    if B[0] != '0':
        if '.' not in B:
            B.append('.')
            B.append('0')
        Bnew = B[:B.index('.')]
        pos = True
    else:
        for i in B[2:]:
            count += 1
            if i != '0':
                break
        Bnew = B[2:2+count]
        pos = False
    digits = len(Bnew)
    if pos == True:
        return '{:.0f} +/- {:1.0e}' . format(round(float(a),-(digits-1)),float(b))
    else:
        if len(A[A.index('.'):]) < len(B[B.index('.'):]):
            STR = '0'*(len(B[B.index('.')+1:]) - len(A[A.index('.')+1:]))
            Amod = [am for am in str(round(float(a),digits))]
            if len(Amod[A.index('.')+1:]) == digits:
                STR = ''
            return f'{round(float(a),digits)}'+STR+' $\pm$ {:1.0e}' .format(float(b), digits = digits)
        else:
            return '{:.>{digits}g} $\pm$ {:1.0e}' .format(round(float(a),digits),float(b), digits = digits)

############################################## Code-User interactions ########################################
Q1 = input('Select one by writting the number of the option: (1) Linear Fit or (2) Error processing\n')
if Q1 == '1':
    pass
    print('You have selected Linear Fit (Y = M*X + N)')
    Xdata = input('\nInsert data points for variable X (Example: 1,3.0,4,...): ')
    for x in Xdata.split(','):
            try:
                float(x)
            except:
                raise TypeError('You must insert numbers.')
    XD = [float(i) for i in Xdata.split(',')]
    Ydata = input('\nInsert data points for variable Y (Example: 1,3.0,4,...): ')
    for y in Xdata.split(','):
            try:
                float(y)
            except:
                raise TypeError('You must insert numbers.')
    YD = [float(i) for i in Ydata.split(',')]
    if len(XD) != len(YD):
        raise ValueError('X and Y must contain the same number of elements or data points.')
    elif len(XD) <=2:
        raise ValueError('You need more than 2 data points.')
    Xsum = sum(XD)
    Ysum = sum(YD)
    Ndp = len(XD)
    X2sum = sum([xd**2 for xd in XD])
    Det = Ndp*X2sum-Xsum**2
    if Det == 0:
        raise ValueError('These data cannot be fitted. Maybe, you have repeated some data points.')
    XYsum = sum([xd*yd for xd,yd in zip(XD,YD)])
    Slope = (Ndp*XYsum-Xsum*Ysum)/Det
    Intercept = (X2sum*Ysum-Xsum*XYsum)/Det
    SumDiff = sum([(yd-Intercept-Slope*xd)**2 for xd,yd in zip(XD,YD)])
    Sfit2 = SumDiff/(Ndp-2)
    UncSlope = sqrt(Sfit2/Det*Ndp)
    UncInter = sqrt(Sfit2/Det*X2sum)
    Ymean = mean(YD)
    Xmean = mean(XD)
    R2 = 1-SumDiff/sum([(yd-Ymean)**2 for yd in YD])
    import matplotlib.pyplot as plt
    import numpy as np
    Xvalues = np.arange(min(XD),max(XD))
    plt.plot(XD, YD, 'o', color = 'black')
    plt.plot(Xmean,Ymean, 'o', color = 'red')
    plt.plot(Xvalues, Slope*Xvalues+Intercept, color ='blue', label = 'Fit')
    plt.xlabel('$X$')
    plt.ylabel('$Y$')
    plt.legend(loc=0)
    plt.title('Linear fit: $R^2$ ='+str(round(R2,2))+', Slope = '+Exper_round(str(Slope),str(UncSlope))+', Intercept = '+Exper_round(str(Intercept),str(UncInter)))
    plt.show()
elif Q1 == '2':
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
        Dat = input('\nInsert data points for variable (Example: 1,3.0,4,...): ')
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
                R = Exper_round(str(mean(Ld)),str(thumb_rules()))
                print('\nResult: '+R.replace('$\pm$',' +/- '))
            else:
                R = Exper_round(str(mean(Ld)),str(sqrt(variance_m(len(Ld), Ld)+thumb_rules()**2)))
                print('Result (best estimated): '+R.replace('$\pm$',' +/- '))
    else:
        raise ValueError("It's a Yes or No question.")
    M1 = False
    if UPL:
        IV = {}
        iv = input('\nInsert symbols of each measured variable (Example: x,y,z,...): ')
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
            Rep = Exper_round(str(f.doit().subs({Symbol(v): IV[v][0] for v in IV})),str(uncertainty(variance_1m(f, IV, U))))
            print('Result: '+Rep.replace('$\pm$',' +/- '))
        else:
            if n == 0 or n == 1:
                raise ValueError(
                    'The number of measurements must be bigger than 1.')
            Unc1 = {}
            for a in iv.split(','):
                D = input('\nInsert data points for variable ' +
                        a + ' (Example: 1,3.0,4,...): ')
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
                    "\nWhat do you want to calculate?\n[1] Combined Variance (no correlated)\n[2] Combined Variance (correlated)\n[3] Uncertainty\n[4] Uncertainty (under correlations)\nType an option: "))
            except:
                raise TypeError('You must insert an integer.')
            for opt in range(1, 5):
                if ChC != n:
                    if n == 4:
                        raise ValueError('Type an existing option.')
                    pass
                else:
                    break
            if ChC == 1:
                print('Result: '+variance_c_nocor(f, IV))
            elif ChC == 2:
                print('Result: '+variance_c_cor(f, IV))
            elif ChC == 3:
                print(40*'-')
                Res1 = Exper_round(str(f.doit().subs({Symbol(v): mean(IV[v]) for v in IV})),str(uncertainty(variance_c_nocor(f, IV))))
                print('\nResult: '+Res1.replace('$\pm$',' +/- '))
            else:
                print(40*'-')
                Res2 = Exper_round(str(f.doit().subs({Symbol(v): mean(IV[v]) for v in IV})),str(uncertainty(variance_c_cor(f, IV))))
                print('\nResult (best estimated with correlations): '+Res2.replace('$\pm$',' +/- '))
else:
    raise ValueError('Only 1 or 2 can be chosen.')
