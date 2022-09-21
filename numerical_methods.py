from math import *
import numpy as np
import scipy.misc
from sympy import diff, Symbol, sympify
from matplotlib import pyplot as plt
# This is a compilation of different numerical methods

try:
    e = int(input('Choose method to use:\n\n*Nonlinear Equations:\n[1] Bisection \n[2] Newton \n[3] Secant'
              '\n\n*Interpolation:\n[4] Polynomial \n[5] Lagrange \n\n*Differential Equations:\n[6]'
              'Trapezoid \n[7] Simpson \n[8] Euler \n[9] Taylor \n[10] Verlet \n\nEnter option number: '))
except:
    raise TypeError('Insert an integer.')

for n in range(1, 11):
    if e != n:
        if n == 10:
            raise ValueError('Type an existing option.')
        pass
    else:
        break


def declare_func(p):
    if p == 1:
        ft = input('\nType the equation f such that f(x): ')
        def func(x): return eval(ft)
        return func
    else:
        ft = input('\nType the equation f such that f(x,y): ')
        def func(x, y): return eval(ft)
        return func


def derivative(f, x, h=1E-10):
    # central method
    return (f(x+h)-f(x-h))/(2*h)


def factorial(n):
    if type(n) != int:
        raise TypeError('n must be integer.')
    if n == 0 or n == 1:
        return 1
    else:
        F = 1
        for i in range(1, n+1):
            F *= i
        return F


def nth_der(m, f1, f2, x, y):
    Terms = []
    XS, YS = Symbol('x'), Symbol('y')
    for i in range(m+1):
        Terms.append(factorial(m)/factorial(i)/factorial(m-i) *
                     diff(f1, XS, m-i).doit().subs({XS: x, YS: y})*diff(f2, XS, i).doit().subs({XS: x, YS: y}))
    return sum(Terms)


if e == 1:

    # Bisection method
    # define the x-interval
    try:
        a0 = float(input('\nInsert lower value of x: '))
        b0 = float(input('\nInsert higher value of x: '))
    except:
        raise ValueError('You must insert numbers.')
    if b0 < a0:
        raise ValueError(
            'The last inserted value must be bigger than '+str(a0))
    f = declare_func(1)
    if f(a0)*f(b0) < 0:
        pass
    else:
        raise ValueError('Choose a valid interval such that f(a)*f(b) < 0.')

    def bisection(a, b, func, eps=1E-8):
        while abs(a-b) > eps:
            m = (a+b)/2
            if f(a) * f(m) <= 0:
                b = m
            else:
                a = m
        return m

    print('\nThe solution is {:.8f}'. format(bisection(a0, b0, f)))

elif e == 2:

    # Newton Method
    f = declare_func(1)
    # define interval
    try:
        a = float(input('\nInsert the lower value of interval: '))
        b = float(input('\nInsert the highest value of interval: '))
        if b < a:
            raise ValueError(
                'The last inserted value must be bigger than '+str(a))
        # define initial guess
        x0 = float(input('\nInsert an initial guess: '))
    except:
        raise ValueError('You must insert numbers.')
    if x0 < a:
        raise ValueError('Initial guess must be bigger than '+str(a))

    def newton(x, f, eps=1E-8):
        while abs(f(x)) > eps:
            x -= f(x)/derivative(f, x)
        if x < a:
            print('There is no solution in the interval ' +
                  '['+str(a)+','+str(b)+']')
        else:
            return x
    if newton(x0, f) != None:
        print('The solution for f(x) is {:.8f}'. format(newton(x0, f)))

elif e == 3:

    # Secant Method
    f = declare_func(1)
    # define interval
    try:
        a = float(input('\nInsert the lower value of interval: '))
        b = float(input('\nInsert the highest value of interval: '))
    except:
        raise ValueError('You must insert numbers.')
    if b < a:
        raise ValueError('The last inserted value must be bigger than '+str(a))

    def secant(x1, x0, f, eps=1E-8):
        Wh = True
        while abs(f(x1)) > eps and abs(f(x0)) > eps:
            x1 -= (x1 - x0)/(f(x1) - f(x0))*f(x1)
            x0 = x1 - (x1 - x0)/(f(x1) - f(x0))*f(x1)
            if x0 > b or x1 < a:
                print('There is no solution in the interval ' +
                      '['+str(a)+','+str(b)+']')
                Wh = False
                break
        if Wh:
            if abs(f(x1)) < eps:
                return x1
            else:
                return x0

    if secant(a, b, f) != None:
        print('The solution for f(x) is {:.8f}'. format(secant(a, b, f)))

elif e == 4:

    # Polynomial Interpolation
    option = int(input(
        'What do you want to do? (1) Interpolation or (2) Integration Approximation \n'))
    if option != 1 and option != 2:
        raise ValueError('Choose between 1 and 2.')
    if option == 1:
        # For interpolating
        '''Given a data array, we must to associate it a polynomial function'''
        data = input(
            '\nInsert data points (separating by semicolon) (Example: x1,y1;x2,y2;...): ')

        def poly_interp(data):
            try:
                coord = [k for k in data.split(';')]
                d = len(coord)-1
                print(
                    'For interpolating, a {:d}-degree polynomial function will be used.'. format(d))
                points = [c.split(',') for c in coord]
            except:
                raise ValueError('You must insert correctly the data points.')
            Xs = [float(ls[0]) for ls in points]
            Ys = [float(ls[1]) for ls in points]
            obj = {}
            for x in Xs:
                obj['LS'+str(Xs.index(x))] = []
                for i in range(len(points)):
                    if i == 0:
                        obj['LS'+str(Xs.index(x))].append(1)
                    else:
                        obj['LS'+str(Xs.index(x))].append(x**i)
            LS = np.array([obj[L] for L in obj])
            RS = np.array(Ys)
            Sols = np.linalg.solve(LS, RS)
            print(Sols)
            Xn = np.arange(min(Xs)-1, max(Xs)+1, 1E-4)
            Yn = []
            for p in Xn:
                S = []
                for g in range(d+1):
                    S.append(Sols[g]*p**g)
                Yn.append(sum(S))
            plt.plot(Xs, Ys, 'ro', Xn, Yn, '-b')
            plt.xlabel('$x$')
            plt.ylabel('$y$')
            plt.show()
        poly_interp(data)

    else:    
        # For polynomial integration
        pass

elif e == 5:

    # Lagrange Interpolation
    '''Given a data array, we must to associate it a polynomial function'''
    data = input(
        '\nInsert data points (separating by semicolon) (Example: x1,y1;x2,y2;...): ')

    def lag_interp(data):
        try:
            coord = [k for k in data.split(';')]
            d = len(coord)-1
            print(
                'For interpolating, a {:d}-degree polynomial function will be used.'. format(d))
            points = [c.split(',') for c in coord]
        except:
            raise ValueError('You must insert correctly the data points.')
        nodes = [float(ls[0]) for ls in points]
        Ys = [float(ls[1]) for ls in points]

        def L(x, i):
            L = 1
            for e in nodes:
                if e != nodes[i]:
                    L *= (x-e)/(nodes[i]-e)
            return L

        S = []
        obj = {}
        Xn = np.arange(min(nodes)-1, max(nodes)+1, 1E-4)
        for i in range(d+1):
            obj['T'+str(i)] = []
            for x in Xn:
                obj['T'+str(i)].append(Ys[i]*L(x, i))
        Poly = []
        for l in range(len(Xn)):
            Poli = []
            for i in range(d+1):
                Poli.append(obj['T'+str(i)][l])
            Poly.append(sum(Poli))
        plt.plot(nodes, Ys, 'ro', Xn, Poly, '-b')
        plt.xlabel('$x$')
        plt.ylabel('$y = p(x)$')
        plt.show()
    lag_interp(data)

elif e == 6:

    # Trapezoid Method
    f = declare_func(1)
    try:
        a = float(input('\nInsert the lower value of interval: '))
        b = float(input('\nInsert the highest value of interval: '))
        if b < a:
            raise ValueError(
                'The last inserted value must be bigger than '+str(a))
        N = float(input(
            '\nInsert the number of subdivisions to consider (higher N, more precision): '))
    except:
        raise ValueError('You must insert numbers.')

    def trapezoid(a, b, N):
        subInt = []
        step = (b-a)/N
        for i in np.arange(a, b, step):
            subInt.append((f(i)+f(i+step))/2*step)
        return sum(subInt)

    print('The resut of integration is {:f}'. format(trapezoid(a, b, N)))

elif e == 7:

    # Simpson Method
    f = declare_func(1)
    try:
        a = float(input('\nInsert the lower value of interval: '))
        b = float(input('\nInsert the highest value of interval: '))
        if b < a:
            raise ValueError(
                'The last inserted value must be bigger than '+str(a))
        N = float(input(
            '\nInsert the number of subdivisions to consider (higher N, more precision): '))
    except:
        raise ValueError('You must insert numbers.')

    def simp(f, a, b, N):
        subInt = []
        step = (b-a)/N
        for i in np.arange(a, b, step):
            subInt.append(step/6*(f(i)+4*f((2*i+step)/2)+f(i+step)))
        return sum(subInt)

    print('The resut of integration is {:f}'. format(simp(f, a, b, N)))

elif e == 8:

    # Euler Integration Method
    '''Given a differential equation y'(x) = f(x,y)
    with initial condition y(x0) = y0'''
    f = declare_func(2)
    try:
        x0 = float(input('\nInsert initial value of x: '))
        y0 = float(input('\nInsert initial value of y: '))
        b = float(input('\nInsert the highest value of x: '))
        if b < x0:
            raise ValueError('The value must be bigger than '+str(x0))
    except:
        raise ValueError('You must insert numbers.')

    def euler(f, x0, y0, b, step=1E-3):
        Sol = []
        X = []
        X.append(x0)
        Sol.append(y0)
        for x in np.arange(x0, b+step, step):
            if x != x0:
                X.append(x)
                Sol.append(Sol[-1] + step*f(x, Sol[-1]))
        plt.plot(X, Sol, '-b')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.show()
    euler(f, x0, y0, b)

elif e == 9:

    # Taylor Integration Method
    '''Given a differential equation y^(m)(x) = f(x,y)
    where m is the m-th partial derivative respect to x'''
    f = declare_func(2)
    try:
        m = int(input('\nInsert derivative order: '))
        x0 = float(input('\nInsert initial value of x: '))
        y0 = float(input('\nInsert initial value of y: '))
        b = float(input('\nInsert the highest value of x: '))
        if b < x0:
            raise ValueError('The value must be bigger than '+str(x0))
    except:
        raise ValueError('You must insert numbers.')

    def taylor(m, f, x0, y0, b, step=1E-3):
        Y = []
        Y.append(y0)
        X = np.arange(x0, b+step, step)

        ML = input('\nCan the declared function be expressed as a function product such that'
                   'f(x,y) = f1(x,y) times f2(x,y)? (Y) or (N)\n')
        if ML != 'Y' and ML != 'N' and ML != 'y' and ML != 'n':
            raise ValueError("It's a YES or NO question.")
        elif ML == 'Y' or ML == 'y':
            x, y = Symbol('x'), Symbol('y')
            f1 = input('\nType first function: ')
            f2 = input('\nType second function: ')
            if not 'x' in f1 and not 'x' in f2:
                raise TypeError('Functions must be depend on "x".')
            elif sympify(f1)*sympify(f2) != f(x,y):
                raise TypeError('Bad factorization. f1 * f2 != f')
            leib_met = True
        else:
            leib_met = False
            if m % 2 != 0:
                o = m
            else:
                o = m + 1

        for x in X:
            if x == x0:
                continue
            else:
                T = []
                for n in range(m):
                    if leib_met:
                        T.append(step**(n+1)*nth_der(n, f1, f2, x, Y[-1]) /
                                 factorial(n+1))
                    else:
                        T.append(step**(n+1)*scipy.misc.derivative(f, x,
                                 step, n, args=(Y[-1],), order=o)/factorial(n+1))
                Y.append(Y[-1] + sum(T))
        plt.plot(X, Y, '-b')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.show()
    taylor(m, f, x0, y0, b)

elif e == 10:

    # Verlet Integration Method
    '''Given a differential equation y''(x) = f(x,y)
    with initial condition y(x0) = y0'''
    f = declare_func(2)
    try:
        x0 = float(input('\nInsert initial value of x: '))
        y0 = float(input('\nInsert initial value of y: '))
        y0p = float(input("\nInsert initial value of y': "))
        b = float(input('\nInsert the highest value of x: '))
        if b < x0:
            raise ValueError('The value must be bigger than '+str(x0))
    except:
        raise ValueError('You must insert numbers.')

    def verlet(f, x0, y0, y0p, b, step=1E-3):
        Sol = []
        X = []
        Yp = []
        Sol.append(y0)
        X.append(x0)
        Yp.append(y0p)
        for x in np.arange(x0, b+step, step):
            if x != x0:
                X.append(x)
                Sol.append(Sol[-1] - Yp[-1]*step - 0.5 *
                           f(X[-1], Sol[-1])*step**2)
        plt.plot(X, Sol, '-b')
        plt.xlabel('$x$')
        plt.ylabel('$y$')
        plt.show()
    verlet(f, x0, y0, y0p, b)
