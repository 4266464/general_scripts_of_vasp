from sympy import *
import numpy

x, y, z = symbols('x y z')
# x = symbols('x', positive = True)
# from sympy.abc import x, y
# vars = symbols('x_1:5')

# 展开函数
expand((x + 1)**2)
# 构造分数
Rational(1, 2)
# 表达式
cos(x) + 1
# 注意和 simplify 区分
sympify('x**2 + 2*x + 1')


# 带入值
((x + 1)**2).subs(x, 0)
# float求值
pi.evalf(3)
# numpy化
lambdify(x, sin(x), 'numpy')(numpy.pi / 3)

# 多项式
simplify(sin(x)**2 + cos(x)**2)
factor(x**3 - x**2 + x - 1)
collect(x*y + x - 3 + 2*x**2 - z*x**2 + x**3, x)
cancel((x**2 + 2*x + 1)/(x**2 + x))
apart((4*x**3 + 21*x**2 + 10*x + 12)/(x**4 + 5*x**3 + 5*x**2 + 4*x))

# diff and integrate
diff(x**4, x, 3)
((x + 1)**2).diff(x, 2)
exp(x*y*z).diff(x)

integrate(cos(x), x)
integrate(exp(-x), (x, 0, oo))
integrate(exp(-x**2 - y**2), (x, -oo, oo), (y, -oo, oo))

limit(sin(x)/x, x, 0)
limit(1/x, x, 0, '+')

sin(x).series(x, 0, 4)

# solve
solveset(Eq(x**2 - x, 0), x, domain = S.Reals)

f = symbols('f', cls = Function)
dsolve(Eq(f(x).diff(x, 2) - 2*f(x).diff(x) + f(x), sin(x)), f(x))

# matrix
Matrix([[1, -1], [3, 4], [0, 2]])
Matrix([1, 2, 3])
Matrix([[1], [2], [3]]).T

eye(4)
zeros(4)
ones(4)
diag(1, 2, 3, 4)

M = Matrix([[1, 3], [-2, 3]])
M**2
M**-1
M.det()
M.eigenvals()
M.eigenvects()

# plot
from sympy.plotting import plot
plot(x**2, (x, -2, 2))
plot_implicit(Eq(x**2 + y**2, 1))

from sympy.plotting import plot3d
plot3d(x*exp(-x**2 - y**2), (x, -3, 3), (y, -2, 2))

# latex output
print(latex(integrate(sqrt(x), x)))