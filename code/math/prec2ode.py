class Poly:
    def __init__(self, initial_poly=None) -> None:
        if initial_poly is None:
            self.coef = [0]
        else:
            self.coef = initial_poly[:]
    
    def __add__(self, o):
        ret = Poly()
        n = max(len(o.coef), len(self.coef))
        ret.coef = [0] * n
        for i in range(n):
            ret.coef[i] = (self.coef[i] if i < len(self.coef) else 0) + (o.coef[i] if i < len(o.coef) else 0)
        return ret
    
    def shift(self):
        ret = Poly()
        ret.coef = [0] + self.coef[0:]
        return ret
    
    def deriv(self):
        ret = Poly()
        ret.coef = self.coef[1:]
        for i in range(len(ret.coef)):
            ret.coef[i] *= (i + 1)
        return ret

    def __str__(self) -> str:
        ret = ""
        for i, c in enumerate(self.coef):
            if i > 0:
                ret += f" + {c}x^{i}"
            else:
                ret += f"{c}"
        return ret

class DiffExpr:
    def __init__(self, initial=None) -> None:
        if initial is None:
            self.coef = []
        else:
            self.coef = initial[:]
    
    def __str__(self):
        ret = ""
        for i, p in enumerate(self.coef):
            if i > 0:
                ret += f" + ({p})f^{{{i}}}"
            else:
                ret += f"({p})f"
        return ret if ret != "" else "nothing"

    def __add__(self, o):
        ret = DiffExpr()
        n = max(len(o.coef), len(self.coef))
        ret.coef = [Poly()] * n
        for i in range(n):
            ret.coef[i] = (self.coef[i] if i < len(self.coef) else Poly()) + (o.coef[i] if i < len(o.coef) else Poly())
        return ret

    def shift(self):
        ret = DiffExpr(self.coef)
        for i in range(len(ret.coef)):
            ret.coef[i] = ret.coef[i].shift()
        return ret
    
    def deriv(self):
        ret = DiffExpr()
        ret.coef = [Poly()] * (len(self.coef) + 1)
        for i in range(len(self.coef)):
            ret.coef[i] = ret.coef[i] + self.coef[i].deriv()
            ret.coef[i + 1] = ret.coef[i + 1] + self.coef[i]
        return ret

class PRec:
    def __init__(self, initial=None) -> None:
        if initial is None:
            self.coef = []
        else:
            self.coef = initial[:]

    def convert_to_diffexpr(self):
        ret = DiffExpr()
        for i in range(len(self.coef)):
            for j in range(len(self.coef[i].coef)):
                now = DiffExpr([Poly([self.coef[i].coef[j]])])
                for _ in range(i):
                    now = now.shift()
                
                for _ in range(j):
                    now = now.deriv().shift()
                # print(i, j, now)
                ret = ret + now
        return ret

            

# x = DiffExpr([Poly([0, 0, 1])]).deriv().shift()

prec1 = PRec([Poly([0, 0, 0, -1, 3, -3, 1]), Poly([1, -3, 3, -1]), Poly([-1])]).convert_to_diffexpr()
prec2 = PRec([Poly([0, 0, 0, 1]), Poly([-1])]).convert_to_diffexpr()
# x = Poly().shift().shift()
print(prec2)