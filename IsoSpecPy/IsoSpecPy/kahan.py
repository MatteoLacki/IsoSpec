


class SSummator:
    def __init__(self):
        self.partials = []  # sorted, non-overlapping partial sums

    def add(self, x):
        i = 0
        for y in self.partials:
            if abs(x) < abs(y):
                x, y = y, x
            hi = x + y
            lo = y - (hi - x)
            if lo:
                self.partials[i] = lo
                i += 1
            x = hi
        self.partials[i:] = [x]
    def get(self):
        return sum(self.partials)


class Summator:
    def __init__(self, keep_partials = False):
        self.sum = 0.0
        self.c = 0.0
        self.keep_partials = keep_partials
        self.partials = []

    def add(self, what):
        y = what - self.c
        t = self.sum + y
        c = (t - self.sum) - y
        self.sum = t
        if self.keep_partials:
            self.partials.append(self.sum)

    def get(self):
        return self.sum

    def get_partials(self):
        return self.partials

class PosNegSummator:
    def __init__(self):
        self.pos = Summator()
        self.neg = Summator()
    def add(self, what):
        if what >= 0:
            self.pos.add(what)
        else:
            self.neg.add(-what)
    def get(self):
        return self.pos.get() - self.neg.get()



