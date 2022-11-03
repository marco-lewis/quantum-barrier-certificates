from complex import Complex

def ComplexVector(name, sz, offset = 0):
        cv = []
        for i in range(offset, offset + sz):
            token = name + '__' + str(i)
            z = Complex(token)
            cv.append(z)
        return cv