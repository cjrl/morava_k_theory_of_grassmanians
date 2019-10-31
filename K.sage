load("gr.sage")

class K:
    def __init__(self, n,Gr):
        self.n = n
        self.Gr = Gr

    def d(self,element):
        return self.Gr.normal_form(self.Gr.Q(self.n)(element))

    def homology_at(self,q):
        # first compute cycles Z
        print "computing homology at q=" + str(q) 
        r = 2^(self.n + 1) - 1
        B1 = self.Gr.additive_basis_for_qth_cohomology(q)
        print "computing d..."
        image = [self.d(b) for b in B1]
        print "d computed"
        print image
        B2 = self.Gr.additive_basis_for_qth_cohomology(q+r)
        print "B2",B2

        cols = []
        for x in image:
            # write x as a vector in basis B2
            cols.append([x.monomial_coefficient(self.Gr.normal_form(b)) for b in B2])
        M = Matrix(FiniteField(2),cols).T
        print "d matrix:"
        print M

        Z = M.right_kernel()
        print Z
        print "kernel dim", Z.rank()
        print "kernel done"

        # now compute boundaries
            
        B3 = self.Gr.additive_basis_for_qth_cohomology(q-r)
        image = [self.d(b) for b in B3]
        image = [x for x in image if x != 0]
        
        # B1.reverse() # do we need this because of the order of the matrix with sage?
        image = [[x.monomial_coefficient(self.Gr.normal_form(b)) for b in B1] for x in image]
        print "normal_forms done"
        B = Z.span(image)
        print "span done"
        print B 

        Q = (Z / B).rank()

        print "quotient done"
        print Q
        return Q

    def rank(self):
        if self.n == 0:
            if self.Gr.d == 1:
                return Mod(self.Gr.c,2)
            if self.Gr.d == 2:
                m = self.Gr.c + 2
                return floor(m / 2)
            
        return sum([self.homology_at(q) for q in range(1,self.Gr.top_cohomological_dim + 1)])

