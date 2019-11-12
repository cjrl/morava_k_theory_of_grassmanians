load("gr.sage")

class K:
    def __init__(self, n,Gr):
        self.n = n
        self.Gr = Gr

    def d(self,element):
        #  RP = Gr_1
        if self.Gr.d == 1:
            print element.degree()
            if Mod(element.degree(),2) == 0:
                return 0
            # print("test")
            # print 2*2^self.n - 1
            # print(self.Gr.R.gens()[0]^(2*2^self.n - 1))
            # print "mul by", self.Gr.R.gens()[0]^(2*2^self.n-1)
            return self.Gr.normal_form(element * self.Gr.R.gens()[0]^(2*2^self.n-1))
            # return self.Gr.normal_form(element^(2^self.n))

        # Gr_d, d > 1
        return self.Gr.normal_form(self.Gr.Q(self.n)(element))

    def homology_at(self,q):
        # first compute cycles Z
        print "computing homology at q=" + str(q) 
        r = 2^(self.n + 1) - 1
        print "r = ", r
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
            try:
                cols.append([x.monomial_coefficient(self.Gr.normal_form(b)) for b in B2])
            except AttributeError:
                # if x == B2[0]g
                # cols.append([1 for b in B2])
                if x == B1[0]:
                    cols = [[1]]
                else:
                    cols = [[0]]
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
        print image
        try:
            image = [[x.monomial_coefficient(self.Gr.normal_form(b)) for b in B1] for x in image]
        except AttributeError:
            # image = [[(lambda t : 1 if t == x else 0)(b) for b in B1] for x in image]
            if x == B1[0]:
                image = [[1]]
            else:
                image = [[0]]
            print "error handled"
        print image
        print "normal_forms done"

        if Z.dimension() == 0:
            return 0
        
        B = Z.subspace(image)
        print "span done"
        print B 

        Q = (Z / B).rank()

        print "quotient done"
        print Q
        return Q

    def test_d_on_all_gens_of_Gr(self):
        for H in self.Gr.H:
            for q in H:
                print self.test_d_on_element(q)
                if self.test_d_on_element(q) == False:
                    return False
        return True
    
    def test_d_on_element(self,element):
        # print element
        # print self.d(element)
        # print self.d(self.d(element))
        return self.d(self.d(element)) == 0

    def rank(self):
        # dim_Q K(0)(-) Formula
        if self.n == 0:
            m = self.Gr.d + self.Gr.c
            n = 0

            if Mod(self.Gr.d,2) == 0:
                k = self.Gr.d / 2
                if Mod(m,2) == 0:
                    n = m/2
                else:
                    n = (m-1)/2
                print "n=",n
                return binomial(n,k)
            else: # Gr_d d odd
                k = (self.Gr.d - 1)/2
                if Mod(m,2) == 1:
                    n = (m-1)/2
                    print "n=",n
                    return binomial(n,k)
                else:
                    n = (m - 2)/2
                    print "n=",n
                    return 2*binomial(n,k)
            
            
        # the +1 is for the added basepoint
        return 1+sum([self.homology_at(q) for q in range(1,self.Gr.top_cohomological_dim + 1)])

    def third_page_ranks(self):
        return [self.homology_at(q) for q in range(1,self.Gr.top_cohomological_dim + 1)]

    def minimum_possible_rank(self):
        if self.n == 0:
           return self.rank()
        E = self.third_page_ranks()
        E_prime = [0 for x in E] 
        k = 2
        while self.k_differential_possibly_nontrivial(k):
            # print self.k_differential_possibly_nontrivial(k)
            degree = self.kth_diff_degree(k)
            print degree
            E_prime = E
            for i,z in enumerate(E):
                # print E3
                # print "degree", degree
                # killed = min(E3[i],E3[i+degree])
                # E3[i] -= killed
                # E3[i+degree] -= killed
                if i+degree < len(E) and E[i+degree] != 0:
                    E_prime[i] = 0
                    E_prime[i+degree] -= E[i]
            E = E_prime
            k += 1
            # print k
        return sum(E)+1
    
    def k_differential_possibly_nontrivial(self,k):
        kth_diff_degree = self.kth_diff_degree(k)
        # print self.Gr.top_cohomological_dim
        # print kth_diff_degree
        if self.Gr.top_cohomological_dim >= kth_diff_degree:
            return True
        return False

    def kth_diff_degree(self,k):
        return k*(2^(self.n+1)-2)+1
    
    # rough approximation
    def possilbe_higher_diffs(self):
        # kth possible k*(2^(n+1)-2) + 1
        second_possible_diff_degree = 2*(2^(self.n+1)-2)+1
        print second_possible_diff_degree
        if self.Gr.top_cohomological_dim >= second_possible_diff_degree:
            return True
        return False
    

