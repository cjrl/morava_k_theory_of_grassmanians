load("gr.sage")

class K:
    def __init__(self, n,Gr):
        self.n = n
        self.Gr = Gr

    def d(self,element):
        """
        This is the first non-trivial differential in the AHSS for K(n)(-)
        
        # Example 1
        # For K(1) d = Q(1) = Sq2Sq1 + Sq2Sq1
        # hence d(w1) = w1^4
        # in Gr(2,c) with c >= 4 otherwise w1^4 is not a generator for 4th homology
        sage: g = Gr(2,4)
        sage: k = K(1,g) 
        sage: k.d(g.W[1])
        w1^4
        """
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
        r = 2^(self.n + 1) - 1
        qth_M = self.qth_differential_as_matrix(q)

        q_minus_r_basis = self.Gr.additive_basis_for_qth_cohomology(q-r)
        q_minus_r_image = [self.Gr.coordinate_vector_in_qth_basis(self.d(b),q) for b in q_minus_r_basis]
        
        Z = qth_M.right_kernel()
        B = Z.subspace(q_minus_r_image)
        Q = (Z / B)

        return Q.rank()

    def qth_differential_as_matrix(self,q):
        """
        We know the first non-trival differential for K(n)(-) is given by Q(n)(-)
        This will write the d^(q,*) as a matrix 
        
        # Example 1 (cont)
        # In Gr(2,4) we have that d(w1) = w1^4
        # H(Gr(2,4);Z/2) has generators w1w2^2,w1^2w2, and w1^4
        # hence d(w1) as a matrix should be [0,0,1]
        sage: g = Gr(2,4)
        sage: k = K(1,g) 
        sage: k.qth_differential_as_matrix(1)
        [0]
        [0]
        [1]
        """
        r = 2^(self.n + 1) - 1
        q_basis = self.Gr.additive_basis_for_qth_cohomology(q)
        # q_plus_r_basis = self.Gr.additive_basis_for_qth_cohomology(q+r)

        image_under_d = [self.d(b) for b in q_basis]
        
        columns = []
        for x in image_under_d: 
            columns.append(self.Gr.coordinate_vector_in_qth_basis(x,q+r))
        return Matrix(FiniteField(2),columns).T

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
            # E_prime = E.copy()
            for i,z in enumerate(E):
                # print E3
                # print "degree", degree
                # killed = min(E3[i],E3[i+degree])
                # E3[i] -= killed
                # E3[i+degree] -= killed
                if i+degree < len(E) and E[i+degree] != 0:
                    print "i", E[i]
                    print "i+d", E[i+degree]
                    E_prime[i+degree] -= E[i]
                    E_prime[i] = 0
                    print "i'", E_prime[i]
                    print "i+d'", E_prime[i+degree]
            print E_prime
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
    

