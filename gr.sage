class Gr:
   
    # only works for odd d right now
    def __init__(self, d,c):
        "Grassmanian of d-planes in \R^(d+c)"
        self.d =d 
        self.c = c
        self.ws = ["w"+str(i) for i in range(1,d+1)]
        self.wbars = ["wbar" + str(i) for i in range(1,c+1)]
        self.S = PolynomialRing(FiniteField(2),",".join(self.ws)+","+",".join(self.wbars))
        
        # Constructs Cohomology Ring 
        self.R = PolynomialRing(FiniteField(2),",".join(self.ws))
        self.W = [1]+[self.R.gens()[i - 1] for i in range(1,self.d + 1)]
        self.sq_dict = {}

        print "Defining Relations..."
        self.relations = [self.R(self.wbar_expression_in_w(i)) for i in range(c+1,c+1+2)]
        print self.relations
        print "Forming the ideal..."
        self.ideal = self.R.ideal(self.relations)
        print "Taking the quotient..."
        self.Mod2CohomologyRing = self.R.quotient(self.ideal)

        print "Detereming Cohomology..."
        self.top_cohomological_dim = self.d*self.c
        self.H = [self.additive_basis_for_qth_cohomology(q) for q in range(self.top_cohomological_dim+1)]

        # print sum([len(y) for y in self.H])

    def w(self,i):
        if i >= len(self.W):
            return 0
        return self.W[i]
        
    def get_w(self,i):
        if i > self.d:
            return 0
        return self.S.gens()[i-1]

    def get_wbar(self,i):
        if i > self.c:
            return 0
        return self.S.gens()[i-1+self.d]
        
    def additive_basis_for_qth_cohomology(self,q):
        S = range(0,self.c+1)
        exponents = cartesian_product([S for i in range(self.d)]).list()
        exponents = filter(lambda x: sum(x) <= self.c, exponents)
        exponents = filter(lambda x: sum([(i+1)*x[i] for i in range(0,len(x))]) == q, exponents)
        return [self.make_word_with_exponents(ex) for ex in exponents]

    def make_word_with_exponents(self,expontents):
        word = 1
        for i,e in enumerate(expontents):
            word *= self.w(i+1)^e
        return word

    def normal_form(self,expression):
        return self.ideal.reduce(self.R(expression))
        
    def wbar_expression(self,k):
        expression = self.get_w(k)
        for i in range(1,k):
            expression += self.get_w(i)*self.get_wbar(k-i)
        return expression

    def wbar_expression_in_w(self,k):
        
        relations = {}
        for i in range(1,k):
            if self.get_wbar(i) != 0:
                relations[self.get_wbar(i)] = self.wbar_expression(i).substitute(relations)
        # print relations
        expression = self.wbar_expression(k).substitute(relations)
        
        # for i in range(1,k): 
        #     expression = expression.substitute(relations)
        return expression

    # Sq^i(w1^aw2^b)
    # Formulas for Sq^1, Sq^2 form the paper
    # ON GROEBNER BASES AND IMMERSIONS OF GRASSMANN MANIFOLDS G2,n
    # ZORAN Z. PETROVIC and BRANISLAV I. PRVULOVIC ́

    def commutator(self,f,g):
        return lambda x : f(g(x)) + g(f(x))

    def Q(self,n):
        # print n
        if n == 0:
            return lambda x : self.Sq(1,x)
        # print n
        return self.commutator(lambda x : self.Sq(2^(n),x),self.Q(n-1))

    def Sq(self,i,expression):
        if (i,expression) in self.sq_dict:
            return self.sq_dict[(i,expression)]
        answer = self.SqBackend(i,expression)
        self.sq_dict[(i,expression)] = answer
        return answer
    
    def SqBackend(self,i,expression):
        expression = self.R(expression)
        
        if i == 0:
            return expression 
        
        # Sq is linear 
        if len(list(expression)) > 1:
            summation = 0
            for term in list(expression):
                term = term[0]*term[1]
                summation += self.Sq(i,term)
            return summation 

        word = expression

        if word == 0:
            return 0

        if word == 1: # i > |word| since i in dim 0
            return 0

        degree = sum([(x+1)*e for x,e in enumerate(word.exponents()[0])])
        
        # prit "test"
        # print i
        # print degree

        if i > degree:
            return 0

        sub_one = {}
        for w in self.W:
            sub_one[w] = 1

        coefficient = word.subs(sub_one)

        if word in self.W:
            return self.wu_formula(i,self.W.index(word))

        # if i == 1:
        #     return coefficient*(a+b)*self.w1^(a+1)*self.w2^b
        # if i == 2:
        #     return coefficient*(b*self.w1^a*self.w2^(b+1) + binomial(a+b,2)*self.w1^(a+2)*self.w2^b)
        # if b == 0:
        #     return coefficient*binomial(a,i)*self.w1^(a+i)

        x,y = self.split_off_letter(word)

        # print "split",x, y 

        answer = 0
        for j in range(i+1):
            answer = self.R(answer)
            # print "sq " + str(i-j) +" " + str(x)+ " sq " + str(j) + " " + str(y)
            s1 = self.R(self.Sq(i-j,x))
            s2 = self.R(self.Sq(j,y))
            answer += s1*s2
        return answer

    def split_off_letter(self,word):
        exponents = list(word.exponents()[0])
        # print exponents
        first_nonzero = 0
        while exponents[first_nonzero] == 0:
            first_nonzero += 1
        letter = self.W[first_nonzero+1]
        # print letter
        exponents[first_nonzero] = exponents[first_nonzero] - 1
        # print "test"
        new_word = 1
        for i in range(1,len(self.W)):
            # print i
            new_word *= self.W[i]^exponents[i-1]
            
        return letter, new_word

    def wu_formula(self,i,j):
        value = 0
        for t in range(0,i+1):
            value += binomial(j+t-i-1,t)*self.w(i-t)*self.w(j+t)
        return value
    
