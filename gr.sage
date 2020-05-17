class Gr:
    """Grassmanian of d-planes in \R^(d+c)

    sage: g = Gr(2,4) 
    sage: len(g.H) - 1 == g.top_cohomological_dim
    True
    
    # Check that the generated ideal yields a ring of the right dimension
    sage: d,c = 2,4
    sage: Gr(d,c).ideal.vector_space_dimension() == binomial(d+c,d)
    True
    
    sage: d,c = 3,6
    sage: Gr(d,c).ideal.vector_space_dimension() == binomial(d+c,d)
    True

    sage: d,c = 2,1
    sage: Gr(d,c).ideal.vector_space_dimension() == binomial(d+c,d)
    True
    
    """
    def __init__(self, d,c):
        self.d =d 
        self.c = c
        self.ws = ["w"+str(i) for i in range(1,d+1)]
        self.wbars = ["wbar" + str(i) for i in range(1,c+1)]
        self.S = PolynomialRing(FiniteField(2),",".join(self.ws)+","+",".join(self.wbars))
        
        # Constructs Cohomology Ring 
        self.R = PolynomialRing(FiniteField(2),",".join(self.ws),order = "deglex")
        self.W = [1]+[self.R.gens()[i - 1] for i in range(1,self.d + 1)]
        self.sq_dict = {}

        # print "Defining Relations..."
        self.relations = [self.R(self.wbar_expression_in_w(i)) for i in range(c+1,c+d+1)]
        # print self.relations
        # print "Forming the ideal..."
        self.ideal = self.R.ideal(self.relations)
        # print "Taking the quotient..."
        self.Mod2CohomologyRing = self.R.quotient(self.ideal)

        # print "Detereming Cohomology..."
        self.top_cohomological_dim = self.d*self.c
        self.H = [self.additive_basis_for_qth_cohomology(q) for q in range(self.top_cohomological_dim+1)]

        # print sum([len(y) for y in self.H])

    def coordinate_vector_in_qth_basis(self,word,q):
        """
        sage: g = Gr(3,4)
        sage: word = g.W[1]
        sage: g.coordinate_vector_in_qth_basis(word,1)
        [1]

        sage: word = g.W[1]^4 + g.W[2]^2
        sage: g.coordinate_vector_in_qth_basis(word,4)
        [1, 0, 0, 1]
        """
        # print word
        basis = self.additive_basis_for_qth_cohomology(q)
        if word == 0:
            return [0 for b in basis]
        return [word.monomial_coefficient(self.normal_form(b)) for b in basis]
        return 0

    def q_degree_element_from_vector(self,vector,q):
        basis = self.additive_basis_for_qth_cohomology(q)
        return sum([vector[i]*basis[i] for i in range(len(vector))])
        
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
        """
        sage: g = Gr(3,5)
        sage: g.additive_basis_for_qth_cohomology(4) == [g.W[2]^2,g.W[1]*g.W[3],g.W[1]^2*g.W[2],g.W[1]^4]
        True
        """
        S = range(0,self.c+1)
        exponents = cartesian_product([S for i in range(self.d)]).list()
        exponents = filter(lambda x: sum(x) <= self.c, exponents)
        exponents = filter(lambda x: sum([(i+1)*x[i] for i in range(0,len(x))]) == q, exponents)
        return [self.make_word_with_exponents(ex) for ex in exponents]

    def make_word_with_exponents(self,expontents):
        """
        sage: g = Gr(3,10)
        sage: g.make_word_with_exponents([1,2,4]) == g.W[1]*g.W[2]^2*g.W[3]^4
        True
    
        sage: g.make_word_with_exponents([])
        1
        """
        word = 1
        for i,e in enumerate(expontents):
            word *= self.w(i+1)^e
        return word

    def normal_form(self,expression):
        return self.ideal.reduce(self.R(expression))
        
    def wbar_expression(self,k):
        """
        sage: g = Gr(2,5)
        sage: g.wbar_expression(1)
        w1
        
        sage g.wbar_expression(2)
        w2 + w1*w2
        """
        expression = self.get_w(k)
        for i in range(1,k):
            expression += self.get_w(i)*self.get_wbar(k-i)
        # print expression
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


    def commutator(self,f,g):
        return lambda x : f(g(x)) + g(f(x))

    def Q(self,n):
        # print n
        if n == 0:
            return lambda x : self.Sq(1,x)
        # print n
        return self.commutator(lambda x : self.Sq(2^(n),x),self.Q(n-1))

    def Sq(self,i,expression):
        """
        sage: g = Gr(2,4) 

        # Squaring Condition
        sage: w1 = g.W[1]
        sage: g.Sq(1,w1) == w1^2
        True
        
        # Unstable Condition
        sage: g.Sq(2,w1) == 0
        True

        # Additive Homomorphism
        sage: w2 = g.W[2]
        sage: g.Sq(1,w1+w2) == g.Sq(1,w1)+g.Sq(1,w2)
        True

        # Cartan Formula
        sage: Sq = g.Sq
        sage: Sq(3,w1^2*w2) == Sq(0,w1^2)*Sq(3,w2)+Sq(1,w1^2)*Sq(2,w2)+Sq(2,w1^2)*Sq(1,w2)+Sq(3,w1^2)*Sq(0,w2)
        True
        
        """
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
        if len(self.R.gens()) != 1 and len(list(expression)) > 1:
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

        
        degree =  word.exponents()[0] # for univariant 
        
        if len(self.R.gens()) > 1:
            degree = sum([(x+1)*e for x,e in enumerate(word.exponents()[0])])
        

        if i > degree:
            return 0

        sub_one = {}
        for w in self.W:
            sub_one[w] = 1

        coefficient = word.subs(sub_one)

        if word in self.W:
            return self.wu_formula(i,self.W.index(word))


        x,y = self.split_off_letter(word)

        # step = "Sq("+str(i)+","+str(expression)+") = "
        # step_answer = "Sq("+str(i)+","+str(expression)+") = "

        answer = 0
        for j in range(i+1):
            answer = self.R(answer)
            # step += "Sq(" + str(i-j) +"," + str(x)+ ")Sq(" + str(j) + "," + str(y) + ")+"
            s1 = self.R(self.Sq(i-j,x))
            s2 = self.R(self.Sq(j,y))
            # step_answer += str(s1)+"*"+str(s2)+"+"
            answer += s1*s2

        # print step[:-1]
        # print step_answer[:-1]
        return answer

    def split_off_letter(self,word):
        """
        sage: g = Gr(4,5)
        sage: word = g.W[1]*g.W[2]
        sage: g.split_off_letter(word)
        (w1, w2)

        sage: word = g.W[1]^2*g.W[2]*g.W[4]
        sage: g.split_off_letter(word)
        (w1, w1*w2*w4)

        sage: word = g.W[1]
        sage: g.split_off_letter(word)
        (w1, 1)

        sage: word = 1
        sage: g.split_off_letter(word)
        (1, 1)
        """
        if word == 1:
            return (1,1)
        
        exponents = list(word.exponents()[0])
        
        first_nonzero = 0
        while exponents[first_nonzero] == 0:
            first_nonzero += 1

        letter = self.W[first_nonzero+1]
        exponents[first_nonzero] = exponents[first_nonzero] - 1
        new_word = 1

        for i in range(1,len(self.W)):
            new_word *= self.W[i]^exponents[i-1]
            
        return letter, new_word

    def wu_formula(self,i,j):
        value = 0
        # step = "WU: Sq("+str(i)+","+"w"+str(j)+") = "
        for t in range(0,i+1):
            value += binomial(j+t-i-1,t)*self.w(i-t)*self.w(j+t)
            # print t
            # print value
            # step += str(binomial(j+t-i-1,t))+"w"+str(i-t)+"*w"+str(j+t)+"+"
        # print step[:-1] + " = " + str(value)
        return value

    def element_degree(self,element):
        if element == 0:
            return 0
        return self.word_degree(list(element)[0][1])

    def word_degree(self,word):
        degree =  word.exponents()[0] # for univariant 
        if len(self.R.gens()) > 1:
            degree = sum([(x+1)*e for x,e in enumerate(word.exponents()[0])])
        return degree

    def all_elements_in_degree_n(self,n):
        all_elements = Combinations(self.H[n])
        element_list = []
        for list_of_gens in all_elements:
            element = sum(list_of_gens)
            element_list.append(element)
        return element_list
