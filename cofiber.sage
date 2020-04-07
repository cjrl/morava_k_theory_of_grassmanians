
# B = K(n,Gr(2,2^(n+1)+2*t+1)).third_page_ranks()
# print Matrix([[i for i in range(1,29)],A,B+[0,0]])

def wrap_index(arr):
     indices = range(0,len(arr))
     return Matrix([indices,arr])

def pad_list_to_length(arr,length):
    if length > len(arr):
        return arr + [0 for i in range(length - len(arr))]
    return arr

def wrap_lists_with_indices(lists):
    pad_length = max([len(x) for x in lists])
    padded = [pad_list_to_length(x,pad_length) for x in lists]
    indices = [range(0,pad_length)]
    return Matrix(indices + padded)

def new_generators(c):
     g = Gr(2,c)
     h = Gr(2,c+1)
     H = [[x for x in gens if i >= len(g.H) or x not in g.H[i]] for i,gens in enumerate(h.H)] 
     return H

def kernel_elements(c):
     g = Gr(2,c)
     H = new_generators(c)
     K = []
     for gens in H:
          K.append([x+g.normal_form(x) for x in gens])
     return K

def new_generators_d(d,c):
     g = Gr(d,c)
     h = Gr(d,c+1)
     H = [[x for x in gens if i >= len(g.H) or x not in g.H[i]] for i,gens in enumerate(h.H)] 
     return H

def kernel_elements_d(d,c):
     g = Gr(d,c)
     H = new_generators_d(d,c)
     K = []
     for gens in H:
          K.append([x+g.normal_form(x) for x in gens])
     return K

def conjectured_kernel_elements(c):
     h = Gr(2,c+1)
     K = [[] for i in range(c+1)]
     for i in range(c-1,2*c+1):
          K.append([h.normal_form(h.W[1]^(i-(c-1))*h.wbar_expression_in_w(c+1))])
     return K

def homology_of_cofiber(c,n):
    diff = 2^(n+1)-1
    g = Gr(2,c)
    h = Gr(2,c+1)
    Z = [x[0] if len(x) == 1 else 0 for x in kernel_elements(c) ]
    im = [h.normal_form(g.Q(n)(x[0])) if len(x) == 1 else 0 for x in kernel_elements(c) ]
    for i,x in enumerate(im):
        if i+diff < len(im) and im[i] == Z[i+diff]:
             Z[i] = 0
             Z[i+diff] = 0
    return Z



# n = 2
# i = 5
# k = 2^(n+1)+2*i+1
# d = 2^(n+1)+1
# # n = 1
# # k = 2^(n+1)+2
# m = k
# diff = 2^(n+1)-1
# c = k-2
# # print "c =",k-2
# # print "skipped =",k+2^(n+1)+8

# h = Gr(2,c+1)
# g = Gr(2,c)

def preimage_under_Qn(n,element,g):
     diff = 2^(n+1)-1
     # if element == 0:
     #      return 0
     source_degree = g.element_degree(element) - diff
     if source_degree <= -1:
          return []
     preimage = [x for x in elements_in_g_H_n(g,source_degree) if g.normal_form(g.Q(n)(x)) == g.normal_form(element)]
     # return elements_in_g_H_n(g,source_degree)
     return preimage

def elements_in_g_H_n(g,n):
     gens = g.H[n] 
     all_elements = Combinations(gens)
     element_list = []
     for list_of_gens in all_elements:
          element = sum(list_of_gens)
          element_list.append(element)
     return element_list

def image_under_Qn(g,n):
     return [g.normal_form(g.Q(n)(x)) for x in elements_in_g_H_n(g,n)]

# m into m+1
# this is the degree in which the lowest new generator in m appears
# c + 1 
def cofiber_lowest_degree(m):
     c = m - 2
     return c + 1

def cofiber_heighest_degree(m):
     return 2*cofiber_lowest_degree(m)

def preimage_degree_for_Qn(degree,n):
     diff = 2^(n+1)-1
     return degree - diff

# def alpha(m,n):
#      degree = preimage_degree_for_Qn(cofiber_lowest_degree(m),n)+1
#      g = Gr(2,m-2)
#      return sum(g.H[degree])

# def beta(m,n):
#      degree = preimage_degree_for_Qn(cofiber_lowest_degree(m),n)+3
#      g = Gr(2,m-2)
#      return sum(g.H[degree])


# def alpha_beta_test(i,n):
#      m = 2^(n+1)+i
#      c = m-2
#      h = Gr(2,c+1)
#      g = Gr(2,c)

#      print [preimage_under_Qn(n,h.W[1]*x,h) for x in terms_in_wbar_formula(c+1,h)]
     # first_gen = kernel_elements(c)[cofiber_lowest_degree(m)+1][0]
     # alpha_target = h.normal_form(h.Q(n)(alpha(m,n)))
     # alpha_outcome = (first_gen == alpha_target)
     # print first_gen,h.normal_form(h.Q(n)(alpha(m,n))) 
     # print "alpha: ",alpha_outcome
     # print alpha_target
     # print first_gen 

     # first_gen = homology_of_cofiber(c,n)[cofiber_lowest_degree(m)+1]
     # # print homology_of_cofiber(c,n)
     # print first_gen
     # alpha_target = h.normal_form(h.Q(n)(alpha(m,n)))
     # # alpha_outcome = h.R(beta(m,n)) in [g.normal_form(x) for x in preimage_under_Qn(n,beta_target,h)]
     # # print h.normal_form(h.Q(n)(alpha(m,n)))
     # alpha_outcome = alpha_target == first_gen 
     # print "alpha: ", alpha_outcome
     # print [g.normal_form(x) for x in preimage_under_Qn(n,first_gen,h)]


     # print g.H[preimage_degree_for_Qn(cofiber_lowest_degree(m)+1,n)]

     # second_gen = homology_of_cofiber(m,n)[cofiber_lowest_degree(m)+3]
     # beta_target = h.normal_form(h.Q(n)(beta(m,n)))
     # # beta_outcome = h.R(beta(m,n)) in [g.normal_form(x) for x in preimage_under_Qn(n,beta_target,h)]
     # # print h.normal_form(h.Q(n)(beta(m,n)))
     # beta_outcome = (beta_target == second_gen)
     # print "beta: ", beta_outcome
     # print [g.normal_form(x) for x in preimage_under_Qn(n,second_gen,h)]
     # # print h.W[1]^preimage_degree_for_Qn(cofiber_lowest_degree(m)+3,n)
     # # print g.normal_form(h.W[1]^preimage_degree_for_Qn(cofiber_lowest_degree(m)+3,n))
     # # print [g.normal_form(x) for x in preimage_under_Qn(n,second_gen,h)]
     # # print g.H[preimage_degree_for_Qn(cofiber_lowest_degree(m)+3,n)]
     # # print beta(m,n)
     # # print homology_of_cofiber(m,n)[cofiber_lowest_degree(m)+3]
     # return alpha_outcome and beta_outcome


# def multiple_alpha_beta_tests():
#     for i in range(1,10):
#         for n in range(1,5):
#             print 2*i,n
#             outcome = alpha_beta_test(2*i,n)
#             print outcome
#             # if not outcome:
#             #      return False

# def preimages():
#     for i in range(1,10):
#         for n in range(1,5):
#             print 2*i,n
#             m = 2^(n+1)+2*i
#             c = m - 2
#             h = Gr(2,c+1)
#             print h.Q(n)(h.wbar_expression_in_w(c+1))==h.W[1]^(2^(n+1)-1)*h.wbar_expression_in_w(c+1)


def conj_gen(a,b,h):
     return h.W[1]^(a-1)*h.W[2]^(b-2^n+1)

def difference():
    preimage = preimage_under_Qn(n,h.W[1]*h.wbar_expression_in_w(c+1),h)
    print preimage
    return [g.normal_form(x) - h.wbar_expression_in_w(h.element_degree(preimage[0])) for x in preimage]

def recursive_wbar_test(a):
     h = Gr(2,a+2)
     return h.wbar_expression_in_w(a) == h.W[1]*h.wbar_expression_in_w(a-1)+ h.W[2]*h.wbar_expression_in_w(a-2)

def wbar_formula_conjecture(n,l):
     m = 2^(n+1)+2*l
     d = 2^(n+1)+1
     c = m - 2
     h = Gr(2,c+1)
     goal = h.normal_form(h.W[1]*h.wbar_expression_in_w(c+1))
     # conj = h.W[1]^(2*l-1)*h.wbar_expression_in_w(d)
     # conj += h.W[1]^(2*l-3+d-2)*h.W[2]^2+h.W[1]^(d-1)*h.W[2]^l

     conj = h.W[1]^(2*l-2)*h.wbar_expression_in_w(d)
     conj += h.W[1]^(2*l-4)*h.W[2]^2*h.wbar_expression_in_w(d-2)
     conj += h.W[2]^l*h.wbar_expression_in_w(d-2)
     conj *= h.W[1]
     print goal
     print h.normal_form(conj)


def terms_in_wbar_formula(l,h):
     pairs = []
     for i in range(0,l+1):
          for j in range(0,l+1):
               if i + 2*j == l:
                    if binomial(i+j,i)%2 == 1:
                        pairs.append(h.W[1]^i*h.W[2]^j)
     return pairs

def all_d(n,l):
     d = 2^(n+1)+1
     m = 2^(n+1)+2*l
     h = Gr(2,m-2)
     pairs = []
     for i in range(0,l+1):
          for j in range(0,l+1):
               if i + 2*j == 2*l-2:
                        pairs.append(h.W[1]^i*h.W[2]^j*h.wbar_expression_in_w(d))
     # return pairs
     return [preimage_under_Qn(n,h.W[1]*x,h) for x in pairs]

# find . -name "*.sage" -print | etags -l "python" -
# ctags --language-force=python --python-types=+l -e "cofiber.sage" 

# def commuting_sq_test(c):
#      g = Gr(2,c)
#      h = Gr(2,c+1)

#      n = 1

#      for gens in h.H:
#           for x in gens:
#                print g.normal_form(h.Q(n)(x)) == g.Q(n)(g.normal_form(x))


# for i in range(0,10):
#      m = 2^(n+1)+2*i +1
#      c = m - 2
#      g = Gr(2,m-2)
#      h = Gr(2,m+1-2)
#      print g.normal_form(h.W[1]*h.wbar_expression_in_w(c+2))

def groeb_basis_conj(c):
    h = Gr(2,c+1)
    g = Gr(2,c)

    conjecture = []
    for i in range(0,(c+1)+1):
         conjecture.append(h.normal_form(h.W[1]^i*h.wbar_expression_in_w(c+1)))
    print g.ideal.groebner_basis()
    print conjecture
    return g.ideal.groebner_basis() == conjecture

     

# this c is for m, h is m+1
def identify_as_groeb_element(c,element,h):
    parent_g = Gr(2,c+2) 
    for i in range(0,(c+2)+1):
        groeb = h.W[1]^i*h.wbar_expression_in_w(c+2)
        groeb = parent_g.normal_form(groeb)
        # print element
        if element == groeb:
            return "w1^"+str(i)+"wbar(c+2)"
    return None

def expand_in_terms_of_groebner(poly,ideal):
    gb = ideal.groebner_basis()
    gb = list(gb)
    steps = [gb[:i] for i in range(len(gb)+1)]
    expansion = []
    remainder = 0
    for step in steps:
        # print "poly:", poly.reduce(step)
        if len(step) > 0:
            quotient = (poly - poly.reduce(step)) / step[-1]
            # print "quotient:", quotient
            # print poly,"=(",quotient,")*(",step[-1],")","+(",poly.reduce(step),")"
            expansion.append([quotient,step[-1]])
            remainder = poly.reduce(step)
            # print step[-1]
        poly = poly.reduce(step)
    return expansion,remainder

n = 1 
i = 1

m = 2^(n+1)+2*i
d = 2^(n+1)+1
c = m-2
h = Gr(2,c+1)
g = Gr(2,c)
k = K(n,g)

gb = h.ideal.groebner_basis()

print groeb_basis_conj(c)

print "Possible Preimage Generators: "
degree = preimage_degree_for_Qn(c+1,n)
possible_pairs =[x for x in h.H[degree]]
print possible_pairs  
print len(possible_pairs)

print "\n\n a odd b even:"

# degree = preimage_degree_for_Qn(2*(c+1),n)-(2^(n+1)+1)

# a = 1
# b = (degree - a)/2

# print "a=",a
# print "b=",b
# expansion = expand_in_terms_of_groebner(h.W[1]^(a+d-2)*h.W[2]^b,h.ideal)

# for i,line in enumerate(expansion[0]):
#      # print line
#      # print line[0],"w1^"+str(i)+"wbar(c+2)" 
#      print line[0],"g"+str(i)
#      # print identify_as_groeb_element(c,line[1],h)
# print expansion[1]

# for a and b odd

# a = 0
# b = 1

# print a
# print b

# print "w2^(b-1)*(wbar(d)+w1^d) expansion:"
# expansion = expand_in_terms_of_groebner(h.W[2]^(b-1)*(h.wbar_expression_in_w(d)+h.W[1]^d),h.ideal)

# for i,line in enumerate(expansion[0]):
#      # print line
#      # print line[0],"w1^"+str(i)+"wbar(c+2)" 
#      print line[0],"g"+str(i)
#      # print identify_as_groeb_element(c,line[1],h)
# print expansion[1]

# expansion = expand_in_terms_of_groebner(h.Q(n)(h.W[1]^a*h.W[2]^b),h.ideal)

# # for i,line in enumerate(expansion[0]):
# #      # print line
# #      print line[0],"w1^"+str(i)+"wbar(c+2)" 
# #      # print identify_as_groeb_element(c,line[1],h)
# # print expansion[1]

# new_list = []
# for i in range(0,len(possible_pairs)):
#     print ""
#     print i
#     terms = expand_in_terms_of_groebner(h.W[2]^(b-1-i)*(h.wbar_expression_in_w(d)+h.W[1]^d),h.ideal)
#     print (a,b)
#     new_line = []
#     for line in terms[0]:
#         # print line 
#         line.pop()
#         new_line.append(line)
#     new_line = [1 if x != 0 else 0 for x in flatten(new_line)]
#     new_list.append(new_line)

# def pretty_print(data):
#     for i,line in enumerate(data):
#         line = "".join([str(x) for x in line])
#         line = str.replace(line,"0"," ")
#         print line


# pretty_print(new_list)

# this looks good
def g4_conj(n,l):
     m = 2^(n+1)+2*l+1
     d = 2^(n+1)+1
     c = m-2
     h = Gr(2,c+1)
     gb = h.ideal.groebner_basis()

     print h.W[2]^(c+1-2^(n+1)+2)*(h.wbar_expression_in_w(d)+h.W[1]^d)
     print gb[4+2*l]

# for i,gens in enumerate(k.detailed_homology()):
#      if i == c+1-2^(n+1)+1:
#           print "c+1-diff",i,gens,[h.normal_form(h.Q(n)(x)) for x in gens]
#      else:
#           print i,gens,[h.normal_form(h.Q(n)(x)) for x in gens]

alpha = lambda c,n,g: g.W[1]*g.wbar_expression_in_w(c+1-2^(n+1)+1)
beta = lambda i,g: g.W[2]^(2*i+1)

def test_conj():
     for i in range (1,5):
          for n in range(1,5):
               m = 2^(n+1)+2*i
               d = 2^(n+1)+1
               c = m-2
               h = Gr(2,c+1)
               g = Gr(2,c)

               print i,n
               # print h.Q(n)(h.W[1]*h.wbar_expression_in_w(c+1-2^(n+1)+1))
               print h.Q(n)(alpha(c,n,h)) == h.W[1]*h.wbar_expression_in_w(c+1)
               print h.normal_form(h.Q(n)(beta(i,h))) == h.normal_form(h.W[1]^(2*(i+1))*h.wbar_expression_in_w(c+1))
               
               # print h.Q(n)(g.W[1]*g.wbar_expression_in_w(c+1-2^(n+1)+1)) == h.W[1]^6*h.wbar_expression_in_w(c+1)


test_conj()            
