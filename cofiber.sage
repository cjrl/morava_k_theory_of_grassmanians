
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


def P_k(k):
    H = [1 for x in range(0,2*k)]
    for i in range(0,k): 
        H[i] = 0
    return H 

def Qn_homology_P_k(n,k):
    # P_k = RP^(2k-1)/RP^k-1
    diff_degree = 2^(n+1)-1
    H = P_k(k) 
    for i,x in enumerate(H):
        if i % 2 == 1:
             if H[i] == 1 and i + diff_degree < len(H):
                H[i] = 0
                H[i+diff_degree] = 0
    return H[1:]

# print wrap_index(Qn_homology_P_k(2,13))

# n = 2
# k = 2^(n+1)+3

# print "k =",k
# print "skipped =",k+2^(n+1)

# A = K(n,Gr(2,k+1)).third_page_ranks()
# B = K(n,Gr(2,k)).third_page_ranks()
# shifted = [0,0,0]+Qn_homology_P_k(n,k)
# print wrap_lists_with_indices([A,shifted,B])

# n = 1
# k = 2^(n+1)+2

# print "k =",k
# print "skipped =",k+2^(n+1)

# A = K(n,Gr(2,k+2)).third_page_ranks()
# B = [0,0]+K(n,Gr(2,k)).third_page_ranks()
# shifted = [0,0]+Qn_homology_P_k(n,k)
# print wrap_lists_with_indices([A,B])

# n = 2
# k = 2^(n+1)

# print "k =",k
# print "skipped =",k+2^(n+1)

# A = K(n,Gr(2,k+1)).third_page_ranks()
# B = K(n,Gr(2,k)).third_page_ranks()
# shifted = Qn_homology_P_k(n,k)
# print wrap_lists_with_indices([A,shifted,B,Qn_homology_P_k(n,k)[3:]])
# print sum(A)
# print sum(B)
# print sum(shifted)

# n = 2
# k = 2^(n+1)+3

# print "k =",k
# print "skipped =",k+2^(n+1)

# A = K(n,Gr(2,k+1)).third_page_ranks()
# B = K(n,Gr(2,k)).third_page_ranks()
# shifted = [0,0,0]+Qn_homology_P_k(n,k)
# print wrap_lists_with_indices([A,B,Qn_homology_P_k(n,k)])

# n = 1
# k = 2^(n+1)+2

# print "k =",k
# print "skipped =",k+2^(n+1)+4

# g = Gr(2,k+1-2)
# A = K(n,Gr(2,k+1-2)).third_page_ranks()
# B = K(n,Gr(2,k-2)).third_page_ranks()
# # shifted = Qn_homology_P_k(n,k)[2^(n+1)-2:]
# shifted = Qn_homology_P_k(n,k)
# # shifted = Qn_homology_P_k(n,k)[-2:]
# # shifted = [0 for i in range(2^(n+1)-1)]+Qn_homology_P_k(n,k)
# # shifted2 = shifted + [0 for i in range(20)]
# C = [A[i]+shifted2[i] for i in range(0,len(A))]
# print wrap_lists_with_indices([shifted,A,B])

# print sum(A)
# print sum(B)
# # print sum([x for i,x in enumerate(shifted) if B[i]!=0])
# print(sum(shifted))



# g = Gr(2,k+1)
# h = Gr(2,k)
# index = 2*(k+1)-2^(n+1)+1
# print index
# print h.H[index]
# classes = g.H[index]
# print classes
# print [g.normal_form(g.Q(n)(x)) for x in classes]
# x = h.H[index][-1]
# print x


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

def conjectured_kernel_elements(c):
     h = Gr(2,c+1)
     K = [[] for i in range(c+1)]
     for i in range(c-1,2*c+1):
          K.append([h.normal_form(h.W[1]^(i-(c-1))*h.wbar_expression_in_w(c+1))])
     return K

# print new_generators(6-2)

# c = 5
# g = Gr(2,c)
# h = Gr(2,c+1)

# print kernel_elements(c)
# print conjectured_kernel_elements(c)
# print [h.normal_form(g.Q(1)(x[0])) if len(x) == 1 else 0 for x in kernel_elements(c) ]

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


# print [1 if x != 0 else x for x in homology_of_cofiber(10,3)]

n = 2
i = 5
k = 2^(n+1)+2*i+1
d = 2^(n+1)+1
# n = 1
# k = 2^(n+1)+2
m = k
diff = 2^(n+1)-1
c = k-2
# print "c =",k-2
# print "skipped =",k+2^(n+1)+8

h = Gr(2,c+1)
g = Gr(2,c)

# A = K(n,Gr(2,k+1-2)).third_page_ranks()
# B = K(n,Gr(2,k-2)).third_page_ranks()
# shifted = [1 if x != 0 else x for x in homology_of_cofiber(k-2,n)][2^(n+1)-1:] # odd shift
# shifted = [1 if x != 0 else x for x in homology_of_cofiber(k-2,n)] # even shift
# print wrap_lists_with_indices([shifted,A,B])

# print "C("+str(k)+") =", sum(shifted)
# print "Gr("+str(k+1)+") =", sum(A)
# print "Gr("+str(k)+") =", sum(B)

# print homology_of_cofiber(k-2,n)

# print (new_generators(k-2))
# print [h.normal_form(h.Q(n)(x)) for x in homology_of_cofiber(k-2,n)]
# index = 2*(c+1)-2*diff
 # index = 2
# print g.H[index]
# print new_generators(k-2)
# print [h.normal_form(h.Q(n)(x)) for x in g.H[index]]
# print sum([h.normal_form(h.Q(n)(x)) for x in g.H[index]])
# # print homology_of_cofiber(k-2,n)

# print kernel_elements(c) == conjectured_kernel_elements(c)

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

def alpha(m,n):
     degree = preimage_degree_for_Qn(cofiber_lowest_degree(m),n)+1
     g = Gr(2,m-2)
     return sum(g.H[degree])

def beta(m,n):
     degree = preimage_degree_for_Qn(cofiber_lowest_degree(m),n)+3
     g = Gr(2,m-2)
     return sum(g.H[degree])

# print homology_of_cofiber(k-2,n)
# # print [h.element_degree(x) for x in homology_of_cofiber(k-2,n)]
# preimage = [preimage_under_Qn(n,x,h) for x in homology_of_cofiber(k-2,n)]
# for i,x in enumerate(preimage):
#      if x != 0:
#           print x
#           H = h.H[len(h.H)-len(preimage)-diff+i]
#           print H
#           # print g.H[len(h.H)-len(preimage)-diff+i]
#           print [y for y in H if y in x]
#           # print [g.normal_form(y) for y in H if y in x]
#           # print H[-1]
#           # print [h.normal_form(h.Q(n)(y)) for y in x]
#           print homology_of_cofiber(k-2,n)[i]
#           # con_preimage = h.H[len(h.H)-len(preimage)-diff+i][-1] 
#           # print "conj preimage", "H: " ,con_preimage
#           print ""
#           # diff = 2^(n+1)-1
#         #   print x
#         #   # print "H: " ,h.H[len(h.H)-len(preimage)-diff+i] 
#         #   # print "H: " ,[g.normal_form(x) for x in h.H[len(h.H)-len(preimage)-diff+i] ]
#         # print homology_of_cofiber(k-2,n)[i]
#         # # print "preimage: ", x, "->", h.normal_form(h.Q(n)(x))# , homology_of_cofiber(k-2,n)[i]
#         # # print "conj preimage", "H: " ,h.H[len(h.H)-len(preimage)-diff+i][-1] 
#         # # print "preimage: ", g.normal_form(x)

# cofib = homology_of_cofiber(k-2,n)
# x = cofib[13]
# print x
# preimage = preimage_under_Qn(n,x,h)
# print preimage
# print [h.normal_form(h.Q(n)(x)) for x in preimage]


def alpha_beta_test(i,n):
     m = 2^(n+1)+i
     c = m-2
     h = Gr(2,c+1)
     g = Gr(2,c)

     print [preimage_under_Qn(n,h.W[1]*x,h) for x in terms_in_wbar_formula(c+1,h)]
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


def multiple_alpha_beta_tests():
    for i in range(1,10):
        for n in range(1,5):
            print 2*i,n
            outcome = alpha_beta_test(2*i,n)
            print outcome
            # if not outcome:
            #      return False

def preimages():
    for i in range(1,10):
        for n in range(1,5):
            print 2*i,n
            m = 2^(n+1)+2*i
            c = m - 2
            h = Gr(2,c+1)
            print h.Q(n)(h.wbar_expression_in_w(c+1))==h.W[1]^(2^(n+1)-1)*h.wbar_expression_in_w(c+1)


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

def commuting_sq_test(c):
     g = Gr(2,c)
     h = Gr(2,c+1)

     n = 1

     for gens in h.H:
          for x in gens:
               print g.normal_form(h.Q(n)(x)) == g.Q(n)(g.normal_form(x))


# for i in range(0,10):
#      m = 2^(n+1)+2*i +1
#      c = m - 2
#      g = Gr(2,m-2)
#      h = Gr(2,m+1-2)
#      print g.normal_form(h.W[1]*h.wbar_expression_in_w(c+2))

# this c is for m, h is m+1
def identify_as_groeb_element(c,element,h):
    for i in range(1,(c+2)+1):
        groeb = h.W[1]^i*h.wbar_expression_in_w(c+2)
        groeb = Gr(2,c+3).normal_form(groeb)
        print element
        if element == groeb:
            return "w1^"+str(i)+"wbar(c+2)"
    return None
