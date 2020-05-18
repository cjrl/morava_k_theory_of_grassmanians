
def qth_matrix_of_inclusion(q,gr_m,gr_m_plus_one):
    q_basis = gr_m_plus_one.additive_basis_for_qth_cohomology(q)
    image = map(gr_m.normal_form,q_basis)
    columns = []
    for x in image:
        columns.append(gr_m.coordinate_vector_in_qth_basis(x,q))
    return Matrix(FiniteField(2),columns).T



# cofiber = [[] for x in range(0,d*(c+1)+1)]
# for q in range(d*(c+1)+1):
#     basis_of_kern = qth_matrix_of_inclusion(q,gr_m,gr_m_plus_one).right_kernel().basis()
#     for b in basis_of_kern:
#         # print gr_m_plus_one.q_degree_element_from_vector(b,q)
#         cofiber[q].append(gr_m_plus_one.q_degree_element_from_vector(b,q))
#     # print ""
# print cofiber

# print gr_m.ideal.groebner_basis()

def get_cofiber(gr,gr_m_plus_one):
    d = gr.d
    c = gr.c
    cofiber = [[] for x in range(0,d*(c+1)+1)]
    for q in range(d*(c+1)+1):
        basis_of_kern = qth_matrix_of_inclusion(q,gr_m,gr_m_plus_one).right_kernel().basis()
        for b in basis_of_kern:
            # print gr_m_plus_one.q_degree_element_from_vector(b,q)
            cofiber[q].append(gr_m_plus_one.q_degree_element_from_vector(b,q))
        # print ""
    return cofiber

def get_margolis_of_cofiber(gr,n):
    gr_m_plus_one = Gr(gr.d,gr.c+1)
    cofiber = get_cofiber(gr,gr_m_plus_one)
    margolis_H = []
    for q in range(gr.d*(gr.c+1)+1):
        Z = qth_diff_as_matrix(cofiber,gr_m_plus_one,q,n).right_kernel()
        image = []
        if q-(2^(n+1)-1) >= 0:
            image = map(gr_m_plus_one.Q(n),cofiber[q-(2^(n+1)-1)])
            image = map(gr_m_plus_one.normal_form,image)
            # print image
            image = map(lambda b : write_in_cofiber_basis(b,q,cofiber,gr_m_plus_one),image)
            # print image
            # image = map(Z.coordinate_vector,image)
            # image = map(vector,image)
            # print Z
            # print image
            # print map(gr_m_plus_one.Q(n),image)
        B = Z.subspace(image)
        Q = Z / B
        lifted_basis = map(Q.lift,Q.basis())
        margolis_H.append(map(lambda b : vector_to_element(b,q,cofiber),lifted_basis))
    return margolis_H
        
def write_in_cofiber_basis(element,q,cofiber,gr_m_plus_one):
    columns = []
    for b in cofiber[q]:
        # print b
        columns.append(gr_m_plus_one.coordinate_vector_in_qth_basis(b,q))
    v = Matrix(gr_m_plus_one.coordinate_vector_in_qth_basis(element,q)).T
        
    M = Matrix(columns).T
    M = M.augment(v,subdivide=True)
    # print M,"\n"
    M = M.echelon_form()
    # print M
    # print"new",list(M.columns()[-1])
    return list(M.columns()[-1][:len(cofiber[q])])

    # bound = len(cofiber[q])
    # coefficients = []
    # for i in range(0,2^bound):
    #     binary = map(int,list(bin(i)[2:]))
    #     while len(binary)<bound:
    #         binary = [0]+binary
    #     coefficients.append(binary)
    # for coeff in coefficients:
    #     in_basis = sum([coeff[i]*x for i,x in enumerate(cofiber[q])])
    #     if element == in_basis:
    #         print "old",coeff
    #         return coeff
    # return None

def vector_to_element(vec,q,cofiber):
    return sum([vec[i]*x for i,x in enumerate(cofiber[q])])
    

def qth_diff_as_matrix(cofiber,gr_m_plus_one,q,n):
    q_basis = gr_m_plus_one.additive_basis_for_qth_cohomology(q)
    # print cofiber[q]
    image = map(gr_m_plus_one.Q(n),cofiber[q])
    image = map(gr_m_plus_one.normal_form,image)
    # print "c",cofiber[q]
    # print "im",image
    # vec = map(lambda x : gr_m_plus_one.coordinate_vector_in_qth_basis(x,q+2^(n+1)-1),image)
    # print vec
    # print map(lambda z : gr_m_plus_one.q_degree_element_from_vector(z,q+2^(n+1)-1),vec) 
    columns = []
    for x in image:
        columns.append(gr_m_plus_one.coordinate_vector_in_qth_basis(x,q+2^(n+1)-1))
    return Matrix(FiniteField(2),columns).T

def groeb_cofiber(gr_m):
    groeb = gr_m.ideal.groebner_basis()
    cofiber = [[] for i in range(gr_m.d*(gr_m.c+1)+1)]
    for element in groeb:
        cofiber[gr_m.element_degree(element)].append(element)
    return cofiber

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
# d = 5
# n = 2
# l = 1
# l_copy = l

# m = 2^(n+1)+2*l
# c = m - d

# gr_m = Gr(d,c)
# # gr_m_plus_one = Gr(d,c+1)
# # cofiber = get_margolis_of_cofiber(gr_m,n)
# cofiber_margolis = map(len,get_margolis_of_cofiber(gr_m,n))


# K_gr_m = K(n,gr_m).third_page_ranks()
# K_gr_m_plus_one = K(n,gr_m_plus_one).third_page_ranks()

# print "C(m) Gr(d,m) Gr(d,m+1)"
# print wrap_lists_with_indices([map(len,get_cofiber(gr_m,gr_m_plus_one)),map(len,gr_m_plus_one.H),map(len,gr_m.H)])
# print "\nMargolis of C(m) Gr(d,m) Gr(d,m+1)"
# print wrap_lists_with_indices([cofiber_margolis,K_gr_m_plus_one,K_gr_m])

# print(sum(K_gr_m))
# print (sum(K_gr_m_plus_one))
# print sum(cofiber_margolis)
# print (sum(K_gr_m_plus_one)) + sum(cofiber_margolis)

# print wrap_lists_with_indices([map(len,get_cofiber(gr_m,gr_m_plus_one)),cofiber_margolis,[0 for i in range(8)]+K(n,Gr(d-1,c+1)).third_page_ranks()])

# print "cofiber"
# print sum(cofiber_margolis)
# print "K"
# print K(n,Gr(d-1,c+1)).rank()

def elements_in_g_H_n(g,n):
     gens = g.H[n] 
     all_elements = Combinations(gens)
     element_list = []
     for list_of_gens in all_elements:
          element = sum(list_of_gens)
          element_list.append(element)
     return element_list

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

def is_zero_map(n,q,gr_m_plus_one,cofiber):
    degree = q+2^(n+1)-1
    if degree >= len(cofiber):
        return True
    print degree
    return []==flatten(map(lambda e: preimage_under_Qn(n,e,gr_m_plus_one),cofiber[degree]))

def compute_diff(d,n,l):
    m = 2^(n+1)+2*l
    c = m - d

    gr_m = Gr(d,c)
    cofiber_margolis = map(len,get_margolis_of_cofiber(gr_m,n))
    print "Cofiber" , sum(cofiber_margolis)
    print "Gr_{d-1}" , K(n,Gr(d-1,c+1)).rank()
