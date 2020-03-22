load("K.sage")

def show_groebner_steps(poly,ideal):
    gb = ideal.groebner_basis()
    gb = list(gb)
    steps = [gb[:i] for i in range(len(gb))]
    for step in steps:
        print step
        print poly.reduce(step)

# Check that kernel is working propertly 

g = Gr(2,4)
k = K(1,g)
for i,gens in enumerate(g.H):
    print "N",i
    # M = k.qth_diffeential_as_matrix(i)
    # print M
    for gen in gens:
        print gen,"->",k.unreduced_d(gen)
        print gen,"->",k.d(gen)
        
    # print "\nAll elements of H^"+str(i)
    # all_elements = Combinations(gens)
    # for list_of_gens in all_elements:
    #     element = sum(list_of_gens)
    #     print element,"->",k.unreduced_d(element)
    #     print element,"->",k.d(element)

    # print k.qth_differential_as_matrix(i)
    A = k.qth_differential_as_matrix(i).echelon_form()
    print "rank:",A.rank()
    if i % 2 == 1:
        print ceil((i+1)/2)- floor((i-1)/4)
        print floor((i-1)/4)
    print "nullity:", A.right_nullity()
    print A
    B = k.qth_differential_as_matrix(g.d*g.c -i - 2^(k.n+1) + 1).T.echelon_form() 
    # print A
    # print ""
    # print B
        # if k.unreduced_d(gen) != 0:
        #     show_groebner_steps(k.unreduced_d(gen),g.ideal)
        # print M*Matrix(g.coordinate_vector_in_qth_basis(gen,i)).T
        # print vector(Matrix(g.coordinate_vector_in_qth_basis(gen,i)).T) in M.right_kernel()
    # print k.qth_differential_as_matrix(i).left_kernel()
    print "\n"

# c = 10
# g = Gr(2,32)
# k = K(2,g)
# g2 = Gr(2,c+2)
# k2 = K(2,g2)
# for i,gens in enumerate(g.H):
#     for gen in gens:
#         print gen,",",k.d(gen),",",k2.d(g2.R(gen))

# print "new gens"
# new_gens = [x for x in flatten(g2.H) if x not in flatten(g.H)]
# for gen in new_gens:
#     print gen,",",k2.d(gen)

# for c in range(30,40):
#     g = Gr(2,1+c)
5#5 5 5   k = K(2,g)
#     # print k.rank()
#     # print sum(k.second_page_ranks())
#     print k.rank()

# M = Matrix(QQ,[[1,-1]])

# https://oeis.org/search?q=2+5+9+14+20+27+35+44+54+65+77+90+104+119+135+152+170+189+209+230&language=english&go=Search

# rows = []
# for c in range(0,2):
    

# for i in range(0,15):
#     print Matrix(K(1,Gr(2,10+2*i)).third_page_ranks()) 
    
# def test_conjecture(c):
#     g = Gr(2,c)
#     ranks = K(1,g).third_page_ranks()
#     for j,rank in enumerate(ranks):
#         i = j + 1
#         if rank == 1 and Mod(i,4) != 0 and i != 2 and i != c-1 and i != c + 1 and i != 2*c - 2:
#             print i
#             return False
#     return True
    
