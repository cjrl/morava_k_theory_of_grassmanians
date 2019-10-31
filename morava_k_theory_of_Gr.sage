load("gr.sage")
load("K.sage")

def compute_morava_k_theory_of_Gr(n,d,m):
    c = m - d 

    if d >= m:
        dim = 0
        f = open("/Users/Chris/box/research/general_grassmannian/gr_data/Gr("+
             str(d) +","+ str(m) + ")_"+str(n),"w+")
        f.write("return "+str(dim)+";")       
        return 0
    
    gr = Gr(d,c)
    dim = K(n,gr).rank()
    print dim
    f = open("/Users/Chris/box/research/general_grassmannian/gr_data/Gr("+
             str(d) +","+ str(m) + ")_"+str(n),"w+")
    f.write("return "+str(dim)+";")
    
