load("gr.sage")
load("K.sage")


def gap_Gr(d,m):
    return Gr(d,m-d)

def compute_morava_k_theory_of_Gr(n,d,m):
    c = m - d 

    # degree/
    if d >= m:
        dim = 0
        f = open("/Users/Chris/box/research/general_grassmannian/gr_data/Gr("+
             str(d) +","+ str(m) + ")_"+str(n),"w+")
        f.write("return "+str(dim)+";")       
        f.close()
        g = open("/Users/Chris/box/research/general_grassmannian/gr_data_min/Gr("+
             str(d) +","+ str(m) + ")_"+str(n),"w+")
        g.write("return "+str(dim)+";")
        g.close()
        return dim
    
    gr = Gr(d,c)
    dim = K(n,gr).rank()
    print dim
    f = open("/Users/Chris/box/research/general_grassmannian/gr_data/Gr("+
             str(d) +","+ str(m) + ")_"+str(n),"w+")
    f.write("return "+str(dim)+";")
    f.close()

    print "Calculating Mins"
    g = open("/Users/Chris/box/research/general_grassmannian/gr_data_min/Gr("+
             str(d) +","+ str(m) + ")_"+str(n),"w+")
    g.write("return "+str(K(n,gr).minimum_possible_rank())+";")
    g.close()
    
    return dim
    

N = 10
d = 3
n = 1


def plot_at_n_for_d_c(n,d,N):
    data = [(c,compute_morava_k_theory_of_Gr(n,d,d+c)) for c in range(0,N)]
    color = (randrange(1,256)/256,randrange(1,256)/256,randrange(1,256)/256)
    return list_plot(data,plotjoined=True,legend_label = str(n),color = color,axes_labels=["$c$","rank"],title = "$Gr_"+str(d)+"(\mathbb{R}^{"+str(d)+"+"+str(c)+"})$")+list_plot(data,size = 25,color = color)

def plot_n_range(n1,n2,d,N):
    return sum(plot_at_n_for_d_c(n,d,N) for n in range(n1,n2+1))

def my_list_plot(data,name):
    color = (randrange(1,256)/256,randrange(1,256)/256,randrange(1,256)/256)
    return list_plot(data,plotjoined=False,color = color,legend_label=name)+list_plot(data,size = 25,color = color)
    

# show(list_plot(data,plotjoined=True))

# A = plot_n_range(1,5,1,25)
# B = plot_n_range(1,5,2,25)
# C = plot_n_range(1,5,3,25)
# D = plot_n_range(1,5,4,25)
 
