LoadPackage("ctbllib");

Gr := function(d,m) 
    return Concatenation(["Gr","(",String(d),",",String(m),")"]);
end;;

dim_eigenspace := function(G,chi,lambda)
    local t,normalizing_coeff;
    normalizing_coeff := 1;
    # this is a part of the issue.
    if Indicator(CharacterTable(G),[lambda],2) = [0] then
        normalizing_coeff := 2;
    fi; 
    t := CharacterTable(G);
    return (Sum( [ 1 .. Length(ConjugacyClasses(t)) ], i -> SizesConjugacyClasses( t )[i] * chi[i]
                            * ComplexConjugate( lambda[i] ) ) / Size( t )) / normalizing_coeff;
end;;

# norm := function(G,chi)
#     local sum;
#     sum := 0;
#     for i in [1..Length(chi)] do
#         sum := sum + chi[i]*ComplexConjugate(chi[i]);
#     od;
#     return sum / Order(G);
# end;;

get_indeterminate := function(name,R)
    local x;
    for x in IndeterminatesOfPolynomialRing(R) do
        if String(x) = name then
            return x;
        fi;
    od;
end;;

get_G_fixed_point_formula_for_Gr_d := function(d,chi,irr_list,G)
    local R, indeterminates, i, tau, dim, list, comb, expression, gr_list;
    indeterminates := [];
    
    gr_list := [];
    
    for tau in irr_list do
        dim := tau[1];
        for i in [dim..d] do
            # Display(dim_eigenspace(G,chi,tau));
            Append(indeterminates,[Gr(i,dim*dim_eigenspace(G,chi,tau))]);
            Append(gr_list,[[i,tau]]);
        od;
    od;
    
    # Display(gr_list);
    
    # for t in gr_list do
    #     Display(t[1]);
    # od;
    
    
    indeterminates := DuplicateFreeList(indeterminates);
    gr_list := DuplicateFreeList(gr_list);
    
    # Display(indeterminates);
    
    R := PolynomialRing(Rationals,indeterminates);
    
    list := [];
    
    for i in [1..d] do
        Append(list,Combinations(gr_list,i));
    od;
    
    list := Filtered(list,comb -> Sum(List(comb, x -> x[1])) = d);
    
    expression := 0;
    
    for comb in list do
        expression := expression + 
        Product(List(comb,x -> get_indeterminate(Gr(x[1],x[2][1]*dim_eigenspace(G,chi,x[2])),R)));
    od;
    return [indeterminates,expression,R];
end;;

wait_for_file_to_exist := function(file_name)
    while not IsExistingFile(file_name) do
        Sleep(2);
        Display("waiting...");
    od;
end;;
 
Gr_from_String := function(gr) 
    return SplitString(gr{[4..Length(gr)-1]},",");
end;;

eval_Gr := function(gr,morava_n)
    local file_name, gr_numbers;
    file_name := Concatenation("/Users/Chris/box/research/general_grassmannian/gr_data/",gr,"_",String(morava_n));
    # Display(file_name);
    gr_numbers := Gr_from_String(gr);
    if not IsExistingFile(file_name) then
    Exec(Concatenation(["osascript sage_bridge.scpt ",String(morava_n)," ",gr_numbers[1]," ",gr_numbers[2]]));
    fi;
    wait_for_file_to_exist(file_name);
    return ReadAsFunction(file_name)();
end;;

eval_Gr_min := function(gr,morava_n)
    local file_name, gr_numbers;
    file_name := Concatenation("/Users/Chris/box/research/general_grassmannian/gr_data_min/",gr,"_",String(morava_n));
    # Display(file_name);
    gr_numbers := Gr_from_String(gr);
    if not IsExistingFile(file_name) then
        Exec(Concatenation(["osascript sage_bridge.scpt ",String(morava_n)," ",gr_numbers[1]," ",gr_numbers[2]]));
    fi;
    wait_for_file_to_exist(file_name);
    return ReadAsFunction(file_name)();
end;;

evaluate_G_fixed_points_at_n := function(d,chi,irr_list,G,morava_n)
    local formula_data,evaluation_list;
    formula_data := get_G_fixed_point_formula_for_Gr_d(d,chi,irr_list,G);
    # Display(G);
    Print("Computed at n=",morava_n,"\n");
    Display("Fixed Point Formula:");
    Display(formula_data[2]);
    
    Display(formula_data[1]);
    # Display(List(formula_data[1], x -> get_indeterminate(x,formula_data[3])));
    # Display(List(formula_data[1], gr -> eval_Gr(gr,morava_n)));
    evaluation_list := List(formula_data[1], gr -> eval_Gr(gr,morava_n));

    Display(evaluation_list);
    return Value(formula_data[2],
                 List(formula_data[1], x -> get_indeterminate(x,formula_data[3])),
                 evaluation_list);
end;;

evaluate_min_G_fixed_points_at_n := function(d,chi,irr_list,G,morava_n)
    local formula_data,evaluation_list;
    formula_data := get_G_fixed_point_formula_for_Gr_d(d,chi,irr_list,G);
    # Display(G);
    Print("Min Computed at n=",morava_n,"\n");
    Display(formula_data[2]);
    
    Display(formula_data[1]);
    # Display(List(formula_data[1], x -> get_indeterminate(x,formula_data[3])));
    # Display(List(formula_data[1], gr -> eval_Gr(gr,morava_n)));
    evaluation_list := List(formula_data[1], gr -> eval_Gr_min(gr,morava_n));

    Display(evaluation_list);
    return Value(formula_data[2],
                 List(formula_data[1], x -> get_indeterminate(x,formula_data[3])),
                 evaluation_list);
end;;

test_chi := function(d,G,H,chi,type,drop,coeff,chi_list,chi_H_list)
    local kG,chi_H,kH_n,kH_n_minus_1,condition,folder_name;
    # kG := get_nth_(G,chi,chi_list,type);
    kG := evaluate_G_fixed_points_at_n(d,chi,chi_list,G,type);
    # kG_type_minus_1 := evaluate_G_fixed_points_at_n(d,chi,chi_list,G,type-1);
    
    # Display("=--------=");
    # Display(IsSubgroup(G,H));
    # Display(IsSubgroup(UnderlyingGroup(CharacterTable(G)),H));
    
    chi_H := RestrictedClassFunction(chi,CharacterTable(H));
    
    kH_n := evaluate_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop); 
    kH_n_minus_1 := evaluate_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop-1); 
    
    kH_n_min := evaluate_min_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop); 
    kH_n_minus_1_min := evaluate_min_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop-1); 

    # condition := kH_n_minus_1 < kG and kG <= kH_n;
    condition := kH_n_minus_1_min < kG and kG <= kH_n;
    # condition := condition and (kH_n_minus_1 - kG <= 32);
    # condition := kH_n_minus_1 = kG and kG = 496;
    
    
    # Print("k0: ",kG_type_minus_1,"\n");
    Display(kG);
    Display([kH_n_minus_1,kH_n]);

    Display("Min Ranks:");
    Display([kH_n_minus_1_min,kH_n_min]);
    if condition then
        Display(coeff);
        Display(kG);
        Display([kH_n_minus_1,kH_n]);
        folder_name := String(Concatenation(IdGroup(G),Positions(AllSubgroups(G),H)));
        Exec(Concatenation("mkdir ./output/\"",folder_name,"\""));
        # AppendTo(Concatenation(["./output/",folder_name,"/type_",String(type),"_drop_",String(drop)]),
        #             Concatenation("d=",String(d),String(coeff),"\n",String(kG),"\n",String([kH_n_minus_1,kH_n]),"\n"));
    fi;
    Display(condition);
    return condition;
end;;

record_chi_test_for_graph := function(d,G,H,chi,type,drop,coeff,chi_list,chi_H_list,kG_list,kH_n_list,kH_n_minus_1_list)
    local kG,chi_H,kH_n,kH_n_minus_1,condition,folder_name;
    kG := evaluate_G_fixed_points_at_n(d,chi,chi_list,G,type);
    
    chi_H := RestrictedClassFunction(chi,CharacterTable(H));
    
    kH_n := evaluate_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop); 
    kH_n_minus_1 := evaluate_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop-1); 

    Add(kG_list,kG);
    Add(kH_n_list,kH_n);
    Add(kH_n_minus_1_list,kH_n_minus_1);
end;;

get_real_reps := function(G)
    local list,indicator,real_chi;
    list := Irr(G);
    indicator := Indicator(CharacterTable(G),2);
    real_chi := [];
    for i in [1..Length(list)] do
        if indicator[i] = 1 then
            Append(real_chi,[list[i]]);
        fi;
        if indicator[i]=0 then
            Append(real_chi,[list[i]+ComplexConjugate(list[i])]);
        fi;
    od;
    return DuplicateFreeList(real_chi);
end;;



hunt := function(d,G,H,type,drop,hunt_range) 
    local chi_list, chi_H_list,coeff,chi;
    # chi_list := Irr(G);
    # chi_H_list := Irr(H);
    
    chi_list := get_real_reps(G);
    chi_H_list := get_real_reps(H);
    # chi_H_list := Irr(G);
    
    # # delete this
    # chi_list := Concatenation(Irr(G){[1..8]},Irr(G){[20]});
    # chi_H_list := Irr(H);

    # Display(chi_list);
    # Display(chi_H_list);
    # [1,2,3,16,19]
    for coeff in Tuples(hunt_range,Length(chi_list)) do
        chi := 0;
        for i in [1..(Length(chi_list))] do
            chi := chi + coeff[i]*chi_list[i];
            Display(chi);
            Print(coeff[i],dim_eigenspace(G,chi,chi_list[i]),"\n");
        od;
        chi := Character(CharacterTable(G),chi);
        Display(chi);
        Display("Coefficients");
        Display(coeff);
        if test_chi(d,G,H,chi,type,drop,coeff,chi_list,chi_H_list) then
            break;
        fi;
        # kG := get_nth_k_rank_for_Gr_of_chi(G,chi,chi_list,0);

        # chi_H := restrict_chi_to_H(chi,H,CharacterTable(G));
        # chi_H_list := Irr(H);

        # n := 3;

        # kH_n := get_nth_k_rank_for_Gr_of_chi(H,chi_H,chi_H_list,n);; 
        # kH_n_minus_1 := get_nth_k_rank_for_Gr_of_chi(H,chi_H,chi_H_list,n-1);; 

        # condition := kH_n_minus_1 <= kG and kG < kH_n;

        # if true then
            # Display(coeff);
            # Display(kG);
            # Display([kH_n_minus_1,kH_n]);
            # Display(kH_n_minus_1 <= kG and kG < kH_n);
        # fi;
    od;
Display("done!");
end;;

graph_hunt := function(d,G,H,type,drop)
    local chi_list, chi_H_list,coeff,chi,kG_list,kH_n_minus_1_list,kH_n_list;
    
    chi_list := get_real_reps(G);
    chi_H_list := get_real_reps(H);
    
    kG_list := [];
    kH_n_list := [];
    kH_n_minus_1_list := [];

    for coeff in Tuples([1,2,3],Length(chi_list)) do
        chi := 0;
        for i in [1..(Length(chi_list))] do
            chi := chi + coeff[i]*chi_list[i];
        od;
        chi := Character(CharacterTable(G),chi);
        Display(chi);
        Display("Coefficients");
        Display(coeff);
        record_chi_test_for_graph(d,G,H,chi,type,drop,coeff,chi_list,chi_H_list,kG_list,kH_n_list,kH_n_minus_1_list);
        
        PrintTo("output_graph.sage","load(\"morava_k_theory_of_Gr.sage\")\n","show(my_list_plot(",kG_list,",\"kG\")+my_list_plot(",kH_n_list,",\"kH n\")+my_list_plot(",kH_n_minus_1_list,",\"kH n-1\"))");
    od;
    Display("done!");
end;;

test_coeff := function(d,G,H,type,drop,coeff)
    local chi_list,chi_H_list,chi;
    
    chi_list := get_real_reps(G);
    chi_H_list := get_real_reps(H);
    
    chi:=0;
    for i in [1..(Length(chi_list))] do
        chi := chi + coeff[i]*chi_list[i];
    od;
    chi := Character(CharacterTable(G),chi);
    test_chi(d,G,H,chi,type,drop,coeff,chi_list,chi_H_list);
end;;
# G := DihedralGroup(8);
# chi := 2*Sum(Irr(G))+Irr(G)[5];
# irr_list := Irr(G);

# H := AllSubgroups(G)[3];
# irr_H_list := Irr(H);

# test_chi(1,G,H,chi,0,2,[2,2,2,2,1],irr_list,irr_H_list);


# G := SmallGroup(16,6);
# H := AllSubgroups(G)[3];
# hunt(1,G,H,0,2); # [1..3]


# G = D16 [ 16, 7 ]
# G := SmallGroup(16,7);
# H := AllSubgroups(G)[3];


# chi := 8*get_real_reps(G)[1]+get_real_reps(G)[7];

# chi_list := get_real_reps(G);
# chi_H_list := get_real_reps(H);

# gap> test_chi(3,G,H,chi,0,3,[],chi_list,chi_H_list);

# G := SmallGroup(16,3);
# H := AllSubgroups(G)[1];

# hunt(2,G,H,0,3); #[1..3]

# Coefficients
# [ 8, 8, 8, 8, 8, 8 ]
# Computed at n=0
# 6*Gr(1,8)^2+4*Gr(2,8)+2*Gr(2,16)
# [ "Gr(1,8)", "Gr(2,8)", "Gr(2,16)" ]
# [ 2, 4, 8 ]
# Computed at n=3
# Gr(2,64)
# [ "Gr(1,64)", "Gr(2,64)" ]
# [ 16, 144 ]
# Computed at n=2
# Gr(2,64)
# [ "Gr(1,64)", "Gr(2,64)" ]
# [ 8, 56 ]
# 56
# [ 56, 144 ]
# false
# done!
# gap> 

# This is wrong since empty set
# Coefficients
# [ 1, 2, 2, 2, 1, 1 ]
# Computed at n=0
# 3*Gr(1,1)*Gr(1,2)+3*Gr(1,2)^2+Gr(2,1)+5*Gr(2,2)
# [ "Gr(1,1)", "Gr(2,1)", "Gr(1,2)", "Gr(2,2)" ]
# [ 1, 1, 2, 1 ]
# Computed at n=3
# Gr(2,11)
# [ "Gr(1,11)", "Gr(2,11)" ]
# [ 11, 55 ]
# Computed at n=2
# Gr(2,11)
# [ "Gr(1,11)", "Gr(2,11)" ]
# [ 7, 23 ]
# 24
# [ 23, 55 ]
# [ 1, 2, 2, 2, 1, 1 ]
# 24
# [ 23, 55 ]
# true
# mkdir: ./output/[ 16, 3, 1 ]: File exists
# true
# done!
# gap>

# Coefficients
# [ 7, 7, 7, 7, 7, 7, 7, 7 ]
# Computed at n=0
# 6*Gr(1,7)^2+4*Gr(2,7)+2*Gr(2,28)+2*Gr(2,14)
# [ "Gr(1,7)", "Gr(2,7)", "Gr(2,28)", "Gr(2,14)" ]
# [ 1, 3, 14, 7 ]
# Computed at n=2
# Gr(2,84)
# [ "Gr(1,84)", "Gr(2,84)" ]
# [ 8, 66 ]
# Computed at n=1
# Gr(2,84)
# [ "Gr(1,84)", "Gr(2,84)" ]
# [ 4, 46 ]
# 60
# [ 46, 66 ]
# [ 7, 7, 7, 7, 7, 7, 7, 7 ]
# 60
# [ 46, 66 ]
# true
# mkdir: ./output/[ 16, 3, 1 ]: File exists
# true
# done!
# gap> 

# Coefficients Range [2,3]
# [ 2, 2, 2, 2, 2, 3, 2, 2 ]
# Computed at n=0
# 6*Gr(1,2)^2+2*Gr(2,4)+4*Gr(2,2)+Gr(2,12)+Gr(2,8)
# [ "Gr(1,2)", "Gr(2,2)", "Gr(2,8)", "Gr(2,12)", "Gr(2,4)" ]
# [ 2, 0, 4, 6, 2 ]
# Computed at n=3
# Gr(2,26)
# [ "Gr(1,26)", "Gr(2,26)" ]
# [ 16, 125 ]
# Computed at n=2
# Gr(2,26)
# [ "Gr(1,26)", "Gr(2,26)" ]
# [ 8, 37 ]
# 38
# [ 37, 125 ]
# [ 2, 2, 2, 2, 2, 3, 2, 2 ]
# 38
# [ 37, 125 ]
# true
# mkdir: ./output/[ 16, 3, 1 ]: File exists
# true
# done!
# gap> hunt(2,G,H,0,3);

#
# G := DihedralGroup(8);
# H := AllSubgroups(G)[3];

hunt_for_Gr_with_trival_higher_K2 := function(d,hunt_range)
    local G,H, chi_list, chi_H_list,coeff,chi,such_Gr_list;
    
    G := DihedralGroup(8);
    H := AllSubgroups(G)[1];

    chi_list := get_real_reps(G);
    chi_H_list := get_real_reps(H);
    
    such_Gr_list := [];
    
    
    for coeff in Tuples(hunt_range,Length(chi_list)) do
        chi := 0;
        for i in [1..(Length(chi_list))] do
            chi := chi + coeff[i]*chi_list[i];
            Print(coeff[i],dim_eigenspace(G,chi,chi_list[i]),"\n");
        od;
        chi := Character(CharacterTable(G),chi);
        Add(such_Gr_list,test_chi_for_Gr_higher_diffs(d,G,H,chi,0,3,coeff,chi_list,chi_H_list));
    od;
    
    such_Gr_list := Filtered(such_Gr_list, x -> (x[1] = true));
    such_Gr_list := DuplicateFreeList(such_Gr_list);
    Display(List(such_Gr_list,x->x[2]));
    
    
end;;

test_chi_for_Gr_higher_diffs := function(d,G,H,chi,type,drop,coeff,chi_list,chi_H_list)
    local kG,chi_H,kH_n,kH_n_minus_1,condition;
    # kG := get_nth_(G,chi,chi_list,type);
    kG := evaluate_G_fixed_points_at_n(d,chi,chi_list,G,type);
    
    # Display("=--------=");
    # Display(IsSubgroup(G,H));
    # Display(IsSubgroup(UnderlyingGroup(CharacterTable(G)),H));
    
    chi_H := RestrictedClassFunction(chi,CharacterTable(H));
    
    kH_n := evaluate_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop); 
    kH_n_minus_1 := evaluate_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop-1); 
    
    kH_n_min := evaluate_min_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop); 
    kH_n_minus_1_min := evaluate_min_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop-1); 

    # condition := kH_n_minus_1 <= kG and kG < kH_n;
    condition := (kH_n_minus_1 = kG and kG <= kH_n) and (kG <= kH_n_min) ;
    
    
    if condition then
        Display(coeff);
        Display(kG);
        Display([kH_n_minus_1,kH_n]);

        Display("Min Ranks:");
        Display([kH_n_minus_1_min,kH_n_min]);
        Display(get_G_fixed_point_formula_for_Gr_d(d,chi,chi_H_list,H)[2]);
    fi;
    return [condition,get_G_fixed_point_formula_for_Gr_d(d,chi,chi_H_list,H)[2]];
end;;


# gap> hunt_for_Gr_with_trival_higher_K2(2,[1,2,3,4,5]);
# [ Gr(2,10), Gr(2,12), Gr(2,14), Gr(2,16), Gr(2,18), Gr(2,11), Gr(2,13), 
#   Gr(2,15), Gr(2,17), Gr(2,19), Gr(2,20), Gr(2,21), Gr(2,22), Gr(2,23), 
#   Gr(2,24), Gr(2,25), Gr(2,26), Gr(2,27) ]

# SmallGroup(16,3);
# Coefficients
# [ 2, 2, 2, 2, 2, 3, 3, 3 ]
# Computed at n=1
# Fixed Point Formula:
# 6*Gr(1,2)^2+Gr(2,4)+4*Gr(2,2)+3*Gr(2,6)
# [ "Gr(1,2)", "Gr(2,2)", "Gr(2,4)", "Gr(2,6)" ]
# [ 2, 1, 6, 7 ]
# Computed at n=0
# Fixed Point Formula:
# 6*Gr(1,2)^2+Gr(2,4)+4*Gr(2,2)+3*Gr(2,6)
# [ "Gr(1,2)", "Gr(2,2)", "Gr(2,4)", "Gr(2,6)" ]
# [ 2, 1, 2, 3 ]
# Computed at n=4
# Fixed Point Formula:
# Gr(2,30)
# [ "Gr(1,30)", "Gr(2,30)" ]
# [ 30, 435 ]
# Computed at n=3
# Fixed Point Formula:
# Gr(2,30)
# [ "Gr(1,30)", "Gr(2,30)" ]
# [ 16, 127 ]
# Min Computed at n=4
# Gr(2,30)
# [ "Gr(1,30)", "Gr(2,30)" ]
# [ 30, 435 ]
# Min Computed at n=3
# Gr(2,30)
# [ "Gr(1,30)", "Gr(2,30)" ]
# [ 16, 41 ]
# k0: 39
# 55
# [ 127, 435 ]
# Min Ranks:
# [ 41, 435 ]
# [ 2, 2, 2, 2, 2, 3, 3, 3 ]
# 55
# [ 127, 435 ]
# mkdir: ./output/[ 16, 3, 1 ]: File exists
# true
# done!
# gap> 
