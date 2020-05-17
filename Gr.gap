LoadPackage("ctbllib");

# The type parameter is the Schur indicator of the real rep
Gr := function(d,m,type) 
    return Concatenation(["Gr","(",String(d),",",String(m),",",String(type),")"]);
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

get_real_rep_type := function(G,chi)
    local indicator;
    indicator := Indicator(CharacterTable(G),[chi],2)[1];
    if indicator = 0 then
        return "C";
    elif indicator = 1 then
        return "R";
    fi;
    return "H";
end;;

get_indeterminate := function(name,R)
    local x;
    for x in IndeterminatesOfPolynomialRing(R) do
        if String(x) = name then
            return x;
        fi;
    od;
end;;
distinct_comb := function(comb)
    local chis;
    chis := List(comb, x -> x[2]);
    return IsDuplicateFree(chis);
end;;

get_G_fixed_point_formula_for_Gr_d := function(d,chi,irr_list,G)
    local R, indeterminates, i, tau, dim, list, comb, expression, gr_list;
    indeterminates := [];
    
    gr_list := [];
    
    for tau in irr_list do
        dim := tau[1];
        # Display(tau[1]);
        # Display(d);
        for i in [dim..d] do
            # Display(dim_eigenspace(G,chi,tau));
            if i mod dim = 0 and dim_eigenspace(G,chi,tau) <> 0 then
                # Display(dim_eigenspace(G,chi,tau));
                Append(indeterminates,[Gr(i / dim,dim_eigenspace(G,chi,tau),get_real_rep_type(G,tau))]);
                Append(gr_list,[[i,tau]]);
                # Display(indeterminates);
            fi;
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

    # ensures each rep is contribing one factor in the product
    list := Filtered(list,comb -> distinct_comb(comb));
    
    expression := 0;
    
    for comb in list do
        expression := expression + 
        Product(List(comb,x -> get_indeterminate(Gr(x[1] / x[2][1],dim_eigenspace(G,chi,x[2]),get_real_rep_type(G,x[2])),R)));
    od;
    return [indeterminates,expression,R];
end;;

pad_list_with_zeros := function(list,length)
   while Length(list) < length do
       Add(list,0);
   od;
   return list;
end;;

composistions := function(d,length)
    local comp,final,x;
    comp := Partitions(d);
    # comp := List(comp, x -> pad_list_with_zeros(x,length));
    comp := Filtered(comp, x -> Length(x)=length);
    comp := List(comp,x -> PermutationsList(x));
    final := [];
    for x in comp do
        Append(final,x);
    od;
    return final;
end;; 

wait_for_file_to_exist := function(file_name)
    while not IsExistingFile(file_name) do
        Sleep(2);
        Display("waiting...");
        Display(file_name);
    od;
end;;

gr_formula := function(q,m,n)
    # local k,i;
    # if n = 0 then
    #     if q mod 2 = 0 then
    #         k = q/2;
    #         if m mod 2 = 0 then
    #             i := m/2;
    #         else
    #             i := (m-1)/2;
    #         fi;
    #         return binomial(i,k)+1;
    #     fi;
    # fi;
    if m mod 2 = 0 then
        return gr_even_formula(q,m,n);
    fi;
    return gr_odd_formula(q,m,n);
end;;

gr_even_formula := function(q,m,n)
    local l,rank,part,i,j,k,perm;
    if m < q then
        return 0;
    fi;
    if q <= m and m <= 2^(n+1) then
        return Binomial(m,q);
    fi;

    l := (m-2^(n+1))/2;
    rank := 0;
    for part in Filtered(composistions(q,2), x-> x[1] mod 2=0) do
        i := part[1]/2;
        for perm in composistions(part[2],2) do
            j := perm[1];
            k := perm[2];
            rank := rank + Binomial(l,i)*Binomial(2^n,j)*Binomial(2^n,k);
        od;
    od;
    for part in composistions(q,2) do
        j := part[1];
        k := part[2];
        rank := rank + Binomial(2^n,j)*Binomial(2^n,k);
    od;
    for part in Filtered(composistions(q,2), x-> x[1] mod 2=0) do
        i := part[1]/2;
        j := part[2];
        rank := rank + 2*Binomial(l,i)*Binomial(2^n,j);
    od;
    rank := rank + 2*Binomial(2^n,q);
    if q mod 2 = 0 then
        rank := rank + Binomial(l,q/2);
    fi;
    return rank;
end;;

gr_odd_formula := function(q,m,n)
    local l,rank,part,i,j,k,perm;
    if m < q then
        return 0;
    fi;
    if q <= m and m <= 2^(n+1) then
        return Binomial(m,q);
    fi;

    l := (m+1-2^(n+1))/2;
    rank := 0;
    for part in Filtered(composistions(q,2), x-> x[1] mod 2=0) do
        i := part[1]/2;
        for perm in composistions(part[2],2) do
            j := perm[1];
            k := perm[2];
            rank := rank + Binomial(l,i)*Binomial(2^n-1,j)*Binomial(2^n,k);
        od;
    od;
    for part in composistions(q,2) do
        j := part[1];
        k := part[2];
        rank := rank + Binomial(2^n-1,j)*Binomial(2^n,k);
    od;
    for part in Filtered(composistions(q,2), x-> x[1] mod 2=0) do
        i := part[1]/2;
        j := part[2];
        rank := rank + Binomial(l,i)*Binomial(2^n-1,j);
        rank := rank + Binomial(l,i)*Binomial(2^n,j);
    od;
    rank := rank + Binomial(2^n-1,q);
    rank := rank + Binomial(2^n,q);
    if q mod 2 = 0 then
        rank := rank + Binomial(l,q/2);
    fi;
    return rank;
end;;
 
Gr_from_String := function(gr) 
    return SplitString(gr{[4..Length(gr)-1]},",");
end;;

gr1_formula := function(m,n)
    local l; # Gr(1,m)=RP^(m-1)=RP^(l-1)
    l:= m;
    if m = 0 then
        return 0;
    fi;
    if m = 1 then
        return 1;
    fi;
    if l < 2^(n+1) then
        return l;
    elif (l mod 2) = 0 then
        return 2^(n+1);
    elif (l mod 2) = 1 then
        return 2^(n+1)-1;
    fi;
end;;

gr2_formula := function(m,n)
    local l;
    if n = 0 then
        if (m mod 2) = 0 then
            return m /2;
        else 
            return (m-1)/2;
        fi;
    fi;
    if m = 1 then
        return 1;
    fi;
    if 2 <= m and m <= 2^(n+1) then
        return Binomial(m,2);
    elif (m mod 2) = 0 then
        l := (m - 2^(n+1))/2;
        return 2^(2*n+1)-2^n+l;
    elif (m mod 2) = 1 then
        l := (m - 2^(n+1)+1)/2;
        return 2^(2*n+1)-2^(n+1)-2^n+1+l;
    fi;
end;;

eval_Gr := function(gr,morava_n)
    local file_name, gr_numbers;
    file_name := Concatenation("/Users/Chris/box/research/general_grassmannian/gr_data/",gr,"_",String(morava_n));
    # Display(file_name);
    gr_numbers := Gr_from_String(gr);
    
    if gr_numbers[3] = "C" then
        return Binomial(Int(gr_numbers[2]),Int(gr_numbers[1]));
    fi;
    
    if gr_numbers[3] = "R" then
        # if gr_numbers[1] = "1" then
        #     return gr1_formula(Int(gr_numbers[2]),morava_n);
        # elif gr_numbers[1] = "2" then
        #     return gr2_formula(Int(gr_numbers[2]),morava_n);
        # fi;
        return gr_formula(Int(gr_numbers[1]),Int(gr_numbers[2]),morava_n);
    fi;
    
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
    
    # kH_n_min := evaluate_min_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop); 
    # kH_n_minus_1_min := evaluate_min_G_fixed_points_at_n(d,chi_H,chi_H_list,H,type+drop-1); 

    condition := kH_n_minus_1 < kG and kG <= kH_n;
    # condition := kH_n_minus_1_min < kG and kG <= kH_n;
    # condition := condition and (kH_n_minus_1 - kG <= 32);
    # condition := kH_n_minus_1 = kG and kG = 496;
    
    
    # Print("k0: ",kG_type_minus_1,"\n");
    Display(kG);
    Display([kH_n_minus_1,kH_n]);

    # Display("Min Ranks:");
    # Display([kH_n_minus_1_min,kH_n_min]);
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
