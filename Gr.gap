LoadPackage("ctbllib");

Gr := function(d,m) 
    return Concatenation(["Gr","(",String(d),",",String(m),")"]);
end;;

dim_eigenspace := function(G,chi,lambda)
    local t;
    t := CharacterTable(G);
    return Sum( [ 1 .. Length(ConjugacyClasses(t)) ], i -> SizesConjugacyClasses( t )[i] * chi[i]
                            * ComplexConjugate( lambda[i] ) ) / Size( t );
end;;

get_indeterminate := function(name,R)
    local x;
    for x in IndeterminatesOfPolynomialRing(R) do
        if String(x) = name then
            return x;
        fi;
    od;
end;;

get_G_fixed_point_formula_for_Gr_d := function(d,chi,irr_list,G)
    local R, indeterminates, i, tau, dim, list, comb, expression;
    indeterminates := [];
    
    for tau in irr_list do
        dim := tau[1];
        for i in [dim..d] do
            Append(indeterminates,[Gr(i,dim*dim_eigenspace(G,chi,tau))]);
        od;
    od;
    
    indeterminates := DuplicateFreeList(indeterminates);
    
    R := PolynomialRing(Rationals,indeterminates);
    
    list := [];
    
    for i in [1..d] do
        Append(list,Combinations(irr_list,i));
    od;
    
    Display(list);
    
    list := Filtered(list,comb -> Sum(List(comb, x -> x[1])) = d);
    
    expression := 0;
    
    for comb in list do
        expression := expression + 
        Product(List(comb,x -> get_indeterminate(Gr(d+1-Length(comb),x[1]*dim_eigenspace(G,chi,x)),R)));
    od;
    return [indeterminates,expression];
end;;


G := DihedralGroup(8);
chi := Sum(Irr(G));
irr_list := Irr(G);

H := AllSubgroups(G)[2];
irr_H_list := Irr(H);

chi_H := RestrictedClassFunction(chi,CharacterTable(H));
