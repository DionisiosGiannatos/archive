function output=dissociation_finder(K, init_guess,Is_Acid=true ,c1 =1,c2 =0,c3=10^(-7),Kw=10^(-14))
  #A fun little script I made, to calculate the dissociation of an acid or base by
  #solving the equation through numerical methods.
  
  #calculator of dissasociation and degree of dissassociation through inputs
  #of the reaction constant K and the concentrations of components 1,2 and 3. Inputs considered in M
  # The reaction is written as following: HA (1) + H2O <-> A- (2)+ H3O+ (3)
  #or, B (1) + H2O <-> HB+ (2) + OH- (3)
  #Kw is taken by default to be equal to 10^(-14), can be changed through args.
  
  #Self ionization of water included in c3 for extremely dilute solutions,
  #now it should give the correct result, regardless of concentration
  
  #The calculator assumes, by default, that a base or acid is dissociated in a
  #pure mixture. If not, change c2 and c3 accordingly
  #C2 expresses concentration of Ions of conjugate base, c3 of preexisting pH or pOH.
  
  #parser of Is_Acid. Any value other than 0 is taken as true.
  Is_Acid=logical(Is_Acid);
  
  #checks whether the mixture is a base or not.
  Is_Base=not(Is_Acid)
  
  #some initial vars for the NR method.
  format long e
  a=init_guess;
  i=0;
  a_new=2;
  a_old=1;
  
  #neuton raphson method for finding roots. It is a quick method,
  #therefore if the number of iterations surpasses 200, then it is assumed that
  #a root doesn't exist.
  while (abs(a_new-a_old)/(a_old)>10^-6) && (i<200)
    i=i+1;
    a_old=a;
    f=a^(2)+a*(c2+c3+K)+(c2*c3-K*c1);
    fder=2*a+(c2+c3+K);
    a=a-(f/fder);
    a_new=a;
  endwhile
  
  if i>200 
    #a check to see if NR could actually find a root. If not, it displays an error
    error("A Root couldn't be found through the Newton-Raphson method.")
  endif
  
  #pH Calculation depends on whether the substance was an acid or not.
  #if a base is dissociated, pH needs to be calculated from pKw
  if (Is_Acid == true)
    #c3+a is taken in case mixture exists in a pre-existing pH value.
    #by default, c3 expresses the self-ionization of water(c=10^-7)
    #making this calculator accurate at extremely dilute concentrations.
    pH=-log10(c3+a);
  elseif (Is_Base == true)
    pOH=-log10(c3+a);
    pH=-log10(Kw)-pOH;
  endif
  #degree of dissociation calculator as a percentage. Expresses amount of
  #original substance dissociated.
  deg_dis=(a/c1)*100;
  #output field
  disp("output: pH, Degree of Dissociation, Amount dissociated, Iterations")
  output=[pH,deg_dis,a,i];
endfunction
