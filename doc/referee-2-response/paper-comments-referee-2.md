# Referee 2 comments
## Overview

The authors present a study on gradient based wind-farm layout optimization and comparison of the optimized results with LES. An important novelty claimed by the authors is the fact that optimized layouts are compared to 'high-fidelity' LES over a full wind rose. Unfortunately, there is potentially a major flaw in the presented LES results (see point 1 below). Therefore I do not believe this work can be accepted for publication in its current form

## Detailed comments

1. LES are performed on a domain of 5x5x1 km. This domain is much too small, and given the forced inflow conditions, will lead to domain blockage and an artificial favorable pressure gradient. This is essentially a direct result from Newtons second law when domain boundaries do not allow the momentum of the flow to freely change. When optimizing the layout, and thus changing to total wind-farm thrust, this will then also lead to a larger favorable pressure gradient, hence artificially enhancing the benefits of the layout in the LES. This seems to be exactly what the authors observe when comparing their LES with the optimized wake models. Also the fact that gains that are calculated using the front row turbine's velocities as reference (as in Fig 13) are closer to the wake models than direct power (which is based on inflow reference) points to significant blockage effects. The authors should show in a revised manuscript that their selected domain size is sufficiently big to avoid blockage effects that are of the same order of magnitude as optimization gains. To this end, they should for a selected case show results on different domain sizes, showing that effects become negligible for the final selected size. It is my expectation that that size is considerably bigger than 5x5x1 km. Subsequently, all simulations should be performed on that domain size.
 
2. Overall, the discussion on the LES set-up in 2.4 and 2.5 is too brief. Based on this, the set-up is simply not reproducible. The authors mention buoyancy and Coriolis effects, but it is my impression that the simulations may simply consist of a pressure driven boundary layer. If not, what is the geostrophic wind, what is the stratification profile. In addition, what is the surface roughness, friction velocity, etc. What is the precursor set-up (domain size, grid size, initialization, spin-up time, etc). What is the simulation cost, …

## Smaller comments

1. [x] Line 84 and also later section 2.1: better justify why the near-wake region can not be simply avoided by using a minimum distance constraint in the optimization, e.g. using a constraint that is larger than x_d obtained Eq.8. In case of interest in set-ups where turbines are placed more closely together, the near wake model may need improvement anyway, and the heuristic adaptation proposed may not suffice.
    - A constraint on the minimum turbine distance is included, but the size of the near wake region is variable and the turbine separation constraints are non-linear. SNOPT minimizes infeasibility for non-linear constraints, so infeasible solutions may be attempted during the optimization and the near wake model may be used at some point during the optimizaiton.
    - We have added clarification at 
        - line 85-90
        - line 141
        - line 153-155
2. [x] Throughout the paper: equations are part of the text, and phrases and punctuations should be used accordingly. Please check with other papers and publication standards to see how it is done
- We have made this change. Thank you for the suggestion.
3. [x] Line 143: “To remove the discontinuity” --> please provide a mathematical expression
- Good suggestion. We have added clarification
    - lines 150-157
4. [x] Line 152: If greater accuracy is desired --> speculative. Either provide data that prove this statement or remove
- We have adjusted this statement slightly to clarify that we are not making a claim, but rather providing a suggestion and reference for readers interested in a more accurate near-wake model.
    - see lines 165-167
5. [x] Line 212: sunflower pattern: please provide reference or formula
- Thank you for the suggestion. The formula, along with a citation and explanation, have been added. 
6. [x] Eq 19: provide units. Result is in kWh and not in J.
- Thank you for the suggestion, we have clarified the units.
7. Eq 21: this is not a correct definition of TI. TI is based on magnitude of fluctuating velocity that includes components in all directions
8. [x] Eq 22: to be technically precise, the “i=1…38” should be added in the subscript below “maximize”
- Thank you for pointing this out. To make things cleaner we have just removed the underset (x_i,y_i) because they are redundant with the shown objective inputs.
9. [x] Line 300: forward differentiation is used. What is the advantage of this over using SNOPT without providing the gradients explicitly? I do not believe that this will be significant, since in that case SNOPT constructs gradients based on FD? This can be even less expensive than forward differentiation, and accuracy loss is often not significant (depending on implementation choices). Please discuss in more detail the gains etc. Usually, significant speed-up would only follow from backward differentiation. This should be better substantiated in the manuscript, in particular since, in the abstract, you seem to claim this as an important innovation.
- Thank you for bringing this up. Our explanation was not clear enough. ForwardDiff.jl does not use the forward finite difference method, but rather forward mode algorithmic differentiation (AD). We have added more detail in section 2.7 to clarify this point. AD has been shown to much less computationally costly and more accurate than finite difference methods, such as the one provided in SNOPT.
10. 
    1. [x] Line 305: please discuss the WEC method in more detail. 
        - We have added more detail about the WEC method in section 2.7, as requested.
    2. [x] Also better explain why standard multistart methods do not work? If they are as good, why not use a standard from an optimization library
        - WEC is a continuation method, not a multistart method. However, using it with a multistart method is recommended. We have provided an explanation of our multi-start approach in section 2.7
11. [x] Line 314: what do you mean with 400 optimizations? This is confusing. If I’m not mistaken, you solve Eq 22 only twice. Better clarify/make distinction
- This should be more clear with our additions made to section 2.7
12. [x] Figures in general: please make as much as possible black&white friendly (some plots are not readable when printed in gray scale). For many figures this should be possible without losing attractiveness of the figure in color.
- While we recognize that the color scheme may not be ideal for printing, we would like to leave it as is. The colors were selected for colorblind readers. Since most readers will be seeing the article online, we believe that having a colorblind-friendly palet is more important than having a printer friendly one.
13. Wind directions throughout paper: add degree symbol