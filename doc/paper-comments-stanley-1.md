Alright! This paper is so cool, again good work. I’m glad I was able to participate. Here are a couple of specific comments, of course I leave it up to you if you want to incorporate them or not.
 
 
# Abstract:
X   I added/changed a couple of things. Adopt/modify as you see fit
    
    The physics models used during wind farm layout optimization use simplifying assumptions that can alter the design space **compared to reality and higher-fidelity simulations**. Some characteristics of these simple models may negatively influence the resulting layouts. In this paper, we perform wind farm layout optimization **using a simple engineering wake model** and then simulate the base and optimized layouts using large-eddy simulation (LES) to confirm that the layout was actually improved, **and not just an artifact of the simplifying assumptions in the low-fidelity wind farm simulation**. We begin by describing the physics models used, including changes specific for use with gradient-based optimization. We then compare the simple model’s output to previously published model and LES results. Using the simple models described, we performed gradient-based wind farm layout optimization using exact gradients. We optimized the wind farm twice, with high and low turbulence intensity (TI) respectively. We then recalculated annual energy production (AEP) using LES for the original and optimized layouts in each TI scenario and compared the results. For the high-TI case, the simple model predicted an improvement in AEP of 7.7%, while the LES reported an AEP improvement of 9.3%. For the low-TI case, the simple model predicted an improvement in AEP of 10.0%, while the LES reported an AEP improvement of 10.7%. We concluded that the improvements found by optimizing with the simple model are not just an artifact of the model, but are real improvements. 
 
# Intro:
X    “that could potentially reduce the optimality of the final layout.” - I get what you’re going for, I think it can be more explicit though. The optimizer could exploit weaknesses in the simple wake models to find a layout that is good for the simple wake model but not ACTUALLY good in real life. I see that you say this more clearly at the end of the paragraph, but still think it should be mentioned early on as well.
 
    Great work, I love the intro
 
 
# Section 2
    X Honest question, it may not matter, will people want an explanation of alpha and beta star (line120)? Maybe not, whatever you think
    
    X Another suggestion, depending on what you think: Whenever I’m reading NP my brain goes to NP-hard. Maybe another abbreviation? NPA? 
    - kept NP because it is parallel with the already prevalent BP used for the gaussian wake model

    X Also, in a few places you say NP and some Niayifar and Porté-Agel, I’d make sure you’re consistent.
    - adjusted the wording in a few places to make clear that NP refers to the model, not the poeple
        or papers. The full names are used for reference to the people and papers.
 
 
# Section 3:
 
X “each optimized layout each” sounds weird, should it be IN each optimized layout?
 
 
X “…directional annual energy improvement. How much the…” This sounds weird, maybe make it one sentence? …improvement, which we defined how much the…
 
 
X Figures 13 and 14 need some cleanup. (I can do it if that would be helpful, just send me the code/data. Although the changes I think are relatively small, it may be easier for you to just do) A couple of things, 
    X black borders don’t go all the way around
        - I'll let PJ take a look, but I don't think it is worth fixing (I tried a bit)
    X turbine index numbers on the top run together, negative double digit numbers in the cells run beyond the borders, and negative numbers run together on the color bar. 
    X the titles look weird being so big, but that is probably just a preference thing. 
 
# Then lastly,
 
.   “It is likely that a more refined wind rose would reduce the AEP gains of these optimized layouts.” Maybe opening a can of worms…but is it worth running this? With a more refined wind rose and see what they look like? I would argue that a full analysis of that would be beyond scope, but just running it with the full wind rose and reporting how the AEP changes could be interesting.
 
From one point of view, it could be an interesting sneak peak at work to come (by you or whoever), about the wind rose fidelity. Like, you could say “For example, with the BP wake model, there is an AEP reduction of X% when evaluating with 1/2/5 degree bins compared to the 12 bins used in this study. This exploration of wind rose fidelity is beyond the scope of this paper, but we recognize it as an important area of future work.”
 
That was really spur of the moment, obviously it can probably be worded better. But something to think about.
 
 
Great work, great paper,
 
PJ