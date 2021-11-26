Nice paper!  Good analysis and well written.  Minor comments below:

X   not sure if validated is the right word for the title.  verification isn’t quite right either.  Maybe just Compared?

X   I assume you need to rewrite the abstract?  The 9.9 and 15.6 look like old numbers, plus there are two sets now.

X   “The model and LES results were then compared” - probably can combine with the previous sentence.  felt a bit choppy.

X   lit review and intro was nice.

X   Fig 1, capitalize TI

X   top of pg 7 typo in Hernandez

X   The algorithm for calculating  - 1 sentence paragraph

X   hard maximum function is given in Alg. 1 - you never say what this is. or later in that paragraph when you mention smooth max therer isn’t any motivation for why you need that.

X   same paragraph you say two different approaches, but in next sentence you say three methods.

X   Is Alg 1 necessary to write out? Sounds like it is pretty much the same as the existing approach except for a smooth max?  or is there more changes?

X   2.2.3 this wasn’t clear until I read it again, but when it ended with 1 pt I didn’t get the impression that you went back to 100 pts for comparisons.  Maybe it’s fine as is though.  I don’t have an obvious suggestion.

X   Fig 3 - I’m not clear about what’s going on here.  What is the case w/o local TI?  Do you ever use that?  I don’t remember you talking about it.  If not, why is it in the plot?  If so, how do you justify it as it looks pretty inaccurate?

X   Same figure - probably need a box around the legend.  Often that’s not desirabl, but the labels are large and the plot a bit cramped that the labels are somewhat confusing as data points.

X   You neve say what Fig 3 and Fig 4 aree.  You say 3a and 4a are 100 samples, and 3b, 4b, 1 sample.  But you never say 3 is power along columns and fig 4 is power by direction.

X   The logical flow was a bit confusing.  You said 1 sample resulted in large errors then in next sentence we determined our implementation was sufficient.  those seem contradictory - didn’t seem to follow from above sentence. 

X   m s^{-1} seems pretty awkward as compared to m/s  (lots of instances of this.

X   Fig 5 and 6 could be combined side by side.  Lots of space for relatively simple info.

X   Fig 7 - it’s not clear why these are numbered.  I think I realized later that you use the numbers in your big table plot.  But maybe you should  mention that upfront that you’ll refer to the numberrs later

X   julia should be capitalized

X   For the julia packages without citatinos you should put a footnote with the url (FLOWFarm, SNOW, etc.).

X   However, the final optimized…are within 1 GW h of each other.  I didn’t understand this sentence and how it contrasted from previous.

X   I don’t see the hexagons - some of this discussion wasn’t clear to me.

X   “aligned with the 12 wind directions” - how can they be aligned with 12 directions?

X   “along rows perpendicular to that side’s direct wind” - what does this mean?

X   “turbines with the wind” - turbines when the wind

X   the optimization favored straight lines aligned with the wind - why?  that doesn’t make sense to me.

X   I think there are some good insights in this paragraph describe the optimized results, but I’m not quite following it all. Maybe needs a little more explanation.

X   Fig 10 probably reduce to just the optimized.  thtat’s the only thing you talk about anyway.  Base is already shown and start doesn’t really tell me much and it isn’t discussed.

X  Fig 10, you don’t refer to the subfigs (a, b, c, )

X   “all the turbines in the farm that were not in any other turbines’ wakes” - what does this mean? 

X   Fig 11 - caption doesn’t refer to a, b, c,

X   we can now see that 220 - mention that this is the most probably direction.  You say that later, but that isn’t clear at this point in the paper

X   actual power predictions we significantly - we -> were

X   under predicted - I think this can be one word

X   BP predicted 7.7 for high and 10.2 for low.  SOWA 9.3 and 10.6.  This would be clearer if you listed this by TI case instead of by method.  When I read this the numbers I wanted to compare was BP to SOWFA. In other words 7.7 to 9.3 and 10.2 to 10.6. Comparing 7.7 and 10.2 is irrelevant, but the way it’s written makes those easily to compare. 

X   same comment in conclusioins

X   I’d remove Table 1 and 2.  I don’t see what they add that isn’t in plots.

X   the [individual] wind turbine power errors are - add a word for clarity

X   “secured funding and provided support” - probably not the most relevant wording.  Kind of weird to mention the funding part especially when those who funded it are on the paper.  and support is too vague.  I’d say something like: provided direction, ideas, and feedback or something like that. 