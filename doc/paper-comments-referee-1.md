# Response to Referee 1

- "The manuscript presents a rigorous methodology for wind farm layout optimization using simplified engineering models. The main strength of the study is in comparing their model’s output against large eddy simulation results. The text is very well written, with the methodology and results presented in a clear and concise manner to aid in reproducibility. The following points could improve the readability of the text."
    - Thank you for your positive feedback and suggestions.

- "In section 2.2.3 the author detail that they used only 1 sampling point for velocity calculation during the optimization procedure, and highlight the errors in the 1 point case in Figure 3. Readers would benefit from seeing the impact of this simplification on the optimal results through some test cases where higher sampling points were taken during optimization as well."
    - Excellent suggestion. We have added this study and it looks like one sample is sufficient for optimization, while 20 samples are enough for AEP estimation.
        - See page 17, lines 362-375
        - See page 19, Fig. 12

- "In section 2, it may be beneficial for completeness and readability to expand a bit more upon the wake expansion continuation method and the meanings of the associated relaxation factors rather than just mentioning it as a reference."
    - Thank you for the suggestion. We have added a brief description of WEC. Please see page 14, lines 306-314

- "Figure 12 show significant differences between the BP model and SOWFA, for both the base and optimized results. Earlier in the paper in section 3, the authors showed that their model is able to match reference LES results for the horns rev wind farm with a high degree of accuracy. Could the authors comment on why there are now larger errors when compared to LES for the developed wind farm layout?"
    - The results from the literature were normalized. We never normalized our results.
    - When we used a non-dimensional form, the agreement was better (see fig 14).
    - Clarification has been added in the following places:
        - page 18, line 379
        - page 19, lines 393-394