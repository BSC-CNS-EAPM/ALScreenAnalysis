## ALScreen pseudocode

ScoreFormer learns to predict docking scores from molecular structures. 
In each active-learning cycle, ALScreen docks the top-ranked compounds, 
adds examples with both low and high prediction errors to the training pool, and retrains ScoreFormer.

The pseudocode is available in [Markdown](alscreen_pseudocode.md) and [PDF](ALScreen_pseudocode.pdf).

[![ScoreFormer training](ALScreen_pseudocode_page_1.png)](ALScreen_pseudocode.pdf)

[![ScoreFormer inference](ALScreen_pseudocode_page_2.png)](ALScreen_pseudocode.pdf)

[![ALScreen active-learning cycle](ALScreen_pseudocode_page_3.png)](ALScreen_pseudocode.pdf)
