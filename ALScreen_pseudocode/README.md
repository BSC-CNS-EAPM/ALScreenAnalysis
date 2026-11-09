## ALScreen pseudocode

ScoreFormer learns to predict docking scores from molecular structures. 
In each active-learning cycle, ALScreen docks the top-ranked compounds, 
adds examples with both low and high prediction errors to the training pool, and retrains ScoreFormer.

The pseudocode is also available in Markdown ([alscreen_pseudocode.md](alscreen_pseudocode.md)) and PDF ([Alscreen_pseudocode.pdf](Alscreen_pseudocode.pdf)) formats.

[![ScoreFormer training](Alscreen_pseudocode_page_1.png)](Alscreen_pseudocode.pdf)

[![ScoreFormer inference](Alscreen_pseudocode_page_2.png)](Alscreen_pseudocode.pdf)

[![ALScreen active-learning cycle](Alscreen_pseudocode_page_3.png)](Alscreen_pseudocode.pdf)
