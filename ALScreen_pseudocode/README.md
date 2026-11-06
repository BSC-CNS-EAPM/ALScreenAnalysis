## ALScreen pseudocode

ScoreFormer learns to predict docking scores from molecular structures. 
In each active-learning cycle, ALScreen docks the top-ranked compounds, 
adds examples with both low and high prediction errors to the training pool, and retrains ScoreFormer.

[Open the full pseudocode (PDF)](pseudocode/alscreen_pseudocode.pdf)

[![ScoreFormer training](pseudocode/alscreen_pseudocode-1.png)](pseudocode/alscreen_pseudocode.pdf)

[![ScoreFormer inference](pseudocode/alscreen_pseudocode-2.png)](pseudocode/alscreen_pseudocode.pdf)

[![ALScreen active-learning cycle](pseudocode/alscreen_pseudocode-3.png)](pseudocode/alscreen_pseudocode.pdf)
