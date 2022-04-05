# Vienna-RNAfold
Bachelor Thesis 2022
Using the VienneRNA Package, we predicted local structures for the SARS-CoV-2 RNA sequence. We implemented an algorithm in order to aggregate these local structures such that it optimizes the Free Energy of the molecule. We used experimental chemical probing data (SHAPE) to produce more accurate results. We evaluated our predictions using the experimental SHAPE reactivity data and achieved up to $95\%$ match between predicted unpaired positions and highly reactive nucleotides.

The function \texttt{concatinate_local_fold} in \texttt{localfolds.py} returns a file containting a header name, the input RNA sequence and the predicted overall structure and its global free energy.
