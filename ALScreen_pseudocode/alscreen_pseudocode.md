# ALScreen framework Pseudocode

## ScoreFormer.Train

```text
Input: D: molecules and their docking scores; H: model and parameters, including weighting strength alpha and log-weight cap c; E: maximum training epochs.
Output: f_theta: trained ScoreFormer.

 1  Prepare molecular representations G from D for ScoreFormer.
 2  Split G into training set G_tr and validation set G_val.
 3  Initialize f_0 as a ScoreFormer model with random weights.
 4  f <- f_0. Initialize the optimizer and learning-rate schedule using H.
 5  for e = 1, ..., E do
 6      Set f to training mode.
 7      for each mini-batch B subset of G_tr do
 8          y_hat_i <- f(G_i) for each molecule i in B.
 9          w_i <- exp(min(-alpha * y_i, c)).
10          Loss_B <- SUM(w_i * (y_hat_i - y_i)^2 for i in B) / size(B).
11          Update the weights of f using the gradient of Loss_B and the optimizer.
12      end for
13      Set f to evaluation mode.
14      Compute validation loss and prediction metrics on G_val.
15      Save a checkpoint according to the configured policy.
16      Update the learning-rate schedule.
17      if a configured stopping condition is met then break.
18  end for
19  f_theta <- f with its learned weights theta.
20  return f_theta.
```

## ScoreFormer.Infer

```text
Input: L: compounds to be scored; f_theta: trained model returned by ScoreFormer.Train.
Output: Y_hat: predicted docking scores associated with compound identifiers.

 1  Set f_theta to evaluation mode; disable gradients and dropout.
 2  Y_hat <- empty collection.
 3  for each molecular mini-batch B in L do
 4      G_B <- prepare ScoreFormer molecular representations for B.
 5      y_hat_B <- f_theta(G_B).
 6      Append compound identifiers and y_hat_B to Y_hat.
 7  end for
 8  return Y_hat.
```

## ALScreen: active-learning framework

```text
Input: L: compound library; O_T: docking oracle for target T; n_sub: initial subset size; n_oracle: compounds docked per cycle; n_pool: labeled compounds added per cycle; K: maximum cycles; H, E: training settings.
Output: R: ranked list of prioritized compounds.

 1  S <- RepresentativeSubset(L, n_sub).
 2  D <- O_T(S).
 3  P <- D.  // Initialize labeled training pool.
 4  f_theta <- ScoreFormer.Train(P, H, E).  // ScoreFormer initial training.
 5  U <- L \ S.  // Compounds not yet docked.
 6  for k = 1, ..., K do
 7      if U is empty then break.
 8      Y_hat <- ScoreFormer.Infer(U, f_theta).
 9      G <- TopRanked(Y_hat, n_oracle).  // Lowest predicted scores.
10      D_G <- O_T(G).
11      e_i <- |y_hat_i - y_i| for each compound i in D_G.
12      C <- select n_pool labeled records from D_G, split equally
13          between the lowest and highest absolute prediction errors.
14      P <- P UNION C.  // Retain compound identities and docking labels.
15      U <- U \ G.
16      f_theta <- ScoreFormer.Train(P, H, E).  // ScoreFormer retraining.
17      if the configured scaffold-convergence criterion is met then break.
18  end for
19  R <- RankAscending(ScoreFormer.Infer(L, f_theta)).
20  return R.
```
