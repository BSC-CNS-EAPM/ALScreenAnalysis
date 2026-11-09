# ALScreen Framework Pseudocode

## ScoreFormer.Train

```text
Input:  P = {(m_i, y_i)}: training pool, molecules m_i and their docking scores y_i;
        H: hyperparameters (optimizer, learning-rate schedule, batch size,
           weighting strength alpha, log-weight cap c, early-stopping patience p);
        E: maximum number of training epochs.
Output: f_theta: trained ScoreFormer.

 1  G <- {(G_i, y_i)}, with G_i = Graph(m_i), for each (m_i, y_i) in P.
 2  Split G into a training set G_tr and a validation set G_val.
 3  Initialize theta, the optimizer, and the learning-rate schedule from H.
 4  l* <- +inf;  theta* <- theta;  s <- 0.
 5  for e = 1, ..., E do
 6      Set f_theta to training mode.
 7      for each mini-batch B subset of G_tr do
 8          y_hat_i <- f_theta(G_i) for each i in B.
 9          w_i <- exp(min(-alpha * y_i, c)) for each i in B.  // more negative (better) scores get larger weights
10          L_B <- (1 / |B|) * SUM_{i in B} w_i * (y_hat_i - y_i)^2.
11          Update theta using the gradient of L_B and the optimizer.
12      end for
13      Set f_theta to evaluation mode.
14      Compute validation loss and prediction metrics on G_val.  // monitoring only
15      theta^(e) <- copy(theta);  save checkpoint theta^(e).
16      Step the learning-rate schedule.
17  end for
18  return f_theta^(E).
```

## ScoreFormer.Infer

```text
Input:  Q = {m_i}: molecules m_i to be scored;
        f_theta^(E): trained model returned by ScoreFormer.Train.
Output: Y_hat = {(m_i, y_hat_i)}: predicted docking scores keyed by molecule identifier.

 1  Set f_theta^(E) to evaluation mode (dropout off, gradients disabled).
 2  Y_hat <- empty set.
 3  for each mini-batch B subset of Q do
 4      G_B <- {Graph(m_i) : m_i in B}.
 5      y_hat_B <- f_theta^(E)(G_B).
 6      Y_hat <- Y_hat UNION {(m_i, y_hat_i) for m_i in B}.
 7  end for
 8  return Y_hat.
```

## ALScreen (active-learning virtual screening)

```text
Input:  L: compound library;
        O_T: docking oracle for target T;
        n_sub: initial subset size;
        n_oracle: compounds docked per cycle;
        n_pool: compounds added to the training pool per cycle;
        K: maximum number of cycles;
        H, E: training settings;
Output: R: ranked list of prioritized compounds.

 1  S <- RepresentativeSubset(L, n_sub).
 2  D_0 <- O_T(S).  // initial training pool {(m_i, y_i)}
 3  P_0 <- D_0;  D_all <- D_0.  // training pool; all docked compounds
 4  f^(0) <- ScoreFormer.Train(P_0, H, E).  // initial training
 5  k_final <- 0.
 6  for k = 1, ..., K do
 7      Y_hat_k <- ScoreFormer.Infer(L, f^(k-1)).
 8      G_k <- the n_oracle compounds in L with the lowest y_hat_i.
 9      D_k <- O_T(G_k).
10      D_all <- D_all UNION D_k.
11      e_m <- |y_hat_m - y_m| for each m in G_k.
12      C_k <- the n_pool best- and n_pool worst-predicted compounds of G_k (lowest and highest e_m), with their docking labels.
13      P_k <- P_{k-1} UNION C_k.
14      f^(k) <- ScoreFormer.Train(P_k, H, E).  // retrained from scratch
15      k_final <- k.
16  end for
17  Y_hat_final <- ScoreFormer.Infer(L, f^(k_final)).
18  G_final <- the n_oracle compounds in L with the lowest predicted scores in Y_hat_final.
19  D_final <- O_T(G_final).  // final docking round
20  D_all <- D_all UNION D_final.
21  R <- D_all sorted by ascending docking score y_i.
22  return R.
```
