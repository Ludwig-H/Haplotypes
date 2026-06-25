# Rapport : dynamique de clusters geometrique sur Z pour graphes signes d'haplotypage

## 1. Objectif

On considere un graphe signe pondere dont les sommets sont des reads ordonnes le long d'un chromosome :

```math
0,1,\ldots,R-1.
```

Chaque read porte un spin cache :

```math
\sigma_i \in \{-1,+1\}.
```

Les poids d'aretes encodent des contraintes ferromagnetiques ou antiferromagnetiques :

- si `W_ij > 0`, l'arete prefere `sigma_i = sigma_j` ;
- si `W_ij < 0`, l'arete prefere `sigma_i != sigma_j` ;
- l'intensite de la contrainte est `|W_ij|`.

L'objectif est de construire une dynamique MCMC adaptee a la geometrie unidimensionnelle, plus globale qu'un Metropolis-Hastings single-spin mais avec un cout comparable lorsque le graphe est local.

Le mouvement propose est centre sur un read `r` choisi uniformement. A partir de ce read, on autorise :

1. ne rien changer ;
2. flipper seulement `r` ;
3. flipper le prefixe `0,...,r-1` ;
4. flipper le prefixe `0,...,r` ;
5. flipper le prefixe `0,...,r+1`.

L'idee est naturelle pour l'haplotypage : un flip de prefixe correspond a deplacer une frontiere de phase, c'est-a-dire un switch local dans la reconstruction le long de l'axe genomique.

## 2. Mesure cible

On ecrit l'energie sous la forme "aretes non satisfaites", comme dans le chapitre 11 :

```math
U(\sigma)
=
\sum_{\{i,j\}: W_{ij}>0}
|W_{ij}| \mathbf 1_{\sigma_i \ne \sigma_j}
+
\sum_{\{i,j\}: W_{ij}<0}
|W_{ij}| \mathbf 1_{\sigma_i = \sigma_j}.
```

La posteriore cible est :

```math
\mu(\sigma \mid W)
\propto
\exp(-U(\sigma)).
```

Si `W_ij` est defini comme log-likelihood ratio pairwise,

```math
W_{ij}
=
\log
\frac{
P(\mathrm{obs}_{ij}\mid \sigma_i\sigma_j=+1)
}{
P(\mathrm{obs}_{ij}\mid \sigma_i\sigma_j=-1)
},
```

alors cette energie est equivalente, a constante pres, a :

```math
\mu(\sigma\mid W)
\propto
\exp\left(
\frac12
\sum_{\{i,j\}} W_{ij}\sigma_i\sigma_j
\right).
```

La convention de ce rapport est donc : `W_ij` est le poids signe naturel de l'arete, et l'energie penalise les contraintes non satisfaites par `|W_ij|`.

## 3. Formule generale de variation d'energie

Pour une arete `e={i,j}`, posons :

```math
y_{ij}=W_{ij}\sigma_i\sigma_j.
```

Avec la convention ci-dessus, une arete est satisfaite si et seulement si :

```math
y_{ij}>0.
```

Soit `A` un ensemble de sommets que l'on flippe :

```math
\sigma_i'=
\begin{cases}
-\sigma_i, & i\in A,\\
\sigma_i, & i\notin A.
\end{cases}
```

Seules les aretes coupees par `A` changent de satisfaction. On note :

```math
\delta(A)=\{\{i,j\}\in E: |\{i,j\}\cap A|=1\}.
```

La variation d'energie est :

```math
\Delta U(A)
=
U(\sigma')-U(\sigma)
=
\sum_{\{i,j\}\in \delta(A)}
W_{ij}\sigma_i\sigma_j
=
\sum_{\{i,j\}\in \delta(A)}
y_{ij}.
```

Cette identite est la cle algorithmique. Pour evaluer un mouvement, il suffit de sommer les contributions des aretes traversees par la coupe induite par ce mouvement.

## 4. Variables duales sur Z

La geometrie sur `Z` suggere d'introduire les variables de murs de domaine :

```math
\tau_t=\sigma_t\sigma_{t+1},
\qquad
t=0,\ldots,R-2.
```

Pour `i<j`, on a :

```math
\sigma_i\sigma_j
=
\prod_{t=i}^{j-1}\tau_t.
```

Un flip de prefixe

```math
P_q=\{0,1,\ldots,q\}
```

ne modifie qu'un seul mur :

```math
\tau_q \mapsto -\tau_q.
```

Un flip singleton `{r}` modifie au plus deux murs :

```math
\tau_{r-1}\mapsto -\tau_{r-1},
\qquad
\tau_r\mapsto -\tau_r,
```

en ignorant les murs hors bornes.

Ainsi, dans les variables duales, les mouvements de prefixe sont locaux. Dans les variables de spins, ils sont des flips macroscopiques de blocs contigus.

## 5. Definition propre des mouvements

On definit les prefixes :

```math
P_q=\{0,\ldots,q\}.
```

Par convention :

```math
P_{-1}=\varnothing,
\qquad
P_{R-1}=\{0,\ldots,R-1\}.
```

Le flip de `P_{R-1}` est un flip global. Il ne change aucune energie pairwise car tous les produits `sigma_i sigma_j` sont invariants.

Pour un read `r`, on propose l'un des cinq mouvements :

```math
\varnothing,
\qquad
\{r\},
\qquad
P_{r-1},
\qquad
P_r,
\qquad
P_{r+1}.
```

Aux bords, les mouvements hors domaine sont rabattus de facon explicite :

- `P_{-1}` est le mouvement nul ;
- `P_R` peut etre remplace par `P_{R-1}`, donc par un flip global.

Cette convention evite les cas particuliers dangereux dans l'implementation.

## 6. Noyau Metropolis-Hastings

Le noyau de proposition est independant de la configuration courante :

1. tirer `r` uniformement dans `{0,...,R-1}` ;
2. tirer un type de mouvement dans les cinq choix ;
3. appliquer Metropolis-Hastings avec :

```math
\alpha(\sigma,A)
=
\min(1,\exp(-\Delta U(A))).
```

Si on travaille a temperature inverse `beta`, on utilise :

```math
\alpha_\beta(\sigma,A)
=
\min(1,\exp(-\beta\Delta U(A))).
```

La proposition est symetrique parce que chaque mouvement est une involution et que sa probabilite ne depend pas de l'etat. Donc le noyau est reversible par rapport a :

```math
\mu_\beta(\sigma)\propto \exp(-\beta U(\sigma)).
```

La chaine est aperiodique grace au mouvement nul. Elle est irreductible car les flips singletons `{r}` sont disponibles pour tous les sommets et engendrent tout l'hypercube `{-1,+1}^R`.

## 7. Encodage optimal des aretes

Pour chaque arete `e={i,j}`, on stocke les indices ordonnes :

```text
left[e]  = min(i,j)
right[e] = max(i,j)
W[e]     = W_ij
y[e]     = W_ij * sigma_i * sigma_j
```

On pre-calcule ensuite deux familles de listes.

### 7.1 Listes incidentes

```text
incident[r] = { e : e est incidente au read r }.
```

Elles servent au flip singleton `{r}` :

```math
\Delta U(\{r\})
=
\sum_{e\in incident[r]} y[e].
```

Si le mouvement est accepte, toutes les aretes incidentes changent de signe :

```math
y[e]\leftarrow -y[e],
\qquad
e\in incident[r].
```

### 7.2 Listes de coupe

Pour une coupe entre `q` et `q+1`, on definit :

```text
cross[q] = { e : left[e] <= q < right[e] }.
```

Elles servent au flip de prefixe `P_q` :

```math
\Delta U(P_q)
=
\sum_{e\in cross[q]} y[e].
```

Si le mouvement est accepte :

```math
y[e]\leftarrow -y[e],
\qquad
e\in cross[q].
```

Pour `q=-1` et `q=R-1`, la coupe est vide : le mouvement est respectivement nul ou global.

## 8. Complexite

Le cout exact d'un pas est :

```math
O(|incident[r]|)
```

pour un flip singleton, et :

```math
O(|cross[q]|)
```

pour un flip de prefixe.

On obtient donc un cout `O(1)` si l'on impose une hypothese de congestion bornee :

```math
\sup_r |incident[r]| = O(1),
\qquad
\sup_q |cross[q]| = O(1).
```

Cette hypothese est realiste pour un graphe d'haplotypage local, avec couverture controlee et longueurs de reads bornees. Elle n'est pas vraie en pire cas :

- reads ultra-longs ;
- region de couverture extreme ;
- graphe trop dense ;
- edges ajoutees a longue portee sans controle.

Il faut donc distinguer deux enonces :

1. algorithme exact en `O(|incident|)` ou `O(|cross|)` ;
2. algorithme attendu en `O(1)` sous hypothese biologique/geometrique de congestion locale bornee.

En pratique, le rapport automatique du benchmark devrait mesurer :

```text
max_degree
mean_degree
max_cut_congestion
mean_cut_congestion
quantiles_cut_congestion
```

Ces quantites disent immediatement si la dynamique sera effectivement quasi constante.

## 9. Construction efficace de `cross`

La construction naive consiste a ajouter chaque arete `e=(i,j)` dans tous les `cross[q]` pour :

```math
i\le q < j.
```

Son cout memoire est :

```math
\sum_{e=(i,j)} (j-i).
```

Ce cout est acceptable si le graphe est local dans l'ordre des reads. Sinon, il peut devenir trop grand.

Plan recommande :

1. commencer avec la representation explicite `cross[q]`, car elle donne le meilleur cout par pas ;
2. mesurer `sum_span = sum_e(right[e]-left[e])` ;
3. si `sum_span` est trop grand, basculer vers une representation alternative `O(log R + |cross[q]|)` par structures d'intervalles, en acceptant que le pas ne soit plus strictement `O(1)`.

Dans un premier prototype mathematique, la representation explicite est la plus saine : elle rend les preuves, les tests et les invariants tres simples.

## 10. Mise a jour des spins

On a deux options.

### Option A : maintenir les spins explicitement

Pour un flip singleton, on fait :

```text
sigma[r] *= -1
```

Pour un flip de prefixe, mettre a jour tous les spins du prefixe couterait `O(R)`, ce qui est a eviter.

### Option B : representation paresseuse par murs

On maintient plutot les variables `tau` et une jauge globale `sigma_0`.

Un flip de prefixe `P_q` change seulement `tau[q]`. Un flip global change seulement `sigma_0`.

Pour reconstruire un spin individuel, il faut :

```math
\sigma_i
=
\sigma_0
\prod_{t=0}^{i-1}\tau_t.
```

Si l'on a besoin de requetes frequentes de spins individuels, on peut utiliser un Fenwick tree de parites pour obtenir `sigma_i` en `O(log R)`. Mais pour la dynamique d'energie, il n'est pas necessaire de reconstruire les spins : les valeurs `y[e]` suffisent.

La recommandation est donc :

- maintenir `y[e]` pour l'energie ;
- maintenir une representation paresseuse des spins seulement pour les sorties, diagnostics et correlations.

## 11. Estimation des correlations spin-spin k-hop

On veut estimer :

```math
C_{ij}
=
\mathbb E_\mu[\sigma_i\sigma_j]
```

pour toutes les paires a distance graphe au plus `k`, par exemple `k=4`.

On definit :

```math
\mathcal P_k
=
\{(i,j): i<j,\ d_G(i,j)\le k\}.
```

Il ne faut pas stocker une matrice dense `R x R` si `R` est grand. On stocke une matrice sparse indexee par les paires de `P_k`.

Pour chaque paire `p=(i,j)`, on maintient :

```text
corr_value[p] = sigma_i * sigma_j
corr_sum[p]
last_time[p]
```

L'estimateur empirique apres `T` pas est :

```math
\widehat C_{ij}
=
\frac1T
\sum_{t=0}^{T-1}
\sigma_i^{(t)}\sigma_j^{(t)}.
```

## 12. Accumulation evenementielle des correlations

Mettre a jour toutes les paires a chaque iteration couterait trop cher. On utilise une accumulation par evenements.

Pour chaque paire `p`, `corr_value[p]` reste constant entre deux flips qui separent ses deux extremites. Si `p` change au temps `t`, on ajoute d'abord sa contribution depuis son dernier changement :

```text
corr_sum[p] += corr_value[p] * (t - last_time[p])
corr_value[p] *= -1
last_time[p] = t
```

A la fin de la chaine, on flush toutes les paires :

```text
corr_sum[p] += corr_value[p] * (T - last_time[p])
C[p] = corr_sum[p] / T
```

Les rejets et les mouvements nuls ne modifient aucune paire. Ils sont automatiquement pris en compte par la duree `t - last_time[p]`.

## 13. Listes de paires affectees

On pre-calcule l'analogue des listes d'aretes.

Pour les flips singletons :

```text
pair_incident[r] = { p=(i,j) in P_k : r=i ou r=j }.
```

Pour les flips de prefixe :

```text
pair_cross[q] = { p=(i,j) in P_k : i <= q < j }.
```

Alors :

- si `{r}` est accepte, seules les paires de `pair_incident[r]` changent de signe ;
- si `P_q` est accepte, seules les paires de `pair_cross[q]` changent de signe.

Le cout de mise a jour des correlations est donc :

```math
O(|pair\_incident[r]|)
```

ou :

```math
O(|pair\_cross[q]|).
```

Sous degre borne et `k` fixe, la taille totale de `P_k` est `O(R)`, et ces mises a jour restent locales en moyenne.

## 14. Construction de `P_k`

On construit `P_k` par BFS tronque depuis chaque sommet :

1. pour chaque sommet `i`, lancer une BFS jusqu'a profondeur `k` dans le graphe non oriente ;
2. pour chaque sommet atteint `j>i`, ajouter `(i,j)` a `P_k` ;
3. stocker la distance `d_G(i,j)` si l'on veut stratifier les correlations par distance.

Complexite :

```math
O\left(R \cdot d^k\right)
```

si le degre typique est `d`.

Pour `k=4`, c'est acceptable seulement si le graphe est sparse. Il faut donc mesurer :

```text
number_tracked_pairs
tracked_pairs_per_node
tracked_pairs_cross_congestion
```

## 15. Rapport entre cette dynamique et Swendsen-Wang

Cette dynamique n'est pas Swendsen-Wang au sens strict :

- elle ne gele pas aleatoirement des aretes satisfaites ;
- elle ne recolorie pas des composantes gelees ;
- elle ne produit pas directement un couplage de percolation comme dans le chapitre 11.

Elle est cependant "cluster-like" dans les spins, car elle propose des flips de blocs contigus potentiellement grands.

Son interpretation naturelle est :

- Metropolis-Hastings local dans les variables de murs `tau` ;
- dynamique de switchs de phase sur `Z` ;
- proposition geometrique adaptee aux erreurs de phase en haplotypage.

Cette distinction est importante. Les preuves de stationnarite sont celles de Metropolis-Hastings, pas celles d'Edwards-Sokal ou de Swendsen-Wang.

## 16. Diagnostics indispensables

Pour valider la dynamique, il faut suivre separement :

```text
acceptance_singleton
acceptance_prefix_r_minus_1
acceptance_prefix_r
acceptance_prefix_r_plus_1
acceptance_global
mean_delta_U_by_move
autocorrelation_energy
autocorrelation_domain_walls
effective_sample_size_correlations
```

Il faut aussi profiler :

```text
mean_incident_size
max_incident_size
mean_cross_size
max_cross_size
mean_pair_cross_size
max_pair_cross_size
```

Ces grandeurs disent si le regime observe est vraiment quasi `O(1)`.

## 17. Tests mathematiques minimaux

### Test 1 : variation d'energie

Sur des petits graphes aleatoires, comparer :

```text
DeltaU_fast(A)
```

avec :

```text
U(flip_A(sigma)) - U(sigma)
```

calcule brutalement.

Ce test doit couvrir :

- flips singletons ;
- prefixes internes ;
- prefixe vide ;
- flip global ;
- poids positifs et negatifs ;
- aretes de longue portee.

### Test 2 : invariance globale

Verifier :

```math
U(\sigma)=U(-\sigma).
```

et :

```math
y_{ij}(\sigma)=y_{ij}(-\sigma).
```

### Test 3 : stationnarite sur petit R

Pour `R <= 16`, enumerer les `2^R` configurations et calculer exactement :

```math
\mu(\sigma)
=
\frac{\exp(-U(\sigma))}{Z}.
```

Puis comparer les frequences MCMC aux probabilites exactes.

### Test 4 : detailed balance

Pour des paires de configurations reliees par un mouvement autorise, verifier numeriquement :

```math
\mu(\sigma)K(\sigma,\sigma')
=
\mu(\sigma')K(\sigma',\sigma).
```

### Test 5 : correlations

Comparer l'accumulation evenementielle avec une accumulation naive sur toutes les paires suivies a chaque iteration.

## 18. Plan de route d'implementation

### Phase 1 : specification mathematique

Formaliser dans la documentation du projet :

- energie cible ;
- convention des poids ;
- mouvements autorises ;
- formule de variation d'energie ;
- preuve de reversibilite MH ;
- hypothese de congestion locale pour le `O(1)`.

### Phase 2 : indexation du graphe

Implementer un module d'indexation produisant :

```text
edges.left
edges.right
edges.weight
edges.y
incident
cross
```

Ajouter les statistiques de congestion dans le rapport d'instance.

### Phase 3 : sampler MH geometrique

Implementer :

```text
step()
propose_move(r)
delta_energy(move)
accept_or_reject(move)
apply_move(move)
```

L'implementation doit traiter explicitement :

- prefixe vide ;
- prefixe global ;
- doublons de mouvements aux bords ;
- temperature `beta`.

### Phase 4 : correlations k-hop

Implementer :

```text
build_k_hop_pairs(k)
pair_incident
pair_cross
event_update_pairs(move, time)
finalize_correlations(T)
```

Sorties recommandees :

```text
correlations_khop.tsv
correlations_khop.npz
correlation_summary.json
```

### Phase 5 : tests exacts

Avant tout benchmark biologique, passer les tests sur petits graphes enumerables.

Priorite :

1. `DeltaU_fast == DeltaU_bruteforce` ;
2. detailed balance ;
3. correlations event-driven vs naive ;
4. invariance par flip global.

### Phase 6 : benchmarks de complexite

Sur les instances HAPLO-BENCH, mesurer :

```text
time_per_step
time_per_accepted_step
mean_cross_size
max_cross_size
mean_pair_cross_size
max_pair_cross_size
acceptance_by_move_type
```

L'objectif est de verifier empiriquement que le cout est controle par la geometrie 1D.

### Phase 7 : comparaison avec dynamiques existantes

Comparer au minimum :

- Metropolis single-spin ;
- prefix-MH geometrique ;
- Swendsen-Wang signe, si disponible ;
- eventuellement une dynamique hybride alternant prefix-MH et single-spin.

Metriques :

```text
energy_trace
autocorrelation_energy
autocorrelation_switches
ESS per second
quality of spin-spin correlations
```

## 19. Points de vigilance

### 19.1 Le `O(1)` depend de la geometrie effective

Le graphe est plonge dans `Z`, mais cela ne suffit pas. Il faut que les aretes soient locales dans l'ordre choisi. Si des reads ultra-longs connectent beaucoup de regions, `cross[q]` peut devenir grand.

Conclusion : le `O(1)` est une propriete du couple `(graphe, ordre)`, pas seulement de l'algorithme.

### 19.2 Ordre des reads

L'ordre doit etre choisi soigneusement :

- par position centrale ;
- ou par debut de read ;
- ou par coordonnee de molecule si disponible.

Le meilleur ordre est celui qui minimise la congestion des coupes :

```math
\max_q |cross[q]|.
```

Il faut donc rapporter cette congestion pour plusieurs ordres possibles si necessaire.

### 19.3 Jauge globale

La posteriore est invariante par flip global :

```math
\sigma \mapsto -\sigma.
```

Les correlations `sigma_i sigma_j` sont invariantes par cette jauge, donc elles sont bien definies. En revanche, les moyennes `E[sigma_i]` ne le sont pas sans fixation de jauge.

### 19.4 Paires k-hop

La notion `k-hop` depend du graphe utilise :

- graphe complet des aretes informatives ;
- graphe filtre par poids ;
- graphe local ;
- graphe apres seuillage.

Il faut fixer cette convention dans les sorties pour que les correlations soient interpretables.

### 19.5 Frustration

Dans un graphe frustre, les flips de prefixes peuvent etre tres efficaces pour bouger des blocs, mais ils ne suppriment pas les barrieres d'energie. Il faudra verifier si la dynamique melange bien dans les regions fortement frustrees.

## 20. Conclusion

La dynamique proposee est mathematiquement propre si elle est presentee comme une dynamique Metropolis-Hastings geometrique sur `Z`.

Son point fort est la representation duale :

```math
\tau_t=\sigma_t\sigma_{t+1}.
```

Dans cette representation, les flips de prefixes deviennent locaux, ce qui correspond bien aux erreurs de phase en haplotypage.

L'encodage naturel du graphe repose sur deux listes :

```text
incident[r]
cross[q]
```

Elles donnent les variations d'energie exactes et les mises a jour exactes des poids courants `y[e]`.

La complexite d'un pas est :

```math
O(|incident[r]|)
```

ou :

```math
O(|cross[q]|).
```

Elle devient effectivement `O(1)` sous hypothese de congestion locale bornee, hypothese qu'il faut mesurer et documenter dans chaque instance.

Enfin, les correlations spin-spin a distance `k-hop` doivent etre stockees sparse et mises a jour par accumulation evenementielle. Cette strategie donne exactement la moyenne empirique MCMC tout en evitant de parcourir toutes les paires a chaque iteration.
