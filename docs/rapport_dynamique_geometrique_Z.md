# Rapport : dynamique de clusters géométrique sur $\mathbb{Z}$ pour graphes signés d'haplotypage

## 1. Objectif

On considère un graphe signé pondéré dont les sommets sont des reads ordonnés le long d'un chromosome :

```math
0, 1, \dots, R-1
```

Chaque read porte un spin caché :

```math
\sigma_i \in \{-1, +1\}
```

Les poids d'arêtes encodent des contraintes ferromagnétiques ou antiferromagnétiques :

- Si $W_{ij} > 0$, l'arête préfère $\sigma_i = \sigma_j$ ;
- Si $W_{ij} < 0$, l'arête préfère $\sigma_i \neq \sigma_j$ ;
- L'intensité de la contrainte est $|W_{ij}|$.

L'objectif est de construire une dynamique MCMC adaptée à la géométrie unidimensionnelle, plus globale qu'un échantillonneur Metropolis-Hastings local (single-spin) mais avec un coût de calcul comparable lorsque le graphe est local.

Le mouvement proposé est centré sur un read $r$ choisi uniformément. À partir de ce read, on s'autorise à :

1. Ne rien changer ;
2. Flipper (inverser) seulement le spin de $r$ ;
3. Flipper le préfixe $0, \dots, r-1$ ;
4. Flipper le préfixe $0, \dots, r$ ;
5. Flipper le préfixe $0, \dots, r+1$.

Cette approche est très naturelle pour l'haplotypage : le renversement d'un bloc de spins (le préfixe) correspond à déplacer une frontière de phase, modélisant ainsi de façon réaliste les switchs locaux lors de la reconstruction le long de l'axe génomique.

## 2. Mesure cible

On écrit l'énergie sous la forme "arêtes non satisfaites" :

```math
U(\sigma)
=
\sum_{\{i,j\}: W_{ij}>0}
|W_{ij}| \mathbf{1}_{\sigma_i \neq \sigma_j}
+
\sum_{\{i,j\}: W_{ij}<0}
|W_{ij}| \mathbf{1}_{\sigma_i = \sigma_j}
```

La postérieure cible est :

```math
\mu(\sigma \mid W)
\propto
\exp\bigl(-U(\sigma)\bigr)
```

Si $W_{ij}$ est défini comme le rapport de log-vraisemblance (*log-likelihood ratio*) par paire :

```math
W_{ij}
=
\log
\frac{
P(\mathrm{obs}_{ij}\mid \sigma_i\sigma_j=+1)
}{
P(\mathrm{obs}_{ij}\mid \sigma_i\sigma_j=-1)
}
```

Alors cette énergie est équivalente, à une constante additive près, à :

```math
\mu(\sigma\mid W)
\propto
\exp\left(
\frac{1}{2}
\sum_{\{i,j\}} W_{ij}\sigma_i\sigma_j
\right)
```

La convention de ce rapport est donc : $W_{ij}$ est le poids signé naturel de l'arête, et l'énergie pénalise les contraintes non satisfaites par $|W_{ij}|$.

## 3. Formule générale de variation d'énergie

Pour une arête $e=\{i,j\}$, posons :

```math
y_{ij} = W_{ij}\sigma_i\sigma_j
```

Avec la convention ci-dessus, une arête est satisfaite si et seulement si :

```math
y_{ij} > 0
```

Soit $A$ un ensemble de sommets que l'on flippe (renversement de spins) :

```math
\sigma_i'=
\begin{cases}
-\sigma_i, & i\in A,\\
\sigma_i, & i\notin A
\end{cases}
```

Seules les arêtes coupées par $A$ changent de satisfaction. On note la coupe induite par $A$ :

```math
\delta(A) = \{\{i,j\}\in E: |\{i,j\}\cap A| = 1\}
```

La variation d'énergie résultant du flip est alors :

```math
\Delta U(A)
=
U(\sigma')-U(\sigma)
=
\sum_{\{i,j\}\in \delta(A)}
W_{ij}\sigma_i\sigma_j
=
\sum_{\{i,j\}\in \delta(A)}
y_{ij}
```

Cette identité est la clé algorithmique. Pour évaluer un mouvement, il suffit de sommer les contributions $y_{ij}$ des arêtes traversées par la coupe $\delta(A)$ induite par ce mouvement.

## 4. Variables duales sur $\mathbb{Z}$

La géométrie sur $\mathbb{Z}$ suggère d'introduire les variables de murs de domaine :

```math
\tau_t = \sigma_t\sigma_{t+1},
\qquad
t = 0, \dots, R-2
```

Pour $i < j$, on reconstruit l'interaction par :

```math
\sigma_i\sigma_j
=
\prod_{t=i}^{j-1}\tau_t
```

Un flip de préfixe :

```math
P_q = \{0, 1, \dots, q\}
```

ne modifie qu'un seul mur dans la représentation duale :

```math
\tau_q \mapsto -\tau_q
```

Un flip singleton $\{r\}$ modifie au plus deux murs :

```math
\tau_{r-1}\mapsto -\tau_{r-1},
\qquad
\tau_r\mapsto -\tau_r
```

en ignorant les murs hors bornes.

Ainsi, dans les variables duales, les mouvements de préfixe sont parfaitement locaux. Dans les variables de spins d'origine, ils correspondent à des flips macroscopiques de blocs contigus le long du chromosome.

## 5. Définition propre des mouvements

On définit les préfixes :

```math
P_q = \{0, \dots, q\}
```

Par convention :

```math
P_{-1} = \varnothing,
\qquad
P_{R-1} = \{0, \dots, R-1\}
```

Le flip de $P_{R-1}$ correspond à un flip global de tous les spins. Il ne change aucune énergie par paire (*pairwise*) car tous les produits $\sigma_i \sigma_j$ restent invariants.

Pour un read $r$, on propose l'un des cinq mouvements d'inversion :

```math
\varnothing,
\qquad
\{r\},
\qquad
P_{r-1},
\qquad
P_r,
\qquad
P_{r+1}
```

Aux bords du domaine, les mouvements hors bornes sont rabattus de façon explicite :

- $P_{-1}$ est le mouvement nul ;
- $P_R$ est remplacé par $P_{R-1}$ (le flip global).

Cette convention évite l'apparition de cas particuliers non triviaux et dangereux dans l'implémentation.

## 6. Noyau Metropolis-Hastings

Le noyau de proposition est indépendant de la configuration de spins actuelle :

1. Tirer un read $r$ uniformément dans $\{0, \dots, R-1\}$ ;
2. Tirer un type de mouvement parmi les cinq choix ;
3. Appliquer le critère d'acceptation Metropolis-Hastings :

```math
\alpha(\sigma,A)
=
\min\bigl(1, \exp(-\Delta U(A))\bigr)
```

Si on travaille à température inverse $\beta$, le taux devient :

```math
\alpha_\beta(\sigma,A)
=
\min\bigl(1, \exp(-\beta\Delta U(A))\bigr)
```

La proposition est symétrique parce que chaque mouvement est une involution et que sa probabilité de sélection ne dépend pas de l'état du système. Le noyau est donc réversible par rapport à la mesure :

```math
\mu_\beta(\sigma) \propto \exp\bigl(-\beta U(\sigma)\bigr)
```

La chaîne est apériodique grâce au mouvement nul. Elle est irréductible car les flips singletons $\{r\}$ sont toujours disponibles pour tous les sommets et engendrent l'ensemble de l'hypercube de configuration $\{-1, +1\}^R$.

## 7. Encodage optimal des arêtes

Pour chaque arête $e=\{i,j\}$, on stocke ses attributs de manière orientée :

```python
left[e]  = min(i,j)
right[e] = max(i,j)
W[e]     = W_ij
y[e]     = W_ij * sigma_i * sigma_j
```

On pré-calcule ensuite deux familles de listes d'arêtes.

### 7.1 Listes incidentes

```math
\text{incident}[r] = \{ e \in E : e \text{ est incidente au read } r \}
```

Elles servent à évaluer le flip singleton $\{r\}$ :

```math
\Delta U(\{r\})
=
\sum_{e \in \text{incident}[r]} y[e]
```

Si le mouvement est accepté, les variables de toutes les arêtes incidentes sont mises à jour :

```math
y[e] \leftarrow -y[e],
\qquad
e \in \text{incident}[r]
```

### 7.2 Listes de coupe

Pour une coupe située entre le read $q$ et $q+1$, on définit :

```math
\text{cross}[q] = \{ e \in E : \text{left}[e] \le q < \text{right}[e] \}
```

Elles servent à évaluer le flip de préfixe $P_q$ :

```math
\Delta U(P_q)
=
\sum_{e \in \text{cross}[q]} y[e]
```

Si le mouvement est accepté :

```math
y[e] \leftarrow -y[e],
\qquad
e \in \text{cross}[q]
```

Pour $q = -1$ et $q = R-1$, la coupe est vide (le coût est nul, ce qui est cohérent avec le mouvement nul ou global).

## 8. Complexité

Le coût exact d'une itération de la dynamique est :

```math
\mathcal{O}(|\text{incident}[r]|)
```

pour un flip singleton, et :

```math
\mathcal{O}(|\text{cross}[q]|)
```

pour un flip de préfixe.

On obtient donc un coût en temps de calcul en $\mathcal{O}(1)$ si l'on impose une hypothèse de congestion bornée :

```math
\sup_r |\text{incident}[r]| = \mathcal{O}(1),
\qquad
\sup_q |\text{cross}[q]| = \mathcal{O}(1)
```

Cette hypothèse est réaliste pour un graphe d'haplotypage local, avec une couverture contrôlée et des longueurs de reads bornées. Elle n'est pas garantie en pire cas si l'on observe :

- Des reads ultra-longs connectant des régions éloignées ;
- Des zones de sur-couverture extrême ;
- Un graphe trop dense en arêtes ;
- Des liaisons ajoutées à longue portée sans contrôle de distance.

Il convient de distinguer deux énoncés :

1. L'algorithme exact s'exécute en $\mathcal{O}(|\text{incident}|)$ ou $\mathcal{O}(|\text{cross}|)$ ;
2. L'algorithme attendu s'exécute en $\mathcal{O}(1)$ sous hypothèse biologique et géométrique de congestion locale bornée.

En pratique, le rapport statistique de benchmark d'instance doit mesurer les variables suivantes :

```yaml
max_degree
mean_degree
max_cut_congestion
mean_cut_congestion
quantiles_cut_congestion
```

Ces indicateurs permettent de valider si la dynamique s'exécutera effectivement à coût constant.

## 9. Construction efficace de `cross`

La construction naïve consiste à ajouter chaque arête $e=(i,j)$ dans toutes les listes de coupe `cross[q]` pour :

```math
i \le q < j
```

Son coût mémoire (en nombre d'entrées stockées) est :

```math
\sum_{e=(i,j)} (j-i)
```

Ce coût est acceptable si le graphe est local dans l'ordre des reads. Sinon, il peut devenir prohibitif.

Plan recommandé :

1. Commencer avec la représentation explicite `cross[q]`, car elle donne le meilleur coût opérationnel par pas ;
2. Mesurer la somme des étendues `sum_span` définie par la somme sur toutes les arêtes de la différence entre `right[e]` et `left[e]` ;
3. Si la variable `sum_span` est trop grande, basculer vers une représentation alternative en $\mathcal{O}(\log R + |\text{cross}[q]|)$ par structures d'intervalles (comme un arbre d'intervalles), en acceptant que le pas ne soit plus strictement en $\mathcal{O}(1)$.

Dans un premier prototype mathématique, la représentation explicite est la plus saine : elle rend les épreuves, les tests et les invariants simples à valider.

## 10. Mise à jour des spins

Deux options sont envisageables pour maintenir l'état des spins.

### Option A : maintenir les spins explicitement

Pour un flip singleton, on effectue la mise à jour directe :

```python
sigma[r] *= -1
```

Pour un flip de préfixe, mettre à jour tous les spins du préfixe individuellement coûterait $\mathcal{O}(R)$, ce qui est incompatible avec un coût de pas constant.

### Option B : représentation paresseuse par murs

On maintient plutôt les variables duales $\tau$ et une jauge globale $\sigma_0$.

Un flip de préfixe $P_q$ change seulement $\tau_q \leftarrow -\tau_q$. Un flip global change seulement $\sigma_0 \leftarrow -\sigma_0$.

Pour reconstruire un spin individuel à la demande, il faut calculer :

```math
\sigma_i
=
\sigma_0
\prod_{t=0}^{i-1}\tau_t
```

Si l'on a besoin de requêtes fréquentes de spins individuels, on peut utiliser un arbre de Fenwick (*Fenwick tree*) de parités pour obtenir n'importe quel $\sigma_i$ en $\mathcal{O}(\log R)$. Cependant, pour le calcul de l'énergie de la dynamique, il n'est pas nécessaire de reconstruire les spins : les valeurs `y[e]` suffisent.

La recommandation est donc :

- Maintenir `y[e]` pour l'évaluation rapide de l'énergie ;
- Maintenir une représentation paresseuse des spins uniquement pour les sorties, les diagnostics et le calcul des corrélations.

## 11. Estimation des corrélations spin-spin $k$-hop

On souhaite estimer l'espérance des corrélations :

```math
C_{ij} = \mathbb{E}_\mu[\sigma_i\sigma_j]
```

pour toutes les paires à distance de graphe au plus $k$, par exemple $k=4$.

On définit l'ensemble des paires suivies :

```math
\mathcal{P}_k = \{(i,j) : i < j,\ d_G(i,j) \le k\}
```

Pour éviter d'utiliser une matrice dense $R \times R$, on stocke les valeurs de manière creuse (*sparse*) sous forme de dictionnaire ou de tableau plat indexé par les éléments de $\mathcal{P}_k$.

Pour chaque paire $p=(i,j)$, on maintient les variables suivantes :

```python
corr_value[p] = sigma_i * sigma_j
corr_sum[p]
last_time[p]
```

L'estimateur empirique après $T$ pas de la chaîne est :

```math
\widehat{C}_{ij}
=
\frac{1}{T}
\sum_{t=0}^{T-1}
\sigma_i^{(t)}\sigma_j^{(t)}
```

## 12. Accumulation événementielle des corrélations

Mettre à jour toutes les paires suivies à chaque itération de la dynamique serait trop coûteux. On utilise à la place une accumulation événementielle.

La variable `corr_value[p]` reste constante entre deux flips acceptés qui séparent les deux extrémités de la paire $p$. Si la valeur de $p$ change au temps $t$, on ajoute d'abord sa contribution accumulée depuis son dernier changement :

```python
corr_sum[p] += corr_value[p] * (t - last_time[p])
corr_value[p] *= -1
last_time[p] = t
```

À la fin de l'échantillonnage, on effectue un flush final pour toutes les paires suivies :

```python
corr_sum[p] += corr_value[p] * (T - last_time[p])
C[p] = corr_sum[p] / T
```

Les rejets et les mouvements nuls ne modifient pas l'état des paires. Ils sont automatiquement pris en compte par la durée $t - \text{last\_time}[p]$ lors de la prochaine modification, où la variable `last_time` désigne le dernier pas de temps mis à jour.

## 13. Listes de paires affectées

On pré-calcule l'analogue des listes d'arêtes pour les paires de corrélation.

Pour les flips singletons, la liste des paires incidentes `pair_incident` est définie par :

```math
\mathcal{P}_{\mathrm{incident}}(r) = \{ p = (i,j) \in \mathcal{P}_k : r=i \text{ ou } r=j \}
```

Pour les flips de préfixe, la liste de coupe de paires `pair_cross` est définie par :

```math
\mathcal{P}_{\mathrm{cross}}(q) = \{ p = (i,j) \in \mathcal{P}_k : i \le q < j \}
```

Alors :

- Si $\{r\}$ est accepté, seules les paires de `pair_incident[r]` changent de signe ;
- Si $P_q$ est accepté, seules les paires de `pair_cross[q]` changent de signe.

Le coût de mise à jour des corrélations est donc proportionnel à la taille des listes affectées :

```math
\mathcal{O}(|\mathcal{P}_{\mathrm{incident}}(r)|) \quad \text{ou} \quad \mathcal{O}(|\mathcal{P}_{\mathrm{cross}}(q)|)
```

Sous hypothèse de degré borné et pour $k$ fixé, la taille totale de $\mathcal{P}_k$ est en $\mathcal{O}(R)$, et ces mises à jour restent locales en moyenne.

## 14. Construction de $\mathcal{P}_k$

On construit $\mathcal{P}_k$ par une recherche en largeur (BFS) tronquée depuis chaque sommet :

1. Pour chaque sommet $i$, lancer une BFS jusqu'à la profondeur $k$ dans le graphe non orienté ;
2. Pour chaque sommet atteint $j > i$, ajouter la paire $(i,j)$ à $\mathcal{P}_k$ ;
3. Stocker la distance de graphe $d_G(i,j)$ si l'on souhaite stratifier les corrélations par distance.

La complexité globale de construction est de :

```math
\mathcal{O}\left(R \cdot d^k\right)
```

où $d$ est le degré moyen du graphe. Pour $k=4$, cette complexité n'est acceptable que si le graphe est creux. Il est donc recommandé d'enregistrer et de surveiller les variables de congestion suivantes :

```yaml
number_tracked_pairs
tracked_pairs_per_node
tracked_pairs_cross_congestion
```

## 15. Rapport entre cette dynamique et Swendsen-Wang

Cette dynamique n'est pas une dynamique de Swendsen-Wang au sens strict :

- Elle ne gèle pas aléatoirement des arêtes satisfaites ;
- Elle ne recolorie pas des composantes gelées ;
- Elle ne produit pas directement un couplage de percolation (comme dans la formulation classique de Swendsen-Wang standard).

Elle est cependant apparentée à une dynamique de type "cluster" pour les spins, car elle propose des flips de blocs contigus potentiellement grands.

Son interprétation naturelle est :

- Un échantillonneur de Metropolis-Hastings local dans les variables de murs $\tau$ ;
- Une dynamique de switchs de phase sur $\mathbb{Z}$ ;
- Une proposition géométrique bien adaptée aux erreurs de phase (switchs de phase) rencontrées en haplotypage.

Cette distinction est importante. Les preuves de stationnarité reposent uniquement sur celles de Metropolis-Hastings, et non sur des arguments de couplage géométrique d'Edwards-Sokal ou de Swendsen-Wang.

## 16. Diagnostics indispensables

Pour valider la dynamique, il faut suivre séparément les indicateurs de taux d'acceptation et de statistiques suivants :

```yaml
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

Il faut également profiler les statistiques de congestion suivantes :

```yaml
mean_incident_size
max_incident_size
mean_cross_size
max_cross_size
mean_pair_cross_size
max_pair_cross_size
```

Ces grandeurs confirment si le régime observé s'approche effectivement du comportement théorique en $\mathcal{O}(1)$.

## 17. Tests mathématiques minimaux

### Test 1 : variation d'énergie

Sur de petits graphes aléatoires, comparer la variable `DeltaU_fast(A)` avec le calcul direct `U(flip_A(sigma)) - U(sigma)` évalué de façon brute, où la variable `DeltaU_bruteforce` sert de comparaison brute. Ce test doit couvrir :

- Les flips singletons ;
- Les préfixes internes ;
- Le préfixe vide ;
- Le flip global ;
- Les poids positifs et négatifs ;
- Les arêtes à longue portée.

### Test 2 : invariance globale (symétrie de jauge)

Vérifier :

```math
U(\sigma) = U(-\sigma)
```

et :

```math
y_{ij}(\sigma) = y_{ij}(-\sigma)
```

### Test 3 : stationnarité sur petit $R$

Pour $R \le 16$, énumérer l'ensemble des $2^R$ configurations possibles de spins, calculer exactement la mesure de Boltzmann :

```math
\mu(\sigma)
=
\frac{\exp\bigl(-U(\sigma)\bigr)}{\mathcal{Z}}
```

et comparer les fréquences empiriques obtenues par MCMC aux probabilités exactes.

### Test 4 : balance détaillée (Detailed Balance)

Pour des paires de configurations reliées par un mouvement autorisé, vérifier numériquement la relation :

```math
\mu(\sigma)K(\sigma, \sigma')
=
\mu(\sigma')K(\sigma', \sigma)
```

### Test 5 : corrélations

Comparer l'accumulation événementielle avec une accumulation naïve effectuant le produit direct sur toutes les paires suivies à chaque itération.

## 18. Plan de route d'implémentation

### Phase 1 : spécification mathématique

Formaliser dans la documentation du projet :

- L'énergie cible ;
- La convention des poids ;
- Les mouvements autorisés ;
- La formule de variation d'énergie ;
- La preuve de réversibilité Metropolis-Hastings (MH) ;
- L'hypothèse de congestion locale pour le coût en $\mathcal{O}(1)$.

### Phase 2 : indexation du graphe

Implémenter un module d'indexation produisant :

```python
edges.left
edges.right
edges.weight
edges.y
incident
cross
```

Ajouter les statistiques de congestion dans le rapport d'instance.

### Phase 3 : échantillonneur MH géométrique

Implémenter les routines suivantes :

```python
step()
propose_move(r)
delta_energy(move)
accept_or_reject(move)
apply_move(move)
```

L'implémentation doit traiter explicitement :

- Le préfixe vide ;
- Le préfixe global ;
- Les doublons de mouvements aux bords ;
- Le paramètre de température inverse $\beta$.

### Phase 4 : corrélations $k$-hop

Implémenter :

```python
build_k_hop_pairs(k)
pair_incident
pair_cross
event_update_pairs(move, time)
finalize_correlations(T)
```

Sorties recommandées :

- Fichier des corrélations : `correlations_khop.tsv`
- Matrice des corrélations : `correlations_khop.npz`
- Résumé des corrélations : `correlation_summary.json`

### Phase 5 : tests exacts

Avant tout benchmark biologique, passer les tests unitaires sur de petits graphes énumérables.

Priorités de validation :

1. Les résultats du test rapide et brut de variation d'énergie concordent ;
2. Balance détaillée (*detailed balance*) ;
3. Corrélations événementielles vs naïves ;
4. Invariance par flip global (symétrie de jauge).

### Phase 6 : benchmarks de complexité

Sur les instances issues de HAPLO-BENCH, mesurer les variables temporelles et statistiques :

- Temps par pas : `time_per_step`
- Temps par pas accepté : `time_per_accepted_step`
- Congestion moyenne des coupes : `mean_cross_size`
- Congestion maximale des coupes : `max_cross_size`
- Congestion moyenne des paires suivies : `mean_pair_cross_size`
- Congestion maximale des paires suivies : `max_pair_cross_size`
- Taux d'acceptation par mouvement : `acceptance_by_move_type`

L'objectif est de vérifier empiriquement que le coût d'une mise à jour est contrôlé par la géométrie 1D.

### Phase 7 : comparaison avec les dynamiques existantes

Comparer au minimum :

- L'échantillonneur de Metropolis classique (single-spin) ;
- La dynamique de préfixe-MH géométrique ;
- La dynamique de Swendsen-Wang signée (si disponible) ;
- Éventuellement, une dynamique hybride alternant les propositions de préfixes et les propositions de single-spin.

Métriques de comparaison :

- Trace d'énergie : `energy_trace`
- Autocorrélation d'énergie : `autocorrelation_energy`
- Autocorrélation des switchs : `autocorrelation_switches`
- ESS par seconde : `ESS_per_second`
- Qualité des corrélations : `quality_of_spin_spin_correlations`

## 19. Points de vigilance

### 19.1 Le $\mathcal{O}(1)$ dépend de la géométrie effective

Le graphe est plongé dans $\mathbb{Z}$, mais cela ne suffit pas. Il faut que les arêtes soient locales dans l'ordre choisi. Si des reads ultra-longs connectent beaucoup de régions disjointes, la liste `cross[q]` peut croître de façon importante.

**Conclusion** : le coût en $\mathcal{O}(1)$ est une propriété conjointe du couple `(graphe, ordre)`, et pas seulement de l'algorithme lui-même.

### 19.2 Ordre des reads

L'ordre d'indexation des reads le long de l'axe linéaire doit être choisi soigneusement :

- Par position centrale de la lecture ;
- Ou par début de lecture ;
- Ou par coordonnée de molécule physique si disponible.

Le meilleur ordre est celui qui minimise la congestion maximale des coupes :

```math
\max_q |\text{cross}[q]|
```

Il convient donc de mesurer et rapporter cette congestion pour plusieurs ordres possibles si nécessaire.

### 19.3 Symétrie de jauge globale

La postérieure est invariante par flip global (symétrie $\mathbb{Z}_2$) :

```math
\sigma \mapsto -\sigma
```

Les corrélations $\sigma_i \sigma_j$ sont invariantes sous cette symétrie et sont donc bien définies. En revanche, les espérances individuelles $\mathbb{E}[\sigma_i]$ sont nulles par symétrie à moins de fixer la jauge.

### 19.4 Paires $k$-hop

La notion de $k$-hop dépend du squelette de graphe utilisé :

- Le graphe complet des arêtes informatives ;
- Le graphe filtré par seuil de poids ;
- Le graphe de connectivité locale.

Il faut fixer rigoureusement cette convention dans les métadonnées de sortie pour rendre les corrélations comparables et interprétables.

### 19.5 Frustration

Dans un graphe fortement frustré, les flips de préfixes peuvent être très efficaces pour déplacer des blocs entiers, mais ils ne suppriment pas les barrières d'énergie locales induites par la frustration. Il convient de vérifier si la dynamique se mélange bien dans les régions de forte frustration.

## 20. Conclusion

La dynamique proposée est mathématiquement rigoureuse si elle est formulée comme une dynamique Metropolis-Hastings géométrique sur $\mathbb{Z}$.

Son point fort réside dans la représentation duale des variables de mur :

```math
\tau_t = \sigma_t\sigma_{t+1}
```

Dans cette représentation, les flips de préfixes deviennent locaux, ce qui correspond physiquement aux erreurs de phase (switchs) rencontrées en haplotypage.

L'encodage du graphe repose naturellement sur deux familles de listes :

* `incident[r]` (les arêtes incidentes à chaque read $r$) ;
* `cross[q]` (les arêtes intersectées par la coupe $q$).

Elles permettent de calculer en temps optimal la variation d'énergie $\Delta U$ et de mettre à jour efficacement les poids d'arêtes `y[e]`.

La complexité d'une itération de la dynamique est en :

```math
\mathcal{O}(|\text{incident}[r]|) \quad \text{ou} \quad \mathcal{O}(|\text{cross}[q]|)
```

Elle devient effectivement en $\mathcal{O}(1)$ sous l'hypothèse de congestion locale bornée, caractéristique géométrique qu'il faut systématiquement mesurer et documenter pour chaque jeu de données.

Enfin, l'estimation des corrélations spin-spin à distance $k$-hop est implémentée de façon creuse (*sparse*) et mise à jour via une accumulation événementielle. Cette stratégie garantit l'obtention exacte des moyennes temporelles MCMC tout en évitant le parcours global coûteux de toutes les paires à chaque itération.
