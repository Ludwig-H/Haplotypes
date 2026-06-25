# HAPLO-BENCH

Benchmark synthétique pour l’inférence d’haplotypes sous modèle bayésien signé.

L’objectif est de générer simultanément :

1. des fichiers biologiques standards ou quasi standards ;
2. un graphe signé pondéré exact, calculé à partir des observations et des paramètres du modèle ;
3. une vérité terrain réservée à l’évaluation.

Le benchmark reste volontairement simple : il isole le problème statistique d’inférence d’haplotypes sans reproduire tout le pipeline biologique réel.

---

## 1. Principe général

On considère une paire de chromosomes homologues, discrétisés sur

$$
\\{0,\\ldots,L-1\\}.
$$

Chaque lecture $i$ est définie par :

* une position centrale $x_i$ ;
* une longueur fixée $2\ell+1$ dans le modèle de base ;
* un haplotype d’origine caché

$$
\\Sigma_i\in\\{-1,+1\\}.
$$

La lecture couvre l’intervalle

$$
I_i=[x_i-\ell,x_i+\ell]\cap \\{0,\\ldots,L-1\\}.
$$

L’objectif est de reconstruire les $\Sigma_i$ à partir des lectures bruitées.

---

## 2. Sorties du benchmark

Une instance du benchmark contient :

```text
instance/
    config.yaml
    truth/
        haplotypes.tsv
        read_truth.tsv
        variants.tsv
        read_assignments.tsv
    bio/
        reference.fa
        reads.fq
        reads_1.fq
        reads_2.fq
        alignments.sam
        alignments.bam
        variants.vcf
    graph/
        nodes.tsv
        edges.tsv
        graph.json
        graph.npz
    report/
        summary.json
        summary.md
```

### Fichiers de vérité terrain

Les fichiers dans `truth/` ne doivent pas être utilisés par les algorithmes.

* `haplotypes.tsv` : séquences ou encodages des deux haplotypes.
* `read_truth.tsv` : haplotype d’origine de chaque lecture.
* `variants.tsv` : positions hétérozygotes vraies et allèles associés.
* `read_assignments.tsv` : vecteur ligne contenant les affectations des reads aux communautés (entiers de 0 à K-1, séparés par des tabulations).

### Fichiers biologiques

Les fichiers dans `bio/` permettent l’interopérabilité.

* `reference.fa` : référence utilisée pour les alignements.
* `reads.fq` : lectures single-end.
* `reads_1.fq`, `reads_2.fq` : lectures paired-end.
* `alignments.sam`, `alignments.bam` : alignements exacts ou simulés.
* `variants.vcf` : sites hétérozygotes simulés.

### Graphe signé pondéré

Le fichier principal est :

```text
graph/edges.tsv
```

Il contient le graphe signé pondéré exact induit par le modèle.

---

## 3. Paramètres principaux

```yaml
reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

chromosome:
  boundary: truncate       # truncate | circular

variants:
  model: bernoulli         # bernoulli | renewal | fixed_count
  density: 0.00075
  min_spacing: 0
  type: biallelic_snp_only

reads:
  position_model: uniform  # uniform | poisson_coverage
  haplotype_prior: balanced
  haplotype_probability_plus: 0.5

coverage:
  target: 30
  compute_n_reads: true

noise:
  model: iid_flip
  epsilon: 0.02
  quality_scores: constant

graph:
  min_shared_variants: 1
  edge_rule: likelihood_ratio
  tie_rule: drop           # drop | random | zero_weight

output:
  write_truth: true
  write_graph: true
  write_fastq: true
  write_bam: true
  write_vcf: true
```

Le nombre de lectures ne doit pas être fixé à la main par défaut. Il est calculé à partir de la couverture cible.

Si les lectures ont longueur moyenne $\mathbb E[R]$, alors

$$
n_{\mathrm{reads}} =
\left\lceil
\frac{cL}{\mathbb E[R]}
\right\rceil.
$$

Dans le cas paired-end, si une molécule produit deux lectures de longueur $R$, alors

$$
n_{\mathrm{pairs}} =
\left\lceil
\frac{cL}{2R}
\right\rceil.
$$

Dans le modèle à longueur fixée $R=2\ell+1$,

$$
n_{\mathrm{reads}} =
\left\lceil
\frac{cL}{2\ell+1}
\right\rceil.
$$

---

## 4. Chromosomes humains recommandés

Deux chromosomes sont recommandés pour les premières expériences.

```yaml
chromosome_presets:
  chr1:
    assembly: GRCh38.p14
    L: 248956422

  chr2:
    assembly: GRCh38.p14
    L: 242193529

  chr20:
    assembly: GRCh38.p14
    L: 64444167

  chr22:
    assembly: GRCh38.p14
    L: 50818468
```

Le chromosome 22 est le meilleur choix pour les tests initiaux : il est suffisamment grand pour voir une vraie limite thermodynamique, mais moins pénible à manipuler que les chromosomes plus longs. Rare moment où la biologie consent à être pratique.

---

## 5. Modèle de variants hétérozygotes

Le benchmark utilise par défaut uniquement des SNPs bialléliques hétérozygotes.

Chaque position $(z\in\\{0,\\ldots,L-1\\})$ est hétérozygote avec probabilité

$$
\rho_{\mathrm{het}}.
$$

Dans le modèle Bernoulli homogène :

$$
H_z\sim \text{Bernoulli}(\rho_{\mathrm{het}}).
$$

L’ensemble des sites hétérozygotes est

$$
\\mathcal V=\\{z:\\ H_z=1\\}.
$$

Le nombre de variants hétérozygotes vérifie

$$
|\mathcal V|\sim \text{Bin}(L,\rho_{\mathrm{het}}),
$$

et donc

$$
\mathbb E[|\mathcal V|]=L\rho_{\mathrm{het}}.
$$

### Densités recommandées

```yaml
variant_presets:
  human_sparse:
    density: 0.0005

  human_default:
    density: 0.00075

  human_high:
    density: 0.001

  toy_dense:
    density: 0.002
```

Interprétation approximative :

| Preset          | Densité $\rho_{\mathrm{het}}$ | Distance moyenne |
| --------------- | -----------------------: | ---------------: |
| `human_sparse`  |         $5\cdot 10^{-4}$ | 1 site / 2000 bp |
| `human_default` |       $7.5\cdot 10^{-4}$ | 1 site / 1333 bp |
| `human_high`    |                $10^{-3}$ | 1 site / 1000 bp |
| `toy_dense`     |         $2\cdot 10^{-3}$ |  1 site / 500 bp |

Pour chr22, avec $(L=50,818,468)$, cela donne en moyenne :

| Preset          | Variants attendus |
| --------------- | ----------------: |
| `human_sparse`  |  $(\approx 25,409)$ |
| `human_default` |  $(\approx 38,114)$ |
| `human_high`    |  $(\approx 50,818)$ |
| `toy_dense`     | $(\approx 101,637)$ |

Le preset `toy_dense` est utile pour les petites instances et les tests algorithmiques, mais il est trop informatif pour un modèle humain réaliste.

---

## 6. Génération des deux haplotypes

À chaque site hétérozygote $z\in\mathcal V$, les deux haplotypes portent deux allèles différents.

On encode l’orientation du site par

$$
A_z\in\\{-1,+1\\}.
$$

Pour une lecture issue de l’haplotype $\Sigma_i$, l’allèle vrai au site $z$ est

$$
X_{iz}^{\mathrm{true}}=\Sigma_i A_z.
$$

Ainsi, pour deux lectures $(i,j)$, au même site hétérozygote $z$, on a

$$
X_{iz}^{\mathrm{true}}X_{jz}^{\mathrm{true}} =
\Sigma_i\Sigma_j.
$$

L’orientation $A_z$ disparaît donc dans les comparaisons entre lectures.

---

## 7. Génération des lectures

Dans le modèle de base, chaque lecture a une position centrale $x_i$ tirée uniformément :

$$
x_i\sim \text{Unif}\\{0,\\ldots,L-1\\}.
$$

La lecture couvre

$$
I_i=[x_i-\ell,x_i+\ell]\cap \\{0,\\ldots,L-1\\}.
$$

L’haplotype d’origine est tiré selon

$$
\mathbb P(\Sigma_i=+1)=\pi,
\qquad
\mathbb P(\Sigma_i=-1)=1-\pi.
$$

Par défaut :

$$
\pi=\frac12.
$$

Configuration minimale :

```yaml
reads:
  molecule_model: single_interval
  length_model: fixed
  half_length: 100
  haplotype_prior: balanced
```

---

## 8. Modèle d’observation bruitée

Pour une lecture $i$ couvrant un site hétérozygote $z$, l’allèle observé est

$$
Y_{iz} = X_{iz}^{\mathrm{true}}\xi_{iz} = \Sigma_i A_z\xi_{iz},
$$

où

$$
\xi_{iz}=
\begin{cases}
+1 & \text{avec probabilité } 1-\varepsilon_{iz} \\
-1 & \text{avec probabilité } \varepsilon_{iz}
\end{cases}
$$

Dans le modèle homogène :

$$
\varepsilon_{iz}=\varepsilon.
$$

Si les qualités Phred sont simulées, alors

$$
\varepsilon_{iz}=10^{-Q_{iz}/10}.
$$

Les indels peuvent être simulés dans les fichiers biologiques. Pour la construction du graphe, ils sont convertis en :

* observation manquante ;
* ou erreur effective au site hétérozygote.

Le modèle recommandé pour commencer reste :

```yaml
noise:
  model: iid_flip
  epsilon: 0.02
  missing_probability: 0.0
  quality_scores: constant
```

---

## 9. Intersections entre lectures

Pour deux lectures $(i,j)$, on définit leur région d’intersection :

$$
O_{ij}=I_i\cap I_j.
$$

Sa longueur est

$$
r_{ij}=|O_{ij}|.
$$

Dans le modèle à longueur fixée sans effet de bord, si

$$
\Delta_{ij}=|x_i-x_j|,
$$

alors

$$
r_{ij} =
\max(0,2\ell+1-\Delta_{ij}).
$$

Les sites hétérozygotes partagés par les deux lectures sont

$$
S_{ij} =
\mathcal V\cap O_{ij}.
$$

Le nombre de sites hétérozygotes partagés est

$$
m_{ij}=|S_{ij}|.
$$

Conditionnellement à $r_{ij}$, dans le modèle Bernoulli homogène,

$$
m_{ij}\sim \text{Bin}(r_{ij},\rho_{\mathrm{het}}).
$$

Donc

$$
\mathbb E[m_{ij}\mid r_{ij}] =
\rho_{\mathrm{het}}r_{ij}.
$$

Une arête est créée si

$$
m_{ij}\ge m_{\min}.
$$

La probabilité qu’une paire de lectures soit informative, conditionnellement à $r_{ij}$, est

$$
\mathbb P(m_{ij}\ge m_{\min}\mid r_{ij}) =
\sum_{h=m_{\min}}^{r_{ij}}
\binom{r_{ij}}{h}
\rho_{\mathrm{het}}^{h}(1-\rho_{\mathrm{het}})^{r_{ij}-h}.
$$

Pour $r_{ij}\rho_{\mathrm{het}}$ petit, l’approximation de Poisson donne

$$
m_{ij}\approx \text{Poisson}(\rho_{\mathrm{het}}r_{ij}).
$$

---

## 10. Concordances et différences observées

Pour chaque site partagé $z\in S_{ij}$, on calcule le vote

$$
V_{ijz}=Y_{iz}Y_{jz}.
$$

Comme

$$
Y_{iz}Y_{jz} =
\Sigma_i\Sigma_j\xi_{iz}\xi_{jz},
$$

on a, en absence d’erreur,

$$
V_{ijz}=\Sigma_i\Sigma_j.
$$

On définit :

$$
c_{ij} = |\\{z\in S_{ij}: V_{ijz}=+1\\}|,
$$

$$
d_{ij} = |\\{z\in S_{ij}: V_{ijz}=-1\\}|.
$$

Ainsi,

$$
m_{ij}=c_{ij}+d_{ij}.
$$

* $c_{ij}$ est le nombre de concordances observées.
* $d_{ij}$ est le nombre de différences observées.

---

## 11. Poids exact d’une arête : mode VCF connu

Ce mode est recommandé.

On suppose que les sites hétérozygotes utilisés dans le graphe sont connus via `variants.vcf`. Le poids est alors calculé conditionnellement à ces sites.

Pour un site partagé $z$, la probabilité que le vote soit correct est

$$
q_{ijz} = \mathbb P(V_{ijz}=\Sigma_i\Sigma_j) = (1-\varepsilon_{iz})(1-\varepsilon_{jz}) + \varepsilon_{iz}\varepsilon_{jz}.
$$

La probabilité que le vote soit inversé est

$$
1-q_{ijz} =
(1-\varepsilon_{iz})\varepsilon_{jz}
+
\varepsilon_{iz}(1-\varepsilon_{jz}).
$$

Le poids signé de l’arête $(e={i,j})$ est le log-rapport de vraisemblance :

$$
W_{ij} = \log \frac{ \mathbb P(\\{V_{ijz}\}_{z\in S_{ij}}\mid \Sigma_i\Sigma_j=+1) }{ \mathbb P(\\{V_{ijz}\}_{z\in S_{ij}}\mid \Sigma_i\Sigma_j=-1) }.
$$

Comme les sites partagés sont conditionnellement indépendants,

$$
W_{ij} =
\sum_{z\in S_{ij}}
V_{ijz}
\log
\frac{q_{ijz}}{1-q_{ijz}}.
$$

Dans le cas homogène $\varepsilon_{iz}=\varepsilon$, on a

$$
q=(1-\varepsilon)^2+\varepsilon^2,
$$

et donc

$$
W_{ij} =
(c_{ij}-d_{ij})
\log
\frac{q}{1-q}.
$$

Puisque

$$
m_{ij}=c_{ij}+d_{ij},
$$

on peut aussi écrire

$$
W_{ij} =
(2c_{ij}-m_{ij})
\log
\frac{q}{1-q}.
$$

L’arête est attractive si

$$
W_{ij}>0,
$$

répulsive si

$$
W_{ij}<0.
$$

Le signe observé est

$$
\widehat S_{ij}=\text{sign}(W_{ij}).
$$

Le poids absolu est

$$
|W_{ij}|.
$$

### Remarque importante sur $\rho_{\mathrm{het}}$

Dans ce mode, $\rho_{\mathrm{het}}$ contrôle :

* le nombre de sites hétérozygotes ;
* la densité du graphe ;
* le nombre moyen de votes par arête ;
* donc indirectement la distribution des poids.

Mais conditionnellement aux sites hétérozygotes observés $S_{ij}$, le paramètre $\rho_{\mathrm{het}}$ se simplifie dans le rapport de vraisemblance. Ce n’est pas une erreur : c’est exactement ce que l’on veut dans le mode VCF connu.

---

## 12. Poids avec confiance de génotypage

Si chaque site hétérozygote $z$ possède une probabilité de confiance

$$
\gamma_z =
\mathbb P(z\text{ vraiment hétérozygote}),
$$

on peut utiliser le modèle de mélange suivant.

Si $z$ est vraiment hétérozygote, le vote est informatif comme précédemment.

Si $z$ n’est pas hétérozygote, le vote ne porte pas d’information sur $\Sigma_i\Sigma_j$.

On note $R_{ijz}(v)$ la probabilité d’observer le vote $(v\in\{-1,+1\})$ lorsque le site n’est pas hétérozygote. Dans le modèle symétrique simple :

$$
R_{ijz}(+1)=q_{ijz},
\qquad
R_{ijz}(-1)=1-q_{ijz}.
$$

Alors

$$
\mathbb P(V_{ijz}=v\mid \Sigma_i\Sigma_j=+1) =
\gamma_z
\mathbb P_{\mathrm{het}}(v\mid +1)
+
(1-\gamma_z)R_{ijz}(v),
$$

et

$$
\mathbb P(V_{ijz}=v\mid \Sigma_i\Sigma_j=-1) =
\gamma_z
\mathbb P_{\mathrm{het}}(v\mid -1)
+
(1-\gamma_z)R_{ijz}(v).
$$

Le poids devient

$$
W_{ij} =
\sum_{z\in S_{ij}}
\log
\frac{
\mathbb P(V_{ijz}\mid \Sigma_i\Sigma_j=+1)
}{
\mathbb P(V_{ijz}\mid \Sigma_i\Sigma_j=-1)
}.
$$

Par défaut :

```yaml
variants:
  genotype_confidence_mode: hard
  gamma: 1.0
```

Le mode `hard` correspond à $\gamma_z=1$ pour tous les sites exportés dans le VCF.

---

## 13. Poids en mode sans VCF

Ce mode est optionnel et moins recommandé.

On ne conditionne plus sur les sites hétérozygotes connus. On intègre l’incertitude sur l’hétérozygotie à chaque position de l’overlap.

Pour une position $z\in O_{ij}$, on suppose

$$
\mathbb P(z\text{ hétérozygote})=\rho_{\mathrm{het}}.
$$

Si les lectures proviennent du même haplotype, alors le vrai vote est $+1$, que le site soit hétérozygote ou non.

Donc

$$
\mathbb P(V_{ijz}=+1\mid \Sigma_i\Sigma_j=+1)=q_{ijz},
$$

$$
\mathbb P(V_{ijz}=-1\mid \Sigma_i\Sigma_j=+1)=1-q_{ijz}.
$$

Si les lectures proviennent d’haplotypes différents :

* avec probabilité $1-\rho_{\mathrm{het}}$, le site n’est pas hétérozygote et le vrai vote est $+1$ ;
* avec probabilité $\rho_{\mathrm{het}}$, le site est hétérozygote et le vrai vote est (-1).

Ainsi,

$$
\mathbb P(V_{ijz}=+1\mid \Sigma_i\Sigma_j=-1) =
(1-\rho_{\mathrm{het}})q_{ijz}
+
\rho_{\mathrm{het}}(1-q_{ijz}),
$$

$$
\mathbb P(V_{ijz}=-1\mid \Sigma_i\Sigma_j=-1) =
(1-\rho_{\mathrm{het}})(1-q_{ijz})
+
\rho_{\mathrm{het}}q_{ijz}.
$$

Le poids est alors

$$
W_{ij} =
\sum_{z\in O_{ij}}
\log
\frac{
\mathbb P(V_{ijz}\mid \Sigma_i\Sigma_j=+1)
}{
\mathbb P(V_{ijz}\mid \Sigma_i\Sigma_j=-1)
}.
$$

Ce mode donne un rôle explicite à $\rho_{\mathrm{het}}$ dans chaque terme du poids.

Il est toutefois déconseillé comme mode principal, car les sites homozygotes dominent très vite les reads longs. Le mode recommandé pour l’inférence d’haplotypes reste le mode VCF connu ou VCF simulé.

---

## 14. Mesure de Gibbs associée

Le graphe signé pondéré définit une postérieure de type Ising signé :

$$
\mu(\sigma\mid W)
\propto
\exp
\left(
\sum_{\{i,j\}\in E}
W_{ij}\sigma_i\sigma_j
\right).
$$

Équivalemment, l’énergie est

$$
U(\sigma) =
\sum_{\{i,j\}:W_{ij}>0}
|W_{ij}|\mathbf 1_{\sigma_i\neq\sigma_j}
+
\sum_{\{i,j\}:W_{ij}<0}
|W_{ij}|\mathbf 1_{\sigma_i=\sigma_j}.
$$

Cette forme est compatible avec les dynamiques de type Swendsen-Wang signé.

## 15. Formats du graphe

### `nodes.tsv`

```tsv
node_id  read_id  molecule_id  position  start  end  read_type
0        read_0   mol_0        1532      1432   1632 theory_fixed
```

La vérité cachée ne doit pas apparaître dans `nodes.tsv` par défaut.

### `read_truth.tsv`

```tsv
read_id  molecule_id  true_haplotype
read_0   mol_0        +1
```

### `edges.tsv`

```tsv
source  target  overlap_length  shared_het_sites  concordances  differences  weight  sign
0       17      145             4                 3             1            2.71    +1
```

Colonnes recommandées :

* `source`, `target` : identifiants des lectures ;
* `overlap_length` : longueur de l’intersection $r_{ij}$ ;
* `shared_het_sites` : $m_{ij}$ ;
* `concordances` : $c_{ij}$ ;
* `differences` : $d_{ij}$ ;
* `weight` : $W_{ij}$ ;
* `sign` : $\text{sign}(W_{ij})$.

### `graph.json`

```json
{
  "model": "haplo-bench",
  "reference": {
    "assembly": "GRCh38.p14",
    "chromosome": "chr22",
    "L": 50818468
  },
  "variants": {
    "density": 0.00075,
    "model": "bernoulli"
  },
  "reads": {
    "preset": "pacbio_hifi"
  },
  "graph": {
    "n_nodes": 84698,
    "n_edges": 1234567,
    "weighted": true,
    "signed": true
  }
}
```

---

## 16. Formats biologiques

### FASTA

`reference.fa` contient une pseudo-référence, par exemple $H^+$.

### FASTQ

`reads.fq` contient les lectures bruitées.

Exemple d’identifiant :

```text
@read_000123 pos=45201
```

La vérité cachée ne doit pas être présente dans les identifiants publics.

### SAM/BAM

`alignments.sam` ou `alignments.bam` contient l’alignement simulé.

Dans le modèle sans indels, le CIGAR est simplement :

```text
(2l+1)M
```

Dans les presets ONT avec indels, le CIGAR peut contenir :

* `M` ;
* `I` ;
* `D`.

Le BAM est produit par le simulateur. Il n’est pas nécessaire d’appeler un aligneur externe.

### VCF

`variants.vcf` contient les variants hétérozygotes simulés.

Le mode principal du graphe suppose que ces variants sont connus.

---

## 17. Presets de séquençage

Les presets fixent des paramètres effectifs. Ils ne prétendent pas simuler toute la chimie du séquençage, parce qu’apparemment nous souhaitons encore terminer ce benchmark avant la fin du siècle.

---

### 17.1 `theory_fixed`

Profil minimal, parfaitement aligné avec le modèle mathématique.

```yaml
preset: theory_fixed

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 30
  compute_n_reads: true

reads:
  molecule_model: single_interval
  length_model: fixed
  half_length: 100

noise:
  model: iid_flip
  epsilon: 0.02
  indels: false
  quality_scores: constant

graph:
  min_shared_variants: 1
  edge_rule: likelihood_ratio
```

---

### 17.2 `illumina_pe150`

Profil paired-end court-read.

```yaml
preset: illumina_pe150

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 30
  compute_n_reads: true

reads:
  molecule_model: paired_end
  read_length: 151
  insert_size_model: normal
  insert_size_mean: 350
  insert_size_sd: 50

noise:
  model: phred_effective
  mean_q: 30
  epsilon_effective: 0.001
  indels: false

graph:
  node_model: read_end       # read_end | molecule
  add_pair_link: true
  pair_link_error: 0.001
  min_shared_variants: 1
```

Pour une paire correcte, les deux extrémités proviennent du même haplotype. Si elles sont modélisées comme deux nœuds distincts, on ajoute une arête ferromagnétique interne de poids

$$
J_{\mathrm{pair}} =
\log
\frac{1-\delta_{\mathrm{pair}}}{\delta_{\mathrm{pair}}}.
$$

Par défaut :

$$
\delta_{\mathrm{pair}}=0.001.
$$

---

### 17.3 `pacbio_hifi`

Profil long-read précis.

```yaml
preset: pacbio_hifi

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 30
  compute_n_reads: true

reads:
  molecule_model: single_interval
  length_model: lognormal
  read_length_mean: 18000
  read_length_sd: 4000
  min_read_length: 1000
  max_read_length: 40000

noise:
  model: phred_effective
  mean_q: 33
  epsilon_effective: 0.0005
  indels: false

graph:
  min_shared_variants: 1
  edge_rule: likelihood_ratio
```

Avec $\rho_{\mathrm{het}}=0.00075$, une lecture de 18 kb couvre en moyenne

$$
18000\times 0.00075 = 13.5
$$

sites hétérozygotes.

---

### 17.4 `ont_q20`

Profil long-read ONT simplex Q20+.

```yaml
preset: ont_q20

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 30
  compute_n_reads: true

reads:
  molecule_model: single_interval
  length_model: lognormal
  read_length_n50: 20000
  read_length_sd: 12000
  min_read_length: 500
  max_read_length: 150000

noise:
  model: phred_effective_with_indels
  mean_q: 20
  epsilon_effective: 0.01
  indels: true
  missing_probability: 0.01

graph:
  min_shared_variants: 2
  edge_rule: likelihood_ratio
```

Le seuil `min_shared_variants: 2` est recommandé pour éviter qu’un unique variant bruité crée une arête trop instable.

---

### 17.5 `ont_duplex`

Profil ONT duplex plus précis.

```yaml
preset: ont_duplex

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 30
  compute_n_reads: true

reads:
  molecule_model: single_interval
  length_model: lognormal
  read_length_n50: 20000
  read_length_sd: 10000
  min_read_length: 500
  max_read_length: 150000

noise:
  model: phred_effective_with_indels
  mean_q: 30
  epsilon_effective: 0.001
  indels: true
  missing_probability: 0.005

graph:
  min_shared_variants: 1
  edge_rule: likelihood_ratio
```

---

### 17.6 `ont_ultralong`

Profil ONT ultra-long.

```yaml
preset: ont_ultralong

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 30
  compute_n_reads: true

reads:
  molecule_model: single_interval
  length_model: pareto_lognormal
  read_length_n50: 75000
  read_length_tail_alpha: 1.8
  min_read_length: 1000
  max_read_length: 4000000

noise:
  model: phred_effective_with_indels
  mean_q: 20
  epsilon_effective: 0.01
  indels: true
  missing_probability: 0.02

graph:
  min_shared_variants: 3
  edge_rule: likelihood_ratio
```

Ce preset sert surtout à tester l’effet des longues connexions sur la percolation du graphe.

---

### 17.7 `hybrid_illumina_pacbio`

Profil hybride short reads précis + PacBio HiFi.

```yaml
preset: hybrid_illumina_pacbio

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

mixture:
  enabled: true
  components:
    - name: illumina_pe150
      coverage: 20
    - name: pacbio_hifi
      coverage: 10

graph:
  min_shared_variants: 1
  edge_rule: likelihood_ratio
```

---

### 17.8 `hybrid_illumina_ont`

Profil hybride short reads précis + ONT bruité.

```yaml
preset: hybrid_illumina_ont

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

mixture:
  enabled: true
  components:
    - name: illumina_pe150
      coverage: 20
    - name: ont_q20
      coverage: 10

graph:
  min_shared_variants: 2
  edge_rule: likelihood_ratio
```

---

## 18. Modes de graphe

### 18.1 Tous les overlaps informatifs

```yaml
graph:
  mode: all_informative_overlaps
  min_shared_variants: 1
```

Toutes les paires partageant au moins un site hétérozygote sont reliées.

### 18.2 Filtrage par nombre de variants

```yaml
graph:
  mode: min_shared_variants
  min_shared_variants: 2
```

On conserve seulement les arêtes telles que

$$
m_{ij}\ge m_{\min}.
$$

### 18.3 Filtrage par poids

```yaml
graph:
  mode: min_abs_weight
  min_abs_weight: 1.0
```

On conserve seulement les arêtes telles que

$$
|W_{ij}|\ge w_{\min}.
$$

### 18.4 Graphe local

```yaml
graph:
  mode: local_window
  max_center_distance: 5000
```

On conserve seulement les paires dont les centres sont proches.

### 18.5 Graphe mixte local + longue portée

```yaml
graph:
  mode: local_plus_long_range
  max_center_distance_local: 5000
  keep_long_reads: true
  keep_pair_links: true
```

Ce mode est utile pour les profils hybrides et paired-end.

---

## 19. Rapport automatique

Chaque instance doit produire un rapport court.

```yaml
summary:
  preset: pacbio_hifi
  chromosome: chr22
  chromosome_length: 50818468
  target_coverage: 30
  empirical_coverage: 30.1
  variant_density: 0.00075
  n_variants: 38114
  n_reads: 84698
  n_edges: 1234567
  mean_overlap_length: 9000.2
  mean_shared_variants_per_edge: 6.8
  positive_edges_fraction: 0.51
  negative_edges_fraction: 0.49
  mean_abs_weight: 5.4
```

Le rapport doit permettre de détecter immédiatement une instance absurde.

---

## 20. Commandes attendues

Générer une instance :

```bash
haplo-bench generate --config configs/pacbio_hifi_chr22.yaml --out data/pacbio_chr22_001
```

Valider une instance :

```bash
haplo-bench validate data/pacbio_chr22_001
```

Exporter uniquement le graphe :

```bash
haplo-bench export-graph data/pacbio_chr22_001
```

Exporter les fichiers biologiques :

```bash
haplo-bench export-bio data/pacbio_chr22_001
```

---

## 21. Instance minimale recommandée

```yaml
preset: theory_fixed

reference:
  assembly: GRCh38.p14
  chromosome: chr22
  L: 50818468

variants:
  model: bernoulli
  density: 0.00075

coverage:
  target: 10
  compute_n_reads: true

reads:
  molecule_model: single_interval
  length_model: fixed
  half_length: 500

noise:
  model: iid_flip
  epsilon: 0.02

graph:
  min_shared_variants: 1
  edge_rule: likelihood_ratio

Cette instance est volontairement plus facile qu’un short-read Illumina réaliste : elle sert à déboguer les algorithmes.

---

## 22. Principe final

L’utilisateur choisit un preset :

```yaml
preset: pacbio_hifi
```

ou

```yaml
preset: ont_q20
```

puis ajuste seulement :

```yaml
reference:
  chromosome: chr22
  L: 50818468

variants:
  density: 0.00075

coverage:
  target: 30

graph:
  min_shared_variants: 1
```

Le benchmark doit produire :

1. un BAM ou SAM cohérent ;
2. un VCF cohérent ;
3. un graphe signé pondéré exact ;
4. une vérité terrain séparée ;
5. un rapport automatique lisible.

Le graphe est l’objet probabiliste principal. Les fichiers biologiques sont une interface, pas une punition.