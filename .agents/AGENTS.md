# Instructions de formatage mathématique (LaTeX/KaTeX) pour GitHub Markdown (.md)

Pour garantir que les équations mathématiques insérées dans les fichiers Markdown (.md) soient correctement lues et rendues par le parseur de GitHub (sans interférence avec le formatage Markdown standard comme les italiques ou les échappements), respectez les directives suivantes :

## 1. Isolation des équations de bloc (`$$`)
* Ne placez jamais de texte ou l'équation elle-même sur la même ligne que les délimiteurs `$$`.
* Les délimiteurs `$$` doivent toujours être placés sur des **lignes isolées** pour forcer le parseur Markdown à traiter le contenu comme un bloc brut (Raw Block), désactivant l'analyse des styles Markdown à l'intérieur.
* **Exemple correct :**
  ```markdown
  $$
  U(\sigma) = \sum_{i=1}^N \sigma_i
  $$
  ```

## 2. Accolades littérales en mode mathématique (`{` et `}`)
* N'écrivez pas `\{` et `\}` directement pour rendre des accolades littérales (comme pour définir des ensembles, ex: $\lbrace -1, +1 \rbrace$). Les parseurs Markdown consomment souvent l'antislash d'échappement, laissant KaTeX avec des accolades `{` et `}` orphelines, ce qui lève une erreur de syntaxe.
* Utilisez **exclusivement** les commandes LaTeX textuelles **`\lbrace`** (suivie d'un espace pour ne pas concaténer la commande aux variables qui suivent) et **`\rbrace`**.
* **Exemples corrects :**
  * `$\sigma_i \in \lbrace -1, +1\rbrace$`
  * `$\mathcal{M} = \lbrace \varnothing, \lbrace r\rbrace \rbrace$`
* *Note : Les accolades servant de regroupement structurel ou d'indices (comme `_{ij}` ou `_{k=0}`) ne doivent pas être précédées d'un antislash et restent standard.*

## 3. Retours à la ligne dans les environnements de cas ou matrices (`\\`)
* Pour effectuer un saut de ligne dans un bloc LaTeX de type `cases` ou `align`, le double antislash `\\` classique est réduit à un simple antislash `\` par le parseur Markdown, empêchant le saut de ligne de s'effectuer.
* Utilisez **quatre antislashs (`\\\\`)** à la place du double antislash habituel.
* **Exemple correct :**
  ```markdown
  $$
  \sigma_i' = \begin{cases} -\sigma_i, & i \in A \\\\ \sigma_i, & i \notin A \end{cases}
  $$
  ```
