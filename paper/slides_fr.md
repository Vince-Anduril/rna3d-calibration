# Slides — Présentation 5 min (FR)

> Format : Markdown convertible via `marp` (https://marp.app) ou `pandoc -t beamer`.
> 7 slides, ~40 s par slide.
> Public : projet de cours NLP / transformers — exposé court.

---

## Slide 1 — Titre

**Un transformer d'ARN « voit-il » la fonction sans qu'on la lui apprenne ?**

*Sonder un Language Model d'ARN avec une variante humain/chimpanzé qui change l'activité d'un facteur 4*

Vincent Le Duigou — Albert School Madrid
Projet de cours NLP / Transformers
2026-05-19

---

## Slide 2 — Le contexte (l'ARN comme « langage »)

- Une séquence d'ARN = une phrase de 4 lettres (A, U, G, C).
- Sa fonction biologique dépend de sa **structure 3D**, qui dépend de la **séquence**.
- Depuis 2024, les meilleurs prédicteurs 3D d'ARN sont des **transformers** :
  - **RhoFold+** (*Nature Methods 2024*) — type AlphaFold (Evoformer).
  - **DRfold2** (*PLOS Biol 2025*) — basé sur un **Language Model d'ARN** (RCLM).
- Question classique en NLP : *« le LM a-t-il appris la grammaire, ou juste la statistique de surface ? »*.
- Question équivalente en bio : *« le LM d'ARN a-t-il appris le mécanisme, ou juste la forme moyenne ? »*

---

## Slide 3 — Notre sonde : CPEB3 humain vs chimpanzé

**Une paire unique dans le benchmark Stanford RNA 3D Folding (Kaggle 2025) :**

| | Séquence (69 nt) | Position 30 | Activité |
|---|---|---|---|
| **R1107 (humain)** | …UCGC**A**GCCC… | **A** | 1× |
| **R1108 (chimp)** | …UCGC**G**GCCC… | **G** | **4×** plus rapide |

- **Une seule mutation, multiplicateur d'activité = 4.**
- Skilandat *et al.* (RNA 2016) : la mutation casse un appariement P1/P1.1 dans le bras P1 (résidus 9 et 60).
- **Le LM, qui n'a jamais vu ni cinétique ni biochimie, place-t-il la conséquence structurale au bon endroit ?**

---

## Slide 4 — Le résultat : DRfold2 retrouve l'ancrage P1

> Méthode : on prédit R1107 et R1108 avec DRfold2, on aligne par Kabsch, on regarde où l'écart est le plus grand.

**Top-5 résidus les plus divergents :**
```
{9, 22, 24, 51, 60}
```
Les résidus **9** et **60** = **les deux extrémités du bras P1** — exactement la région que Skilandat 2016 a identifiée biochimiquement.

Le pic n'est PAS au site de la mutation (pos 30) : il est à 21 et 30 résidus de distance, sur les contacts long-range.

> *Le LM, sans supervision fonctionnelle, propage la mutation vers les positions mécanistiquement pertinentes.*

---

## Slide 5 — Le contrôle : c'est spécifique, pas un prior

5 mutations arbitraires (single-nt) ailleurs dans la séquence :

| Mutation | Top-5 divergent | P1 anchor (pos 9 ou 60) ? |
|---|---|---|
| pos **30** (réelle, biologique) | {**9**, 22, 24, 51, **60**} | **✓** |
| pos 5 | {1, 2, 3, 22, 23} | ✗ |
| pos 20 | {23, 24, 26, 48, 50} | ✗ |
| pos 41 | {48, 49, 50, 51, 52} | ✗ |
| pos 55 | {22, 23, 24, 25, 51} | ✗ |
| pos 64 | {50, 64, 65, 66, 67} | ✗ |

**0/5 mutations arbitraires** font apparaître P1. L'ancrage est **spécifique** à la mutation biologique.

---

## Slide 6 — Pourquoi c'est pertinent en NLP / transformers

**1. Le sondage (probing) :**
Comme on demande à BERT « quelles fonctions syntaxiques as-tu apprises sans qu'on te les dise ? », on a demandé à DRfold2 « quels mécanismes as-tu appris sans qu'on te les dise ? ».

**2. La connaissance émergente d'un LM auto-supervisé :**
DRfold2 est entraîné sur **des séquences brutes** d'ARN, sans labels fonctionnels. Il apprend implicitement les règles d'appariement, les motifs structuraux, et — comme on le montre — au moins une *signature mécanistique* (la cascade vers P1).

**3. L'importance de l'architecture :**
RhoFold+ (Evoformer, conçu pour MSA) en mode séquence-seule **ignore l'input** : top-5 identique pour toute mutation. → Les LM « purs » et les transformers MSA-based ne réagissent **pas pareil** aux perturbations.

**4. Analogie NLP :**
- DRfold2 = comme un BERT qui change correctement la prédiction quand on remplace "the *cat* sat" → "the *cats* sat".
- RhoFold+ sans MSA = comme un modèle qui répond la même phrase quelle que soit l'entrée.

---

## Slide 7 — Conclusion et perspectives

**Ce qu'on a montré :**
- Un **language model d'ARN** retrouve, sans supervision fonctionnelle, l'ancrage mécanistique d'une variante humain/chimpanzé documentée biochimiquement.
- Le signal est **spécifique** (0/5 mutations contrôles le reproduisent).
- Tous les LM d'ARN ne réagissent pas pareil — il y a un **angle de calibration différentielle** à explorer.

**Pourquoi ça compte au-delà de la bio :**
Le protocole de probing par micro-perturbation s'applique à **tous les language models de séquences** (NLP comme bio). Une seule lettre changée révèle si le modèle « comprend » ou « moyenne ».

**Coût total expérimental :** ~30 min de GPU RTX 5090 (~0.30 €).
**Code et résultats :** https://github.com/Vince-Anduril/rna3d-calibration

> *« One letter, four-fold function : a transformer probe of mechanistic knowledge in RNA language models. »*

**Merci.**
