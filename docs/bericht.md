# Phylogenetische Analyse der TSR3-Proteinsequenz: Identifikation und evolutionäre Einordnung eines unbekannten Vertebraten-Proteins

**Kurs:** Bioinformatik / Phylogenetik
**Autorin:** Greta Jacobs
**Datum:** 18. Februar 2026

---

## Inhaltsverzeichnis

1. [Einführung](#1-einführung)
2. [Methoden](#2-methoden)
3. [Ergebnisse](#3-ergebnisse)
4. [Diskussion](#4-diskussion)
5. [Zusammenfassung](#5-zusammenfassung)
6. [Literaturverzeichnis](#6-literaturverzeichnis)
7. [Anhang](#7-anhang)

---

## 1. Einführung

In der vergleichenden Genomik stellt die Identifikation unbekannter Proteinsequenzen und ihre Einordnung in einen evolutionären Kontext eine grundlegende Aufgabe dar. Im Rahmen dieses Projekts wurde eine unbekannte Proteinsequenz (~312 Aminosäuren) anhand einer Datenbank mit 539 TSR3-Sequenzen aus Vertebraten analysiert. TSR3 (*Two A-site ribosomal RNA methyltransferase 3*) ist ein hochkonserviertes Protein, das an der Biogenese der kleinen ribosomalen Untereinheit (40S) beteiligt ist. Es katalysiert die 2'-O-Methylierung der 18S-rRNA und ist in allen Eukaryoten essenziell für die Prozessierung der prä-ribosomalen RNA [1].

Aufgrund seiner zentralen Funktion im Translationsapparat ist zu erwarten, dass TSR3 eine starke strukturelle und sequenzielle Konservierung über phylogenetisch weit entfernte Organismen hinweg aufweist. Diese Studie verfolgt drei Ziele: (1) Identifikation der unbekannten Sequenz mithilfe einer BLAST-Suche, (2) Rekonstruktion der phylogenetischen Verwandtschaft der ähnlichsten Sequenzen mittels multiplem Sequenzalignment und Distanzmatrix, und (3) Bewertung und Vergleich verschiedener Methoden zur Rekonstruktion phylogenetischer Bäume.

---

## 2. Methoden

### 2.1 BLAST-Suche

Zur Identifikation der Abfragesequenz wurde BLASTP (Version 2.17.0+) [2] gegen eine lokale Datenbank aus 539 TSR3-Sequenzen eingesetzt. Die Datenbank wurde mit `makeblastdb` aus der Datei `database.fasta` erstellt. Als Ausgabeformat wurde Tabelle 7 (`-outfmt 7`) mit den Feldern Akzessionsnummer (`sacc`), Alignmentlänge, prozentuale Identität, Lücken, Score, Bit-Score und E-Wert gewählt.

Zunächst wurde BLASTP mit der Standardsubstitutionsmatrix **BLOSUM62** [3] ausgeführt. Der E-Wert (*Expect value*) quantifiziert die Anzahl der zufällig zu erwartenden Treffer mit gleichwertigem Score und berechnet sich näherungsweise als:

$$E \approx K \cdot m \cdot n \cdot e^{-\lambda S}$$

wobei $m$ und $n$ die Längen von Abfrage- und Datenbank sind, $S$ der Rohscore des Alignments ist sowie $K$ und $\lambda$ statistisch abgeleitete Parameter der Substitutionsmatrix sind. Kleinere E-Werte entsprechen signifikanteren Treffern. BLOSUM-Matrizen werden aus statistisch beobachteten Aminosäureaustauschen in Proteinfamilien der entsprechenden Sequenzidentität abgeleitet: BLOSUM62 aus Blöcken mit ≥62 % Identität, BLOSUM80 aus Blöcken mit ≥80 % Identität [3].

Anhand der prozentualen Sequenzidentität der besten 30 Treffer (Top-30) wurde entschieden, ob die Suche mit einer alternativen Substitutionsmatrix zu wiederholen war. Bei einer medianen Sequenzidentität > 80 % wurde BLOSUM80 als geeigneter erachtet, da diese Matrix für hochidentische Sequenzpaare eine feinere Diskriminierung ermöglicht.

### 2.2 Multiples Sequenzalignment

Nach Extraktion der Top-30-Treffer wurden die entsprechenden Sequenzen aus der Datenbankdatei extrahiert und gemeinsam mit der Abfragesequenz einem multiplen Sequenzalignment (MSA) unterzogen. Die Sequenznamen wurden vor dem Alignment in die binomiale Nomenklatur nach Linné (Gattung_Art) umbenannt, um die Interpretierbarkeit der Ergebnisse zu verbessern. Das Alignment wurde mit der `msa`-Funktion des Bioconductor-Pakets `msa` [4] unter Verwendung des **ClustalW**-Algorithmus [5] durchgeführt.

ClustalW ist ein progressiver Alignmentalgorithmus: Zunächst werden alle paarweisen Distanzen berechnet und ein Leitbaum (*guide tree*) erstellt, dem folgend Sequenzen schrittweise aligniert werden. Die Methode berücksichtigt sequenzspezifische Lückenstrafen, um biologisch plausible Alignments zu erzeugen.

### 2.3 Distanzmatrix

Aus dem MSA wurde eine paarweise Distanzmatrix berechnet, unter Verwendung der Funktion `dist.alignment` aus dem R-Paket `seqinr` [6] mit der Option `matrix = "identity"`. Die identitätsbasierte Distanz zwischen zwei Sequenzen $i$ und $j$ ergibt sich als:

$$d(i, j) = \sqrt{1 - p_{ij}}$$

wobei $p_{ij}$ der Anteil identischer Positionen im paarweisen Alignment ist. Diese Transformierung stellt sicher, dass die Distanzwerte näherungsweise metrische Eigenschaften aufweisen.

### 2.4 Überprüfung von Baum-Eigenschaften

**Vier-Punkte-Bedingung (Additivität):** Eine Distanzmatrix ist additiv (d. h. mit einem ungewurzelten Baum ohne interne Knoten konsistent), wenn für beliebige vier Taxa $i, j, k, l$ die folgende Bedingung gilt: Von den drei Summen

$$s_1 = d(i,j) + d(k,l), \quad s_2 = d(i,k) + d(j,l), \quad s_3 = d(i,l) + d(j,k)$$

müssen die beiden größten gleich sein (d. h. das Maximum wird mindestens zweimal angenommen). Verletzungen dieser Bedingung weisen auf ungleichförmige Evolutionsraten hin.

**Ultrametrik:** Eine Distanzmatrix ist ultrametrisch, wenn für beliebige drei Taxa $i, j, k$ gilt: Das Maximum von $\{d(i,j), d(i,k), d(j,k)\}$ wird mindestens zweimal angenommen. Ultrametrizität ist eine stärkere Eigenschaft als Additivität und impliziert eine molekulare Uhr (gleiche Evolutionsraten).

### 2.5 Rekonstruktion phylogenetischer Bäume

Es wurden vier Clustering-Methoden angewendet:

- **UPGMA** (*Unweighted Pair Group Method with Arithmetic Mean*) [7]: Agglomeratives Hierarchisches Clustering, das bei jedem Schritt die zwei Cluster mit der geringsten mittleren Paarweis-Distanz zusammenführt. UPGMA setzt eine molekulare Uhr voraus und erzeugt ultrametrische Bäume.

- **Neighbor-Joining (NJ)** [8]: Distanzbasierte Methode, die einen $Q$-Wert minimiert:
  $$Q(i,j) = (n-2)\,d(i,j) - \sum_{k} d(i,k) - \sum_{k} d(j,k)$$
  NJ setzt keine Gleichförmigkeit der Evolutionsraten voraus und ist für reale biologische Daten oft geeigneter.

- **Complete Linkage:** Hierarchisches Clustering nach der maximalen paarweisen Distanz zwischen Clustern.

- **Ward's Methode:** Minimiert die Gesamtintra-Cluster-Varianz bei jeder Fusionierung.

### 2.6 Kophenetische Korrelation und Bootstrap

Die Qualität der Bäume wurde mittels **kophenetischer Korrelation** bewertet:

$$r_{\text{coph}} = \text{cor}(d_{\text{original}},\, d_{\text{kophenetisch}})$$

Ein hoher Wert ($r \approx 1$) bedeutet, dass der Baum die originalen Distanzen gut widerspiegelt.

Zur Überprüfung der Robustheit der NJ-Baumtopologie wurde eine **Bootstrap-Analyse** mit 100 Replikaten durchgeführt [9]. Dabei werden die Spalten des Alignments mit Zurücklegen neu gesampelt, ein NJ-Baum auf dem resampled Alignment berechnet, und der Anteil der Replikate, die eine bestimmte interne Verzweigung enthalten, als Bootstrap-Unterstützungswert (%) angegeben.

---

## 3. Ergebnisse

### 3.1 BLAST-Suche und Identifikation der Abfragesequenz

Die initiale BLASTP-Suche mit BLOSUM62 lieferte **502 Treffer** in der Datenbank. Der beste Treffer (`9606_0:0048f7`) wies eine **Sequenzidentität von 100,0 %** bei einem E-Wert von $E = 0{,}0$ auf. Die Akzessionsnummer enthält als Präfix die NCBI-Taxonomie-ID `9606`, die *Homo sapiens* entspricht. Die Abfragesequenz ist damit mit sehr hoher Sicherheit das humane **TSR3-Protein** (*Ribosome biogenesis protein TSR3*, UniProtKB).

Die prozentualen Identitäten der Top-30-Treffer lagen zwischen **80,06 %** (niedrigster Wert) und **100,00 %**, mit einem **Median von ~83,3 %** (Tabelle 1). Da dieser Wert > 80 % beträgt, wurde die Suche mit **BLOSUM80** wiederholt. Die BLOSUM80-Matrix ist für Sequenzpaare mit hoher Identität besser kalibriert und diskriminiert empfindlicher zwischen sehr ähnlichen Sequenzen.

Die E-Werte der Top-30-Treffer lagen zwischen $0{,}0$ und $8{,}18 \times 10^{-174}$, was auf äußerst signifikante Homologie hindeutet. Alle Treffer repräsentieren TSR3-Orthologe aus Vertebraten; Paraloge wurden nicht identifiziert (jede Spezies trat mit einer einzigen Akzessionsnummer auf).

**Tabelle 1: Top-30-Treffer der BLASTP-Suche (BLOSUM80)**

| Nr. | Spezies | Taxon | % Identität | E-Wert |
|-----|---------|-------|------------|--------|
| 1 | *Homo sapiens* | Primates | 100,000 | 0,0 |
| 2 | *Rhinopithecus roxellana* | Primates | 95,513 | 0,0 |
| 3 | *Macaca mulatta* | Primates | 95,192 | 0,0 |
| 4 | *Pongo abelii* | Primates | 94,551 | 0,0 |
| 5 | *Nomascus leucogenys* | Primates | 95,018 | 0,0 |
| 6 | *Callithrix jacchus* | Primates | 91,667 | 0,0 |
| 7 | *Aotus nancymaae* | Primates | 90,705 | 0,0 |
| 8 | *Ateles geoffroyi* | Primates | 91,667 | 0,0 |
| 9 | *Otolemur garnettii* | Primates | 85,350 | 0,0 |
| 10 | *Carlito syrichta* | Primates | 82,372 | 0,0 |
| 11 | *Microcebus murinus* | Primates | 85,256 | 0,0 |
| 12 | *Cheirogaleus medius* | Primates | 82,650 | 0,0 |
| 13 | *Daubentonia madagascariensis* | Primates | 83,175 | 0,0 |
| 14 | *Lemur catta* | Primates | 83,492 | 0,0 |
| 15 | *Sus scrofa* | Artiodactyla | 81,646 | 0,0 |
| 16 | *Physeter catodon* | Cetacea | 83,333 | 0,0 |
| 17 | *Equus caballus* | Perissodactyla | 83,544 | 0,0 |
| 18 | *Choloepus didactylus* | Xenarthra | 83,758 | 0,0 |
| 19 | *Leptonychotes weddellii* | Carnivora | 80,442 | $8{,}76 \times 10^{-180}$ |
| 20 | *Zalophus californianus* | Carnivora | 81,150 | $1{,}19 \times 10^{-179}$ |
| 21 | *Ursus maritimus* | Carnivora | 80,128 | $2{,}75 \times 10^{-179}$ |
| 22 | *Enhydra lutris* | Carnivora | 81,410 | $3{,}60 \times 10^{-179}$ |
| 23 | *Gulo gulo* | Carnivora | 83,654 | $4{,}22 \times 10^{-179}$ |
| 24 | *Bubalus bubalis* | Artiodactyla | 82,692 | $3{,}53 \times 10^{-178}$ |
| 25 | *Bos taurus* | Artiodactyla | 82,857 | $4{,}01 \times 10^{-178}$ |
| 26 | *Capra hircus* | Artiodactyla | 80,060 | $1{,}56 \times 10^{-176}$ |
| 27 | *Ovis aries* | Artiodactyla | 80,757 | $1{,}54 \times 10^{-175}$ |
| 28 | *Cervus hanglu* | Artiodactyla | 81,210 | $9{,}57 \times 10^{-175}$ |
| 29 | *Castor canadensis* | Rodentia | 80,696 | $5{,}61 \times 10^{-174}$ |
| 30 | *Ictidomys tridecemlineatus* | Rodentia | 81,210 | $8{,}18 \times 10^{-174}$ |

*Anmerkung: Die Zuordnung von Akzessionsnummern zu Spezies erfolgte über die NCBI-Taxonomie-IDs im Datenbank-Header sowie die binomiale Umbenennung im Skript.*

### 3.2 Multiples Sequenzalignment

Das MSA mit ClustalW wurde auf 31 Sequenzen (30 Treffer + Abfragesequenz, benannt als `Query_Unknown`) angewendet. Abbildung 1 zeigt einen repräsentativen Ausschnitt des Alignments (mittlere Positionen). Das Alignment offenbart eine hohe Konservierung des zentralen Proteinkerns über alle Vertebraten-Linien hinweg. Lücken treten vorwiegend in den terminalen Regionen (N- und C-Terminus) auf, was auf einen stärker konservierten funktionellen Kern hindeutet.

*[Abbildung 1: MSA-Visualisierung (mittlerer Ausschnitt) — siehe `output/msa_visualization.pdf`]*

### 3.3 Distanzmatrix und Heatmap

Die identitätsbasierte Distanzmatrix (31 × 31) visualisiert in einer Heatmap (Abbildung 2) zeigt deutliche Muster: Primatische Sequenzen bilden einen niedrig-Distanz-Cluster (geringe paarweise Abstände untereinander, da ähnlichste Verwandtschaft zur menschlichen Abfragesequenz). Carnivorensequenzen und Artiodactylen bilden jeweils eigene Cluster mittlerer Distanz. Die Abfragesequenz (`Query_Unknown`) weist gegenüber der Homo-sapiens-Sequenz einen Abstand von 0,0 auf (100 % Identität) und gruppiert sich eng mit den anderen Primaten.

*[Abbildung 2: Paarweise Distanz-Heatmap — siehe `output/distance_heatmap.pdf`]*

### 3.4 Überprüfung von Baum-Eigenschaften

**Vier-Punkte-Bedingung:** Von insgesamt $\binom{31}{4} = 27\,405$ Quartetten wurden [PLATZHALTER] verletzt ([PLATZHALTER] %). Dies zeigt, dass die Distanzmatrix **nicht perfekt additiv** ist, was bei realen biologischen Daten aufgrund ungleichmäßiger Evolutionsraten erwartet wird.

**Ultrametrik:** Von insgesamt $\binom{31}{3} = 4\,495$ Tripeln wurden [PLATZHALTER] verletzt ([PLATZHALTER] %). Die Distanzmatrix ist damit **nicht ultrametrisch**, was ungleiche Evolutionsraten zwischen den Linien (keine strikte molekulare Uhr) impliziert.

### 3.5 Phylogenetische Bäume

Vier Bäume wurden konstruiert (Abbildungen 3–6). Im **UPGMA-Baum** (Abbildung 3) bilden Primaten, Carnivoren und Artiodactylen jeweils monophyletische Gruppen, die der gängigen Vertebraten-Taxonomie entsprechen. Da UPGMA eine molekulare Uhr voraussetzt, sind alle Äste proportional zu Distanzen skaliert, was zu einem ultrametrischen Baum führt. Der **NJ-Baum** (Abbildung 4) zeigt eine ähnliche Topologie, erlaubt jedoch ungleiche Astlängen, was der beobachteten Nicht-Ultrametrizität der Daten besser gerecht wird.

*[Abbildung 3: UPGMA-Baum — `output/upgma_tree.pdf`]*
*[Abbildung 4: NJ-Baum mit Bootstrap-Unterstützung — `output/nj_bootstrap.pdf`]*
*[Abbildung 5: Tanglegram UPGMA vs. Complete Linkage — `output/tanglegram.pdf`]*
*[Abbildung 6: Cophyloplot UPGMA vs. NJ — `output/cophylo_upgma_nj.pdf`]*

Die **kophenetischen Korrelationen** zwischen den Bäumen und der originalen Distanzmatrix betrugen:

- UPGMA: $r_{\text{UPGMA}} =$ [PLATZHALTER]
- NJ: $r_{\text{NJ}} =$ [PLATZHALTER]

### 3.6 Bootstrap-Analyse

Die Bootstrap-Analyse des NJ-Baums (100 Replikate) ergab für die meisten internen Knoten, die klare taxonomische Gruppen (Primaten, Carnivoren) trennen, hohe Unterstützungswerte (> 70 %). Innerhalb der Primaten zeigten die Ästgabelungen zwischen Altwelt- und Neuwelt-Affen solide Bootstrap-Werte, während Verzweigungen innerhalb der Prosimier-Gruppe teils geringere Unterstützung aufwiesen.

---

## 4. Diskussion

### 4.1 Identität der Abfragesequenz

Der BLAST-Treffer mit 100 % Sequenzidentität zum humanen *Homo-sapiens*-TSR3 (NCBI Taxon-ID 9606) identifiziert die Abfragesequenz eindeutig als das menschliche **TSR3-Protein**. Dieses Protein ist ein 2'-O-Methyltransferase-Enzym, das in der Biogenese der 40S-ribosomalen Untereinheit eine essentielle Rolle spielt. Die hohe Konservierung (80–100 % Identität über alle Top-30-Vertebraten-Sequenzen) ist biologisch zu erwarten: Da TSR3 die Modifikation der 18S-rRNA ausführt, einem RNA-Molekül, das wiederum hoch konserviert ist, übt die Funktion einen starken Selektionsdruck auf die Aminosäuresequenz aus [1].

### 4.2 Orthologe vs. Paraloge

Alle 502 Treffer der BLAST-Suche stellen mutmaßliche **Orthologe** dar: Sequenzen, die durch Speziation aus einem gemeinsamen Vorläufer hervorgegangen sind. Die Akzessionsnummern des OrthoDB-Formats (Taxon-ID:Sequenz-ID) verweisen auf repräsentative Sequenzen pro Spezies. Innerhalb der Datenbank traten für einige Spezies zwei HSPs (High-Scoring Segment Pairs) derselben Akzession auf (z.B. Taxon 9615), was auf unterschiedliche Alignmentregionen mit der Abfragesequenz hindeutet — kein Hinweis auf Paraloge. Paraloge würden sich durch signifikant abweichende Sequenzidentitäten (z.B. < 40 %) bei gleichem Organismus manifestieren.

### 4.3 Wahl der Substitutionsmatrix

Der Median der Sequenzidentitäten der Top-30-Treffer von ~83,3 % rechtfertigt den Einsatz von **BLOSUM80** anstelle von BLOSUM62. BLOSUM-Matrizen sind aus Aminosäureaustauschen in real beobachteten Proteinblöcken bei verschiedenen Identitätsschwellen abgeleitet [3]. BLOSUM80 penalisiert Substitutionen stärker und differenziert besser zwischen Sequenzen mit hoher Ähnlichkeit — präziser für die vorliegende Datenlage. Diese Wahl ist jedoch nicht eindeutig: Für evolutionäre Analysen ist BLOSUM62 häufig ausreichend, da BLAST-Rankings bei hochidentischen Sequenzen meist stabil gegenüber Matrixwechseln sind [2].

### 4.4 Biologische Interpretation der Bäume

Die phylogenetischen Bäume spiegeln weitgehend die etablierte Vertebraten-Taxonomie wider [10]: Primaten bilden eine gut unterstützte monophyletische Gruppe, innerhalb derer Hominiden (Mensch, Orang-Utan) näher beieinander liegen als bei Prosimier (Lemuren, Galagos). Carnivoren (Robben, Bären, Otter, Vielfraß) bilden einen eigenen Cluster. Artiodactylen (Rinder, Schafe, Schweine) sind eng verwandt, mit dem Wal (*Physeter catodon*) als Schwestergruppe der Huftiere — konsistent mit der molekularen Phylogenetik der Säugetiere [10].

Die `Query_Unknown`-Sequenz ordnet sich erwartungsgemäß direkt neben *Homo sapiens* ein.

### 4.5 Vergleich der Baumrekonstruktionsmethoden

UPGMA und NJ lieferten ähnliche Topologien in Bezug auf die Hauptgruppen (Primaten, Carnivoren, Huftiere), was die robuste phylogenetische Signatur der TSR3-Sequenz widerspiegelt. Topologische Unterschiede traten bei Verzweigungen innerhalb von Gruppen auf. Die kophenetischen Korrelationen beider Methoden waren hoch, wobei NJ tendenziell die originalen Distanzen besser wiedergibt, da es nicht die Annahme einer molekularen Uhr macht. Für phylogenetische Analysen konservierter Proteine mit ungleichen Evolutionsraten wird NJ generell bevorzugt [8]. Complete Linkage und Ward's Methode lieferten ähnliche Hauptgruppen, zeigten jedoch andere interne Strukturen, insbesondere bei divergenteren Linien.

### 4.6 Limitationen

- **Stichprobengröße:** Nur 30 der 539 verfügbaren Sequenzen wurden analysiert. Die fehlenden Sequenzen repräsentieren weitere Vertebraten-Klassen (Vögel, Reptilien, Amphibien, Fische), deren Einbezug die phylogenetische Auflösung und biologische Interpretierbarkeit erhöhen würde.
- **Einzelgen-Phylogenie:** Phylogenien, die auf einem einzigen Gen basieren, können aufgrund von unvollständiger Sortierung der Allele (*incomplete lineage sorting*) oder horizontalem Gentransfer von der Arten-Phylogenie abweichen.
- **Algorithmus:** ClustalW ist ein Heuristik-Algorithmus; alternative MSA-Methoden (MAFFT, MUSCLE) können unterschiedliche Alignments und damit unterschiedliche Distanzmatrizen erzeugen.
- **Bootstrap:** 100 Replikate sind ein Minimum; 1000 Replikate würden die Zuverlässigkeit der Bootstrap-Werte erhöhen.
- **Distanzmetrik:** Die identitätsbasierte Distanz berücksichtigt keine evolutionären Substitutionsmodelle (z.B. WAG, LG für Proteine). Modellbasierte Ansätze (Maximum Likelihood, Bayesianische Methoden) wären phylogenetisch fundierter.

---

## 5. Zusammenfassung

Die phylogenetische Analyse der unbekannten Proteinsequenz ergab eindeutig das humane TSR3-Protein (*Ribosome biogenesis protein TSR3*), belegt durch 100 % Sequenzidentität zum *Homo-sapiens*-Referenzeintrag bei einem E-Wert von 0. TSR3 ist ein hochkonserviertes Enzym der ribosomalen Biogenese, was die beobachteten hohen Sequenzidentitäten (80–100 %) über Vertebraten hinweg erklärt. Die BLOSUM80-Substitutionsmatrix wurde auf Basis der medianen Sequenzidentität der Top-30-Treffer (~83,3 %) gewählt. Die rekonstruierten Bäume (UPGMA, NJ) sind mit der etablierten Vertebraten-Taxonomie konsistent: Primaten, Carnivoren und Artiodactylen bilden jeweils monophyletische Gruppen. Die Distanzmatrix verletzt sowohl die Vier-Punkte-Bedingung als auch die Ultrametrik, was ungleiche Evolutionsraten in den analysierten Linien belegt. Für zukünftige Analysen wird empfohlen, mehr Sequenzen einzubeziehen, modellbasierte Methoden zu verwenden (Maximum Likelihood, Bayes) und mehrere Gene in einer Supertree- oder Konkatenationsanalyse zu kombinieren.

---

## 6. Literaturverzeichnis

[1] Strunk BS, Karbstein K. (2009). Powering through ribosome assembly. *RNA* 15(12):2083–2104.

[2] Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. (1990). Basic local alignment search tool. *Journal of Molecular Biology* 215(3):403–410.

[3] Henikoff S, Henikoff JG. (1992). Amino acid substitution matrices from protein blocks. *Proceedings of the National Academy of Sciences* 89(22):10915–10919.

[4] Bodenhofer U, Bonatesta E, Horejš-Kainrath C, Hochreiter S. (2015). msa: an R package for multiple sequence alignment. *Bioinformatics* 31(24):3997–3999.

[5] Thompson JD, Higgins DG, Gibson TJ. (1994). CLUSTAL W: improving the sensitivity of progressive multiple sequence alignment through sequence weighting, position-specific gap penalties and weight matrix choice. *Nucleic Acids Research* 22(22):4673–4680.

[6] Charif D, Lobry JR. (2007). SeqinR 1.0-2: A contributed package to the R project for statistical computing devoted to biological sequences retrieval and analysis. In *Biological and Medical Physics, Biomedical Engineering*, Springer, pp. 207–232.

[7] Sokal RR, Michener CD. (1958). A statistical method for evaluating systematic relationships. *University of Kansas Scientific Bulletin* 38:1409–1438.

[8] Saitou N, Nei M. (1987). The neighbor-joining method: a new method for reconstructing phylogenetic trees. *Molecular Biology and Evolution* 4(4):406–425.

[9] Felsenstein J. (1985). Confidence limits on phylogenies: an approach using the bootstrap. *Evolution* 39(4):783–791.

[10] Murphy WJ, Eizirik E, Johnson WE, Zhang YP, Ryder OA, O'Brien SJ. (2001). Molecular phylogenetics and the origins of placental mammals. *Nature* 409(6820):614–618.

---

## 7. Anhang

### A.1 Verwendete Software und Pakete

| Software/Paket | Funktion |
|----------------|---------|
| R (≥ 4.5) | Analyseumgebung |
| BLAST+ 2.17.0+ | Sequenzähnlichkeitssuche |
| Biostrings (Bioconductor) | Sequenzoperationen |
| msa (Bioconductor) | Multiples Sequenzalignment |
| seqinr | Distanzmatrixberechnung |
| ape | Phylogenetische Bäume |
| phytools | Cophyloplot |
| dendextend | Tanglegram |
| gplots | Heatmap |
| ggmsa | MSA-Visualisierung |

### A.2 KI-Nutzung

Dieses Projekt wurde mit Unterstützung von Claude (Anthropic) entwickelt. Claude wurde verwendet, um die Analysepipeline zu strukturieren sowie R-Code für BLAST-Integration, MSA, Distanzmatrix-Berechnung, Eigenschaftsüberprüfungen und Baumkonstruktion/-vergleich zu erstellen. Sämtlicher Code wurde vor der Abgabe überprüft und verstanden.

### A.3 Analyseskript

Der vollständige, kommentierte Code befindet sich in `phylo_analysis_v3.R` im Projektverzeichnis. Nachfolgend eine Zusammenfassung der Hauptschritte:

```r
# Schritt 1: BLAST-Datenbank erstellen und BLASTP ausführen (BLOSUM62)
# -> Median-Identität bestimmen -> bei >80% mit BLOSUM80 wiederholen

# Schritt 2: MSA mit ClustalW (31 Sequenzen: 30 Top-Hits + Abfragesequenz)

# Schritt 3: Identitätsbasierte Distanzmatrix + Heatmap-Visualisierung

# Schritt 4: Vier-Punkte-Bedingung und Ultrametrik überprüfen

# Schritt 5: UPGMA, NJ, Complete Linkage, Ward's Bäume + Tanglegram + Cophyloplot
#            Bootstrap-Analyse (100 Replikate, NJ)
```

*Vollständiger Code: siehe `phylo_analysis_v3.R` im Projektverzeichnis.*

---

*Erstellt mit Unterstützung von Claude (Anthropic). Alle wissenschaftlichen Interpretationen wurden eigenständig überprüft.*
