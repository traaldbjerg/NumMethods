# CalculNum2223

Ce code, écrit dans le cadre du projet de Calcul Numérique pour le cours MATH-H301 permet de générer et résoudre un problème décrit par l'équation stationnaire de la chaleur, avec pour domaine et conditions au bord celles de la pièce numéro 23.


Ce projet a été réalisé sur un WSL Ubuntu, donc les limitations de RAM n'étaient pas un problème: il est très possible qu'il faille augmenter la quantité de RAM s'il est exécuté sur la machine virtuelle padi-VM sans changer les paramètres de base (comme le nombre d'itérations ou le m).

## Question 1 : Génération du problème

Le fichier prob.c génère le problème qui doit être résolu, qui à cause de la discrétisation peut se résumer par un système linéaire: Ax = b.
Le programme génère la matrice A et le vecteur b sous un format CSR qui permet d'optimiser l'espace mémoire occupé par les variables qui les forment en ignorant toutes les valeurs nulles. Le format CSR s'est avéré compliqué à utiliser, surtout pour le vecteur ja, étant donné que la longueur de la pièce (où plutôt le nombre d'inconnues sur une ligne horizontale) varie, que ce soit à cause de la largeur de la pièce en elle-même ou bien à cause des conditions de Dirichlet qui font perdre des inconnues.

## Question 2 : Calcul du résidu

Cette sous-question est résolue dans le fichier residue.c (qui définit une fonction appelée en bas de main.c). Il suffit simplement de créer explicitement le système A et de soustraire le carré de chaque ligne au carré de la composante correspondante de b, puis de prendre la racine de la somme de ces composantes et de la diviser par la norme de b.

## Question 3 : Affichage

Dans le fichier main.c, une section a été ajoutée afin de créer le fichier out.dat sur base du vecteur x. En effet, ce vecteur est uniquement composé des points dont la température est à résoudre sur base de l'équation de la chaleur. Or, les points de la fenêtre et de la porte ne sont pas des inconnues puisque leur valeur est imposée par les conditions de Dirichlet sur la porte et la fenêtre. Il a donc fallu faire attention à rajouter ces points en particulier avant de fournir out.dat à gnuplot, qui nous affiche sur base de ces données un heatmap.

## Question 4 : Fonction source

Dans un nouveau fichier rho.c, la fonction rho a été déclarée, fonction de x et de y (et d'autres paramètres comme la température du radiateur) et est utilisée dans prob.c (à l'aide d'un pointeur de fonction) pour ajouter aux composantes du vecteur b un terme source quand le point en question se trouve sur un emplacement de radiateur.

L'itération peut être désactivée en assignant iter_max à 0.

        recherche de l'optimum de température du radiateur, recherche linéaire, source_value = 0.000 jusque 999.0 par pas de 1

    radiateur en bas:
        meilleur std_dev/avg: 0.160079 pour source_value = 534.000
        meilleur std_dev: 2.939890 pour source_value = 434.000

    radiateur en haut:
        meilleur std_dev/avg: 0.367950 pour source_value = 174.000
        meilleur std_dev: 5.050443 pour source_value = 67.000

    En se basant sur le critère d'écart-type, il semblerait que le radiateur en bas de la pièce mène à la configuration de température la plus uniforme.

    Visuellement, std_dev/avg semble être un meilleur critère d'uniformité.
    C'est logique étant donné que ce critère ne pénalise pas particulièrement les hautes températures.

Le flux de chaleur sortant par la porte et la fenêtre a également été estimé : le flux a été estimé en intégrant la densité de flux de chaleur (à savoir -k * dT/d[x ou y]) sur la longueur de la porte/fenêtre). La puissance du radiateur, elle, a été obtenue en intégrant k * rho sur la surface effective (c'est-à-dire la surface "vue" par la discrétisation) du radiateur.
Pour un radiateur éteint, le flux de chaleur traversant la porte est d'environ 0.261594 W, tandis que par la fenêtre on a -0.269866 W (calculés à m = 1101). Ces 2 valeurs ne sont pas tout à fait égales mais ceci s'explique simplement par la quantité d'approximations réalisées de par la nature de la discrétisation des équations continues.

Pour un radiateur allumé, pour le rho tel que la distribution de température de la pièce est la plus uniforme (rho = 534.000 K/m2), on voit que la somme des flux de chaleur à travers la porte/fenêtre est très proche de la puissance du radiateur, de l'ordre de 3.3 W, avec un écart de l'ordre de 1% pour des valeurs de m > 150. Ceci est logique : ce code traite de l'équation de la chaleur sous sa forme stationnaire, ce qui implique qu'il ne peut y avoir de variation de température au cours du temps. De ce fait, il ne peut y avoir accumulation de chaleur, et donc toute la puissance dissipée par le radiateur doit être compensée par les pertes au niveau de la porte/fenêtre.

En chipotant quelque peu avec les nombres, on s'aperçoit rapidement que les valeurs de rho telles que la pièce est uniforme sont toujours en dessous de 1000, ce qui justifie le choix de la valeur de m à 166 permettant à la fois une précision acceptable tout en mettant moins de 2 minutes pour générer et résoudre 1000 fois le problème.

## Question 5 : Solveur open-source

La librairie choisie est PETSc, en partie car elle s'avéra avoir des instructions d'installation relativement compréhensibles (bien qu'il ait quand même fallu beaucoup de chipotage).

Le makefile qui compile main.c et toutes ses dépendances a été adapté pour avoir accès aux fonctionnalités de PETSc.

Pour que PETSc s'exécute sur votre machine, placez-vous dans le répertoire petsc et exécutez dans un terminal

    ./configure --with-fc=0 --with-mpi=0

puis

    make PETSC_DIR=/home/alexs/code/cn/project/CalculNum2223/CodePub/petsc PETSC_ARCH=arch-linux-c-debug all

où le chemin pour la variable PETSC_DIR devra probablement être modifié selon l'endroit où vous placez ce projet.
Ensuite, pour que main.c puisse avoir accès à la librairie, il faut indiquer à Linux où trouver la librairie. Exécutez

    export LD_LIBRARY_PATH+=/pathtothisproject/petsc/arch-linux-c-debug/lib/

dans un terminal (il faudra le refaire si vous éteignez votre machine). Sur ma machine, cela donne :

    export LD_LIBRARY_PATH+=/home/alexs/code/cn/project/CalculNum2223/CodePub/petsc/arch-linux-c-debug/lib/

Une fois ceci fait, lancer make dans le répertoire où se trouve main.c devrait fonctionner.

Sans préconditionneur, la comparaison montre que PETSc et UMFPACK aboutissent à des solutions relativement proches, par exemple pour le tout premier point UMFPACK trouve une température de 17.790890, tandis que PETSc trouve une température de 17.783690 (m = 166, rho = 534.00 K/m2). L'affichage des 2 solutions est visuellement identique. Cependant, PETSc (sans préconditionnement custom) est environ 5 fois plus lent. De plus, le résidu atteint par PETSc sans préconditionnement et sans tolérance prédéfinie est mauvais (de l'ordre de 3.93e-6 au lieu de 1.04e-14 pour UMFPACK).

En rajoutant un préconditionneur PCILU (PCJACOBI est sensiblement plus lent) et en mettant une tolérance de 1e-15, le résidu final est comparable à UMFPACK, cependant le temps de calcul devient beaucoup plus long! En effet, pour m = 331 et le radiateur éteint on a :

    Temps de solution, UMFPACK (CPU):   0.3 sec
    Temps de solution, UMFPACK (horloge):   0.3 sec

    Résidu de la solution: 5.7279914362e-15

    Temps de calcul du résidu UMFPACK (CPU):   0.0 sec
    Temps de calcul du résidu UMFPACK (horloge):   0.0 sec

    Temps de solution, PETSc (CPU):  20.9 sec
    Temps de solution, PETSc (horloge):  20.9 sec

    Résidu de la solution PETSc: 6.0074350247e-15

Il est possible que ceci soit dû au fait que PETSc a été compilé sans MPI, qui semble être une composante importante de la philosophie de design de PETSc (à en croire le site).