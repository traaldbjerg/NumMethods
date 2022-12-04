# CalculNum2223

Ce code, écrit dans le cadre du projet de Calcul Numérique pour le cours MATH-H301 permet de générer et résoudre un problème décrit par l'équation stationnaire de la chaleur, avec pour conditions au bord celles de la pièce numéro 23.

## Question 1 : génération du problème

Le fichier prob.c génère le problème qui doit être résolu, qui à cause de la discrétisation peut se résumer par un système linéaire: Ax = b.
Le programme génère la matrice A et le vecteur b sous un format CSR qui permet d'optimiser l'espace mémoire occupé par les variables qui les forment en ignorant toutes les valeurs nulles. Le format CSR s'est avéré compliqué à utiliser, surtout pour le vecteur ja, étant donné que la longueur de la pièce (où plutôt le nombre d'inconnues sur une ligne horizontale) varie, que ce soit à cause de la largeur de la pièce en elle-même ou bien à cause des conditions de Dirichlet qui font perdre des inconnues.

## Question 2 : Calcul du résidu

Cette sous-question est résolue dans le fichier residue.c (qui définit une fonction appelée en bas de main.c). Il suffit simplement de créer explicitement le système A et de soustraire le carré de chaque ligne au carré de la composante correspondante de b, puis de prendre la racine de la somme de ces composantes et de la diviser par la norme de b.

## Question 3 : Affichage

Dans le fichier main.c, une section a été ajoutée afin de créer le fichier out.dat sur base du vecteur x. En effet, ce vecteur est uniquement composé des points dont la température est à résoudre sur base de l'équation de la chaleur. Or, les points de la fenêtre et de la porte ne sont pas des inconnues puisque leur valeur est imposée par les conditions de Dirichlet sur la porte et la fenêtre. Il a donc fallu faire attention à rajouter ces points en particulier avant de fournir out.dat à gnuplot, qui nous affiche sur base de ces données un heatmap.

## Question 4 : Fonction source

Dans un nouveau fichier rho.c, la fonction rho a été déclarée, fonction de x et de y (et d'autres paramètres comme la température du radiateur) et est utilisée dans prob.c (à l'aide d'un pointeur de fonction) pour ajouter aux composantes du vecteur b un terme source quand le point en question se trouve sur un emplacement de radiateur

        recherche de l'optimum de température du radiateur, recherche linéaire, temp_rad = 0.000 jusque 999.0 par pas de 1

    radiateur en bas:
        meilleur std_dev/avg: 0.160079 pour temp_rad = 534.000
        meilleur std_dev: 2.939890 pour temp_rad = 434.000

    radiateur en haut:
        meilleur std_dev/avg: 0.367950 pour temp_rad = 174.000
        meilleur std_dev: 5.050443 pour temp_rad = 67.000

    En se basant sur le critère d'écart-type, il semblerait que le radiateur en haut de la pièce mène à la configuration de température la plus uniforme.

    Visuellement, std_dev/avg semble être un meilleur critère d'uniformité.
    C'est logique étant donné que ce critère pénalise moins les hautes températures.

Le flux de chaleur sortant par la porte et la fenêtre a également été estimé.

