double rho(double x, double y){
    // épaisseur réaliste: radiateur de chez moi: largeur à peu près égale à 10 cm
    double res = 0.0, temp_rad = 50.0; // [K/m2]
    int down = 1; // représente si on prend l'emplacement du radiateur du bas ou du haut
    if (down) {
        if ((0.45 <= y) && (y <= 0.55) && ((1.5 <= x)&&(x <= 4.0))) { // intervalle du radiateur
            res = temp_rad;
        }
    }
    else {
        if ((2.45 <= y) && (y <= 2.55) && ((1.5 <= x)&&(x <= 4.0))) { // intervalle du radiateur
            res = temp_rad;
        }
    }

    //res = 0.0; // debug, permet bêtement de désactiver rho

    return res;
}
