void compute_heat_flux(int m, int q, double **x, double source, double *sum_grad_x, double *sum_grad_y, double *rad_power, double (*source_func)(double, double, double)) {
    // intégrer le gradient de température sur la porte et la fenêtre

    double k = 0.02512; // conductivité de l'air à 10°C et 1 atm [W/(m K)]
    double L = 5.5;
    double h = L/(m-1);

    // flux par la fenêtre, iy == 1; q*3 <= ix; ix <= q*8
    // quel indice de x correspond au point directement au-dessus du coin gauche de la fenêtre?
    // m - q*5 - 1 + q*3 (on ne rajoute pas un +1 supplémentaire car les indices commencent à 0) 
    // gradient selon x pas à considérer pour la fenêtre car 1n et 1x orthogonaux
    // donc pas nécessaire de considérer grad_x aux extrémités de la fenêtre

    *sum_grad_x = 0.0;
    *sum_grad_y = 0.0;

    for (int ind = m - q*5 - 1 + q*3; ind < m + q*3; ind++) {
        // différence avant
        *sum_grad_y += -k * ( (*x)[ind]/2.0 + (*x)[ind + 1]/2.0); //pow((x[ind] - 0.0) / h, 2.0); // on devrait rajouter un signe moins car normale de la surface mais perte de temps de calcul car de toute façon mis au carré
                                         // on ne divise pas par h car on remultiplie par h (élément dl d'intégration)
    }

    // flux par la porte, ix == 1; q*6 <= iy; iy <= q*8
    // quel indice de x correspond au point directement à droite du bord inférieur de la porte?
    // 6 * q * m - 5 * q - 1
    // indice final?
    // 6 * q * m - 5 * q - 1 + (m-1) + (2 * q - 1) * 3 * q

    for (int ind = 6 * q * m - 5 * q - 1; ind < 6 * q * m - 5 * q - 1 + (m-1) + (2 * q - 1) * 3 * q; ind += 3 * q) {
        // différence avant
        
        if (ind == 6 * q * m - 5 * q - 1) {
            *sum_grad_x += -k * (((*x)[ind] - 20.0)/2.0 + ((*x)[ind + m - 1] - 20.0)/2.0);
            ind += - 3 * q + m - 1; // la première fois il faut longer tout le mur supérieur droit avant d'atteindre le prochain point en face de la porte
        } else {
            *sum_grad_x += -k * (((*x)[ind] - 20.0)/2.0 + ((*x)[ind + 3 * q] - 20.0)/2.0); //pow((x[ind] - 20.0) / h, 2.0); // on devrait rajouter un signe moins mais perte de temps de calcul car de toute façon mis au carré
                                          // on ne divise pas par h car on remultiplie par h (élément dl d'intégration)
        }
    }

    // estimer le flux sortant du radiateur
    // rho de valeur constante sur le radiateur -> il suffit de faire rho * S * k pour obtenir les bonnes valeurs

    // intégrer sur le radiateur permet de prendre en compte son épaisseur correctement sans se soucier de si des points tombent pile sur le périmètre du radiateur ou non

    // il faut une méthode des trapèzes 2D -> on interpole sur l'interpolation?
    // rho constante donc uniquement intéressant pour points à l'extrémité

    *rad_power = 0.0;

    for (int iy = 1; iy < q*6; iy++) { // inégalités strictes -> on ne prend pas en compte les points sur le demi-périmètre supérieur droit -> intégrale respectée
        for (int ix = q * 3; ix < q*8; ix++) {
            //if (ix == q*3 || iy)
            *rad_power += (*source_func)(ix*h, iy*h, source) * k * h * h; // h2 pour l'élément de surface dS
        }
    }

}
