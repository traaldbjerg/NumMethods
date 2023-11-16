void restriction(int old_m, double L, int *old_n, int **ia, int **ja, double **a, double **b, double **x, double **r) {
    /*function which restricts the fine residue to a coarser grid*/

    // get the new amount of points, variables and so forth
    int m = (old_m-1)/2 + 1;
    int q = (m-1) / 11; // nombre de fois que m-1 est multiple de 11
    int p = 4 * m - 4; // nombre de points sur le périmètre
    int nb_next_to_wall = 4 * (m-2); // nombre de points adjacents à un mur
    int nb_dir = p; // nombre de points sur une porte/fenêtre
    int invh2 = (m-1)*(m-1) / (L*L); /* h^-2 */
    //*n  = (m-1) * m; /* nombre d'inconnues */
    int *n = m * m // nombre total de points dans le carré
        - (5 * q) * (8 * q) // nombre de points dans le rectangle supérieur droit
        - nb_dir; // number of points on the walls
    printf("Value of n: %d\n", *n); // \n is necessary, otherwise the buffer is not flushed
    int nnz = 5 * (m * m - (5 * q) * (8 * q) - p) // points soumis à l'équation de la chaleur classique sans conditions particulières (on retire le rectangle supérieur droit)
        - 5 * nb_next_to_wall // we also need to remove the points directly in front of a wall (because of the dirichlet condition)
        + 4 * nb_next_to_wall // points directly in front of a wall
        ;
    int i_hole = (m-2)*q*5; // index "of the hole" in the big rectangle

    //restriction of the residue to the coarse grid
    double *r_restr = malloc(*n * sizeof(double)); // restricted residue

    int walls_hit = 0; // when a wall is hit, the next point in a is not skipped, so we keep track of it

    //fill restricted residue vector
    

}

