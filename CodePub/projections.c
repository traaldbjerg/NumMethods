void restriction(int m, int q, int *n, int **ia, int **ja, double **a, double **b, double **x, double **r, double **r_restr) {
    /*function which restricts the fine residue to a coarser grid*/

    // get the new amount of points, variables and so forth
    
    int i_hole = (m-2)*q*6; // index of start of the row "of the hole" in the big rectangle

    //restriction of the residue to the coarse grid
    int jump = 0; // when a row is jumped (wall hit), we need to keep track of it + also keep track of the consecutive points in front of a wall

    printf("This is i_hole = %d\n", i_hole);
    printf("This is m = %d\n", m);
    printf("This is q = %d\n", q);
    printf("This is n = %d\n", *n);
    //fill restricted residue vector
    for (int i = 0; i < *n; i++) {
        (*r_restr)[i] = (*r)[2*i + jump];
        //printf("This is r[%d] = %f\n", 2*i+jump, (*r)[2*i+jump]);
        //printf("This is 2*i + jump = %d\n", 2*i+jump);
        if (((2*i+jump) % (m-2) == 0) && (i != 0) && ((2*i+jump) <= i_hole)) { // if we are in the bigger rectangle of the room
            jump += m-3; // m-2 - 1 because start of row is at index +1 of end of row and not +2
        }
        else if (((2*i+jump-i_hole) % (3*q-1) == 0) && ((2*i+jump-i_hole) != 0) && ((2*i+jump) > i_hole)){
            jump += 3*q-2; // -2 because dirichlet condition and -1 because start of row is at index +1 of end of row and not +2
        }
    }

}

void prolongation(int m, int q, int *n, double **r_coarse, double **r_prol) {
    // interpolates the coarse grid residue to the fine grid

    int i_hole = (m-2)*q*6; // index of start of the row "of the hole" in the big rectangle

    // fill the prolongation of the residue vector
    for (int i = 0; i < *n; i++) {
        (*r_prol)[2*i] = (*r_coarse)[i];
        if (((2*i) % (m-2) == 0) && (i != 0) && ((2*i) <= i_hole)) { // if we are in the bigger rectangle of the room
            (*r_prol)[2*i] += (*r_coarse)[i-1];
        }
        else if (((2*i-i_hole) % (3*q-1) == 0) && ((2*i-i_hole) != 0) && ((2*i) > i_hole)){
            (*r_prol)[2*i] += (*r_coarse)[i-1];
        }
        else {
            (*r_prol)[2*i] += (*r_coarse)[i-1] + (*r_coarse)[i];
        }
        (*r_prol)[2*i] /= 2;
    }

}
