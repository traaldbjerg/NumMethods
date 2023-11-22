void restriction(int m, int q, int *n, int **ia, int **ja, double **a, double **b, double **x, double **r, double **r_restr) {
    /*function which restricts the fine residue to a coarser grid*/

    // get the new amount of points, variables and so forth
    
    int i_hole = (m-2)*q*6; // index of start of the row "of the hole" in the big rectangle

    //restriction of the residue to the coarse grid
    int jump = 0; // when a row is jumped (wall hit), we need to keep track of it + also keep track of the consecutive points in front of a wall

    //printf("This is i_hole = %d\n", i_hole);
    //printf("This is m = %d\n", m);
    //printf("This is q = %d\n", q);
    //printf("This is n = %d\n", *n);
    //fill restricted residue vector
    for (int i = 0; i < *n; i++) {
        (*r_restr)[i] = (*r)[2*i + jump];
        printf("This is r[%d] = %f\n", 2*i+jump, (*r)[2*i+jump]);
        //printf("This is 2*i + jump = %d\n", 2*i+jump);
        if (((2*i+jump) % (m-2) == 0) && (i != 0) && ((2*i+jump) <= i_hole)) { // if we are in the bigger rectangle of the room
            jump += m-3; // m-2 - 1 because start of row is at index +1 of end of row and not +2
        }
        else if (((2*i+jump-i_hole) % (3*q-1) == 0) && ((2*i+jump-i_hole) != 0) && ((2*i+jump) > i_hole)){
            jump += 3*q-2; // -2 because dirichlet condition and -1 because start of row is at index +1 of end of row and not +2
        }
    }

}

void prolongation(int m_coarse, int q, int *n_coarse, double **r_coarse, double **r_prol) {
    // interpolates the coarse grid residue to the fine grid

    int i_hole = (m_coarse-2)*q*6; // index of start of the row "of the hole" in the big rectangle
    int row_jump = 0; // when a row is jumped (wall hit), we need to keep track of it + also keep track of the consecutive points in front of a wall
    int jump = 0; // when a fine row is done, we need to keep track of it so that we know how many fine points have been dealt with already in the fine vector

    // fill the prolongation of the residue vector
    for (int i = 0; i < *n_coarse; i++) {
        if (i < i_hole) {
            if (row_jump % 2 == 0) {
                if (i % (m_coarse-2) == 0) {
                    (*r_prol)[2*i+jump] = 0.5*(*r_coarse)[i]; // nothing more because the wall wich comes right before is dirichlet, coef 0
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
                else if (i % (m_coarse-2) == m_coarse-3) {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                    (*r_prol)[2*i+jump+2] = 0.5*((*r_coarse)[i]); // nothing more because the wall which comes right after is dirichlet, coef 0
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has no coarse points at all

                }
                else {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
            } else { // if we are on a fine row with no coarse points
                for (int j = 0; j < m_coarse-2; j++) {
                    if (j % (m_coarse-2) == 0) { // next to the first wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-m_coarse+3]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-m_coarse+3]);
                    } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+m_coarse-3] + (*r_coarse)[i+j] + (*r_coarse)[i+j-m_coarse+2]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+m_coarse-3]);
                        (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+m_coarse-3]);
                        row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                        jump += 2 * (m_coarse-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                    } else {
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-m_coarse+3] + (*r_coarse)[i+j] + (*r_coarse)[i+j-m_coarse+2]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-m_coarse+3]);
                    }
                }
            }
        // still need to take care of the last fine row before the smaller rectangle
        // after the coutout, the points above are on the dirichlet wall -> residue = 0
        // less terms in the sum
        } else if (i == i_hole) {
            for (int j = 0; j < m_coarse-2; j++) {
                if (j % (m_coarse-2) == 0) { // next to the first wall
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-m_coarse+3]);
                } else if (j > (3*q - 1)) {
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-m_coarse+2] + (*r_coarse)[i+j-m_coarse+3]); 
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-m_coarse+3]);                              
                } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+m_coarse-3] + (*r_coarse)[i+j-m_coarse+2]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+m_coarse-3]);
                    (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j+m_coarse-3]);
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                }
            }
        // fill the vector with the second rectangle
        } else { // if we are in the upper small rectangle
            if (row_jump % 2 == 0) {
                if ((i-i_hole) % (3*q-1) == 0) {
                    (*r_prol)[2*i+jump] = 0.5*(*r_coarse)[i]; // nothing more because the wall wich comes right before is dirichlet, coef 0
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
                else if ((i-i_hole) % (3*q-1) == 3*q-2) {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                    (*r_prol)[2*i+jump+2] = 0.5*((*r_coarse)[i]); // nothing more because the wall which comes right after is dirichlet, coef 0
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has no coarse points at all

                }
                else {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
            } else { // if we are on a fine row with no coarse points
                for (int j = 0; j < 3*q-1; j++) {
                    if (j % (3*q-1) == 0) { // next to the first wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-3*q+2]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-3*q+2]);
                    } else if (j % (3*q-1) == 3*q-2) { // next to the second wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+3*q-2] + (*r_coarse)[i+j] + (*r_coarse)[i+j-3*q+1]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+3*q-2]);
                        (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+3*q-2]);
                        row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                        jump += 2 * (3*q-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                    } else {
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-3*q+2] + (*r_coarse)[i+j] + (*r_coarse)[i+j-3*q+1]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-3*q+2]);
                    }
                }
            }
        }
    }
}
