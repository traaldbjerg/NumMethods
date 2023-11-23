void restriction(int m, int q, int *n, int **ia, int **ja, double **a, double **b, double **x, double **r, double **r_restr) {
    /*function which restricts the fine residue to a coarser grid*/

    // get the new amount of points, variables and so forth
    int m_fine = (m-1)*2 + 1; // number of points in a row of the fine grid
    int q_fine = (m_fine-1)/11; // number of points in a row of the fine grid in the small rectangle
    int i_hole = (m_fine-2)*(q_fine*6-1); // index of start of the row "of the hole" in the big rectangle
    //printf("This is i_hole = %d\n", i_hole); // debug

    //restriction of the residue to the coarse grid
    int jump = m_fine-1; // when a row is jumped (wall hit), we need to keep track of it + also keep track of the consecutive points in front of a wall
                    // starts at m-1 because the first row is jumped, as it only contains fine points and no coarse points
                    // m-1 is the coordinate of the first point to take for the coarse vector

    //fill restricted residue vector
    for (int i = 0; i < *n; i++) {
        (*r_restr)[i] = (*r)[2*i + jump];
        //printf("This is r[%d] = %f\n", 2*i+jump, (*r)[2*i+jump]); // debug
        if (((2*i+jump+2) % (m_fine-2) == 0) && (i != 0) && ((2*i+jump) <= i_hole)) { // if we are in the bigger rectangle of the room
            jump += m_fine-1; // m because m-2 points per row, plus the last of the previous coarse line and the first of the next coarse line
        } 
        else if (((2*i+jump-i_hole+2) % (3*q_fine-1) == 0) && ((2*i+jump-i_hole) != 0) && ((2*i+jump) > i_hole)){
            jump += 3*q_fine; // 3*q-1 per row in the small rectangle, plus the last of the previous coarse line and the first of the next coarse line
        }
    }
}

void prolongation(int m_coarse, int q_coarse, int *n_coarse, double **r_coarse, double **r_prol) {
    // interpolates the coarse grid residue to the fine grid

    int i_hole = (m_coarse-2)*(q_coarse*6-1); // index of start of the row "of the hole" in the big rectangle
    int row_jump = 1; // when a row is jumped (wall hit), we need to keep track of it + also keep track of the consecutive points in front of a wall
                      // starts at 1 because the first row is a fine row with no coarse points
    int jump = 0; // when a fine row is done, we need to keep track of it so that we know how many fine points have been dealt with already in the fine vector

    // fill the prolongation of the residue vector
    for (int i = 0; i < *n_coarse; i++) {
        if (i < i_hole-1) {
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
                if (i == 0) { // first fine row, less terms in the sums -> special case
                    for (int j = 0; j < (m_coarse-2); j++) {
                        if (j % (m_coarse-2) == 0) { // next to the first wall
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1]);
                        } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1]);
                            (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j+1]);
                            row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                            jump += 2 * (m_coarse-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                        } else {
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1]);
                        }
                    } 
                } else { 
                    for (int j = 0; j < (m_coarse-2); j++) {
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
            }
        // still need to take care of the last fine row before the smaller rectangle
        // after the cutout, the points above are on the dirichlet wall -> residue = 0
        // less terms in the sum
        } else if (i == i_hole-1) {
            for (int j = 0; j < (m_coarse-2); j++) {
                if (j % (m_coarse-2) == 0) { // next to the first wall
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-m_coarse+3]);
                } else if (j > (3*q_coarse - 1)) { // next to the upper right wall
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-m_coarse+2] + (*r_coarse)[i+j-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-m_coarse+3]);                              
                } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+m_coarse-3] + (*r_coarse)[i+j-m_coarse+2]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+m_coarse-3]);
                    (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j+m_coarse-3]);
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                }  else {
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-m_coarse+3] + (*r_coarse)[i+j] + (*r_coarse)[i+j-m_coarse+2]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-m_coarse+3]);
                    }
            }
        // fill the vector with the second rectangle
        } else { // if we are in the upper small rectangle
            if (row_jump % 2 == 0) {
                if ((i-i_hole) % (3*q_coarse-1) == 0) {
                    (*r_prol)[2*i+jump] = 0.5*(*r_coarse)[i]; // nothing more because the wall wich comes right before is dirichlet, coef 0
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
                else if ((i-i_hole) % (3*q_coarse-1) == 3*q_coarse-2) {
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
                for (int j = 0; j < 3*q_coarse-1; j++) {
                    if (j % (3*q_coarse-1) == 0) { // next to the first wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-3*q_coarse+2]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j-3*q_coarse+2]);
                    } else if (j % (3*q_coarse-1) == 3*q_coarse-2) { // next to the second wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+3*q_coarse-2] + (*r_coarse)[i+j] + (*r_coarse)[i+j-3*q_coarse+1]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+3*q_coarse-2]);
                        (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i+j+3*q_coarse-2]);
                        row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                        if ((i == *n_coarse-1) && (j == 3*q_coarse-2)) {
                            for (int k = 0; k < 2 * i + jump + 2 * j + 3; k++) {
                                printf("This is r_prol[%d] = %f\n", k, (*r_prol)[k]);
                            }
                        }
                        jump += 2 * (3*q_coarse-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                    } else {
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-3*q_coarse+2] + (*r_coarse)[i+j] + (*r_coarse)[i+j-3*q_coarse+1]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j+1] + (*r_coarse)[i-3*q_coarse+2]);
                    }
                }
            }
        }
    }
}
