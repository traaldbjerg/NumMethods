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
                  // basically counts the total number of j increments
    int swap = 1; // switches a condition off
    //printf("This is i_hole = %d\n", i_hole); // debug

    // fill the prolongation of the residue vector
    for (int i = 0; i < *n_coarse; i++) {
        if (i < i_hole) {
            if (row_jump % 2 == 0) {
                if (i % (m_coarse-2) == 0) {
                    //printf("In prolongation : r_coarse[%d] = %f\n", i, (*r_coarse)[i]);
                    (*r_prol)[2*i+jump] = 0.5*(*r_coarse)[i]; // nothing more because the wall wich comes right before is dirichlet, coef 0
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
                else if (i % (m_coarse-2) == m_coarse-3) {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                    (*r_prol)[2*i+jump+2] = 0.5*((*r_coarse)[i]); // nothing more because the wall which comes right after is dirichlet, coef 0
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has no coarse points at all
                    jump += 1; // last point of the row
                }
                else {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
            } else { // if we are on a fine row with no coarse points
                if (i == 0) { // first fine row, less terms in the sums -> special case
                    for (int j = 0; j < (m_coarse-2); j++) {
                        if (j % (m_coarse-2) == 0) { // next to the first wall
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1]);
                        } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1]);
                            (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j-1+1]);
                            row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                            jump += 2 * (m_coarse-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                            i--; // we have not iterated over the coarse vector, so we need to go back one step, otherwise i gets too big
                        } else {
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1]);
                        }
                    } 
                } else { 
                    for (int j = 0; j < (m_coarse-2); j++) {
                        if (j % (m_coarse-2) == 0) { // next to the first wall
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3]);
                        } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3] + (*r_coarse)[i+j-1] + (*r_coarse)[i+j-1-m_coarse+2]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3]);
                            (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3]);
                            row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                            jump += 2 * (m_coarse-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                            i--; // we have not iterated over the coarse vector, so we need to go back one step, otherwise i gets too big
                        } else {
                            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i-m_coarse+3] + (*r_coarse)[i+j-1] + (*r_coarse)[i+j-1-m_coarse+2]);
                            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i-m_coarse+3]);
                        }
                    }
                }
            }
        // still need to take care of the last fine row before the smaller rectangle
        // after the cutout, the points above are on the dirichlet wall -> residue = 0
        // less terms in the sum
        } else if (swap && (i == i_hole)) {
            //printf("Last fine row before the cutout\n");
            //printf("This is i = %d\n", i);
            //printf("This is i_hole = %d\n", i_hole);
            //printf("\n");
            // PROBLEMS HERE APPARENTLY
            for (int j = 0; j < (m_coarse-2); j++) {
                if (j % (m_coarse-2) == 0) { // next to the first wall
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-m_coarse+3]);
                } else if (j % (m_coarse-2) == m_coarse-3) { // next to the second wall
                    //printf("In the second wall\n");
                    //printf("This is r_coarse[%d] = %f\n", i+j-1-m_coarse+2, (*r_coarse)[i+j-1-m_coarse+2]);
                    //printf("This is r_coarse[%d] = %f\n", i+j-1-m_coarse+3, (*r_coarse)[i+j-1-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1-m_coarse+3] + (*r_coarse)[i+j-1-m_coarse+2]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j-1-m_coarse+3]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j, (*r_prol)[2*i+jump+2*j]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+1, (*r_prol)[2*i+jump+2*j+1]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+2, (*r_prol)[2*i+jump+2*j+2]);
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time
                    jump += 2 * (m_coarse-2) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                    i--;// we have not iterated over the coarse vector, so we need to go back one step, otherwise i gets too big
                    swap = 0; // switch the condition off 
                } else if (j == (3 * q_coarse - 1)) { // next to the obtuse corner
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i-m_coarse+3] + (*r_coarse)[i+j-1] + (*r_coarse)[i+j-1-m_coarse+2]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i-m_coarse+3]);
                } else if (j > (3 * q_coarse - 1)) { // next to the upper right wall
                                                   // this condition was the source of many bugs because it needs to be after the end of row one
                                                   // otherwise the end of wall condition is never reached
                    //printf("In the upper right corner\n");
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1-m_coarse+2] + (*r_coarse)[i+j-1-m_coarse+3]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1-m_coarse+3]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j, (*r_prol)[2*i+jump+2*j]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+1, (*r_prol)[2*i+jump+2*j+1]);
                } else {
                    (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i-m_coarse+3] + (*r_coarse)[i+j-1] + (*r_coarse)[i+j-1-m_coarse+2]);
                    (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i-m_coarse+3]);
                }
            }
        // fill the vector with the second rectangle
        } else { // if we are in the upper small rectangle
            if (row_jump % 2 == 0) {
                if ((i-i_hole) % (3 * q_coarse-1) == 0) {
                    (*r_prol)[2*i+jump] = 0.5*(*r_coarse)[i]; // nothing more because the wall wich comes right before is dirichlet, coef 0
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                    //printf("Hey this works for coarse rows\n");
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump, (*r_prol)[2*i+jump]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+1, (*r_prol)[2*i+jump+1]);
                    //printf("This is i = %d\n", i);
                }
                else if ((i-i_hole) % (3 * q_coarse-1) == 3 * q_coarse-2) {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                    (*r_prol)[2*i+jump+2] = 0.5*((*r_coarse)[i]); // nothing more because the wall which comes right after is dirichlet, coef 0
                    //printf("Hey this works for coarse rows, END OF ROW\n");
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump, (*r_prol)[2*i+jump]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+1, (*r_prol)[2*i+jump+1]);
                    //printf("This is r_prol[%d] = %f\n", 2*i+jump+2, (*r_prol)[2*i+jump+2]);
                    //printf("This is i = %d\n", i);
                    row_jump++; // we have hit a wall, so we need to jump to the next row which has no coarse points at all
                    jump += 1; // last point of the row
                } else {
                    (*r_prol)[2*i+jump] = 0.5*((*r_coarse)[i] + (*r_coarse)[i-1]);
                    (*r_prol)[2*i+jump+1] = (*r_coarse)[i];
                }
            } else { // if we are on a fine row with no coarse points
                // need to take care of the last row of the small rectangle, special case with less terms in the sums OUTSIDE OF THE FOR LOOP
                for (int j = 0; j < 3 * q_coarse-1; j++) {
                    if (j % (3 * q_coarse-1) == 0) { // next to the first wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-3 * q_coarse+2]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-3 * q_coarse+2]);
                        //printf("Hey this works for fine rows\n");
                        //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j, (*r_prol)[2*i+jump+2*j]);
                        //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+1, (*r_prol)[2*i+jump+2*j+1]);
                        //printf("This is i = %d\n", i);
                    } else if (j % (3 * q_coarse-1) == 3 * q_coarse-2) { // next to the second wall
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-3 * q_coarse+2] + (*r_coarse)[i+j-1] + (*r_coarse)[i+j-1-3 * q_coarse+1]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-3 * q_coarse+2]);
                        (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i+j-1-3 * q_coarse+2]);
                        //printf("Hey this works for fine rows, END OF ROW\n");
                        //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j, (*r_prol)[2*i+jump+2*j]);
                        //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+1, (*r_prol)[2*i+jump+2*j+1]);
                        //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+2, (*r_prol)[2*i+jump+2*j+2]);
                        //printf("This is i = %d\n", i);
                        //printf("This is r_coarse[%d] = %f\n", i+j-1-3 * q_coarse+2, (*r_coarse)[i+j-1-3 * q_coarse+2]);
                        row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time 
                        jump += 2 * (3 * q_coarse-1) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
                        i--; // we have not iterated over the coarse vector, so we need to go back one step, otherwise i gets too big
                    } else {
                        (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i-3 * q_coarse+2] + (*r_coarse)[i+j-1] + (*r_coarse)[i+j-1-3 * q_coarse+1]);
                        (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1+1] + (*r_coarse)[i-3 * q_coarse+2]);
                    }
                }
            }
        }
    }

    // do the last row of the fine grid outside of the for loop
    int i = *n_coarse;
    for (int j = 0; j < 3 * q_coarse-1; j++) {
        if (j % (3 * q_coarse-1) == 0) { // next to the first wall
            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1-3 * q_coarse+2]);
            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1-3 * q_coarse+2]);
            //printf("Hey this works for fine rows\n");
            //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j, (*r_prol)[2*i+jump+2*j]);
            //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+1, (*r_prol)[2*i+jump+2*j+1]);
        } else if (j % (3 * q_coarse-1) == 3 * q_coarse-2) { // next to the second wall
            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i+j-1-3 * q_coarse+2] + (*r_coarse)[i+j-1-3 * q_coarse+1]);
            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i+j-1-3 * q_coarse+2]);
            (*r_prol)[2*i+jump+2*j+2] = 0.25 * ((*r_coarse)[i+j-1-3 * q_coarse+2]);
            //printf("Hey this works for fine rows, END OF ROW\n");
            //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j, (*r_prol)[2*i+jump+2*j]);
            //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+1, (*r_prol)[2*i+jump+2*j+1]);
            //printf("This is r_prol[%d] = %f\n", 2*i+jump+2*j+2, (*r_prol)[2*i+jump+2*j+2]);
            row_jump++; // we have hit a wall, so we need to jump to the next row which has coarse points this time 
            jump += 2 * (3 * q_coarse-1) + 1; // we have jumped a row, so we need to add the number of fine points in a row to the jump
            i--; // we have not iterated over the coarse vector, so we need to go back one step, otherwise i gets too big
        } else {
            (*r_prol)[2*i+jump+2*j] = 0.25 * ((*r_coarse)[i-3 * q_coarse+2] + (*r_coarse)[i+j-1-3 * q_coarse+1]);
            (*r_prol)[2*i+jump+2*j+1] = 0.5 * ((*r_coarse)[i-3 * q_coarse+2]);
        }
    }

}
