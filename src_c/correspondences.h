/* Header for the old correspondences code. */

#ifndef CORRESPONDENCES_H
#define CORRESPONDENCES_H

#include "ptv.h"
#include <optv/parameters.h>
#include <optv/calibration.h>

int correspondences_4 (target pix[][nmax], coord_2d geo[][nmax], int num[],
    volume_par *vpar, control_par *cpar, Calibration cals[], n_tupel *con, int match_counts[]);

#endif
