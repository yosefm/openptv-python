/* Definitions for tracking routines. */

#ifndef TRACK_H
#define TRACK_H

#include "tracking_run.h"

tracking_run* trackcorr_c_init();
void trackcorr_c_loop (tracking_run *run_info, int step, int display);
void trackcorr_c_finish(tracking_run *run_info, int step);
void trackback_c();

#endif
