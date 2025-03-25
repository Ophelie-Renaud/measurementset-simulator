#ifndef MAJOR_LOOP_ACTOR_H
#define MAJOR_LOOP_ACTOR_H

#ifdef __cplusplus
extern "C" {
#endif

	#include "preesm.h"

	void iterator(int ITER, int START, OUT int *cycle_out);

	void sink_gains(int NUM_RECEIVERS, IN PRECISION2 *input);

#ifdef __cplusplus
}
#endif


#endif
