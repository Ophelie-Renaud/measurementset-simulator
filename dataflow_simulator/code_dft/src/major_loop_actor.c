
#include "major_loop_actor.h"



void iterator(int ITER, int START, OUT int *cycle_out)
{
	for (int i = 0; i < ITER; i++)
	{
		cycle_out[i] = START + i;
	}
}

void sink_gains(__attribute__((unused)) int NUM_RECEIVERS, __attribute__((unused))  PRECISION2 *input)
{
	return;
}