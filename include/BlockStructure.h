/**********************************************************
 * Piotr T. Różański (c) 2015-2023                        *
 *   Enhanced Matching Pursuit Implementation (empi)      *
 * See README.md and LICENCE for details.                 *
 **********************************************************/
#ifndef EMPI_BLOCK_STRUCTURE_H
#define EMPI_BLOCK_STRUCTURE_H

struct BlockStructure {
	double scale;
	int envelope_length;
	int transform_size;
	double input_shift;
};

#endif //EMPI_BLOCK_STRUCTURE_H
