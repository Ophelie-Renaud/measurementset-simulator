//
// Created by hmiomand on 09/03/2020.
//


#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <limits.h>
#include <float.h>

#ifndef __NVCC__
#include <stdatomic.h>
#endif

#ifndef TEST_FIFOGET_FIFOFUNCTION_H
#define TEST_FIFOGET_FIFOFUNCTION_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef min
#define min(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

#ifndef max
#define max(a, b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#endif

// Probably not the best way to handle this
#if (__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__)
	#ifndef __NVCC__
		#include <x86intrin.h>
		#define BYTE_SWAP64(x)		_bswap64(x)
		#define BYTE_SWAP32(x)		_bswap(x)
		#define BYTE_SWAP16(x)		_rotwl(x, 8)
	#else
		#warning Using generic byteswap implementation.
		#define BYTE_SWAP64(x)		invert_uint64(x)
		#define BYTE_SWAP32(x)		invert_uint32(x)
		#define BYTE_SWAP16(x)		invert_uint16(x)
	#endif
#endif

#if 1
#define FIFO_GET_FUNC(fifo, index)					fifoGet_int32(fifo, index)
#define FIFO_SET_FUNC(fifo, data, index)			fifoSet_int32(fifo, data, index)
#else
#define FIFO_GET_FUNC(fifo, index)					fifoGet_int32_8bytes(fifo, index)
#define FIFO_SET_FUNC(fifo, data, index)			fifoSet_int32_8bytes(fifo, data, index)
#endif

#define FP32_EXPONENT_BIAS	(127)
#define FP32_MANTISSA_SIZE	(23)
#define SCALING_AMPLITUDE_CORRECTION	(1.01f)

#define SIGN_POS	(0)
#define SIGN_NEG	(1)
#define SIGN_UNDEF	(-1)

typedef enum FifoType {NO_TYPE, FxP, FxP_control, cFP, SCALING, SCALING_control} FifoType;

typedef struct StorageFifo {
	uint8_t* pointerToBuffer;
	int32_t storageSize;
	int32_t workingSize;
	uint32_t nb_element;
	FifoType fifoType;

	union {
		int32_t storageDecimalSize; // Use as mantissa size if custom float size ?
		int32_t mantissaSize;		// Might be used in a union with storageDecimalSize to save space
	};
	union {
		float stepSize;
		int32_t exponentSize;
	};
	union {
		float offset;
		int32_t exponentBias;
	};
	int32_t signBit;

	int32_t originalDecimalSize;
	float originalOffset;
	float originalStepSize;

	float maxValue;
	float minValue;

	float currentMaxValue;
	float currentMinValue;

	int32_t currentStorageIntegerSize;

	// For future purposes
//	uint32_t (*ptr_fifoGet_int32) (struct StorageFifo*, int);
//	uint32_t (*ptr_fifoSet_int32) (struct StorageFifo*, int);

	uint32_t indexOffset;
	uint32_t nbInstance;

} StorageFifo;

typedef struct StorageFifoVec {
	uint8_t* pointerToBuffer;
	int32_t storageSize;
	int32_t workingSize;
	uint32_t nb_element;
	FifoType fifoType;

	int32_t* storageDecimalSize; // Use as mantissa size if custom float size ?
	int32_t* mantissaSize;		// Might be used in a union with storageDecimalSize to save space
	int32_t* exponentSize;
	int32_t* exponentBias;
	int32_t* signBit;
//	uint8_t vectorSize;

	float* stepSize;
	float* offset;

	uint8_t nbComponent;
} StorageFifoVec;

int32_t getExponentValue(float value);

uint64_t invert_uint64(uint64_t data_in);
uint32_t invert_uint32(uint32_t data_in);
uint16_t invert_uint16(uint16_t data_in);

uint32_t fifoGet_int32(StorageFifo* fifo, int index);
void fifoSet_int32(StorageFifo* fifo, uint32_t data_in, int index);

uint32_t fifoGet_int32_8bytes(StorageFifo* fifo, int index);
void fifoSet_int32_8bytes(StorageFifo* fifo, uint32_t data_in, int index);

uint32_t fifoGet_int32_old(StorageFifo* fifo, int index);
void fifoSet_int32_old(StorageFifo* fifo, uint32_t data_in, int index);

uint32_t leftSideDataGet_8bits();

void fifoGet(StorageFifo* fifo, void* data_out, int index);
void fifoSet(StorageFifo* fifo, void* data_in, int index);
float fifoGet_float(StorageFifo* fifo, int index);
float fifoSet_float(StorageFifo* fifo, float data_in, int index);
StorageFifo fifoCreate(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize);
StorageFifo fifoCreate_multi(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, int nb_instance);
StorageFifo fifoCreate_slave(StorageFifo* master_AxB, int actor_instance);
void fifoDestroy(StorageFifo* fifo);
void fifoReset(StorageFifo* fifo);
void fifoSimulate_float(int NB_ELEMENT, float* buffer, int storageSize);

// FxP
StorageFifo fifoCreate_FxP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float maxValue);
float fifoGet_float_FxP(StorageFifo* fifo, int index);
void fifoSet_float_FxP(StorageFifo* fifo, float data_in, int index);
void fifoSimulate_float_FxP(int NB_ELEMENT, float* buffer, int storageSize);

void fifoSet_float_FxP_control(StorageFifo* fifo, float data_in, int index);
float fifoGet_float_FxP_override(StorageFifo* fifo, int index, uint32_t storageDecimalSize);
void fifoSet_float_FxP_override(StorageFifo* fifo, float data_in, int index, uint32_t storageDecimalSize);
void fifoReset_FxP(StorageFifo* fifo);

// Scaling
StorageFifo fifoCreate_scaling(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float minValue, float maxValue);
float fifoGet_float_scaling(StorageFifo* fifo, int index);
void fifoSet_float_scaling(StorageFifo* fifo, float data_in, int index);
void fifoSimulate_float_scaling(int NB_ELEMENT, float* buffer, int storageSize);

void fifoSet_float_scaling_control(StorageFifo* fifo, float data_in, int index);
float fifoGet_float_scaling_override(StorageFifo* fifo, int index, float stepSize, float offset);
void fifoSet_float_scaling_override(StorageFifo* fifo, float data_in, int index, float stepSize, float offset);
void fifoReset_scaling(StorageFifo* fifo);


// cFP
StorageFifo fifoCreate_cFP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float maxAbsExponent);
StorageFifo fifoCreate_cFP_better(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float minExponent, float maxExponent);
StorageFifo fifoCreate_cFP_better_sign(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float minExponent, float maxExponent, int sign);
float fifoGet_float_cFP(StorageFifo* fifo, int index);
void fifoSet_float_cFP(StorageFifo* fifo, float data_in, int index);
void fifoSimulate_float_cFP(int NB_ELEMENT, float* buffer, int storageSize);

void fifoSet_float_cFP_control(StorageFifo* fifo, float data_in, int index);

// VECTOR generic
float fifoVecGet_float(StorageFifoVec* fifo, int index, int component);
float fifoVecSet_float(StorageFifoVec* fifo, float data_in, int index, int component);

uint32_t fifoVecGet_int32(StorageFifoVec* fifo, int index, int component);
void fifoVecSet_int32(StorageFifoVec* fifo, uint32_t data_in, int index, int component);
void fifoVecSimulate_float(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent);

// cFP VECTOR
StorageFifoVec fifoVecCreate_cFP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* maxAbsExponent, int nbComponent);
StorageFifoVec fifoVecCreate_cFP_better(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* minExponent, float* maxExponent, int nbComponent);
StorageFifoVec fifoVecCreate_cFP_better_sign(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* minExponent, float* maxExponent, int* sign, int nbComponent);
float fifoVecGet_float_cFP(StorageFifoVec* fifo, int index, int component);
void fifoVecSet_float_cFP(StorageFifoVec* fifo, float data_in, int index, int component);
void fifoVecSimulate_float_cFP(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent);


StorageFifoVec fifoVecCreate(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, uint8_t nbComponent);
void fifoVecDestroy(StorageFifoVec* fifo);

// FxP VECTOR
float fifoVecGet_float_FxP(StorageFifoVec* fifo, int index, int component);
void fifoVecSet_float_FxP(StorageFifoVec* fifo, float data_in, int index, int component);
void fifoVecSet_float_FxP_control(StorageFifoVec* fifo, float data_in, int index, int component);

StorageFifoVec fifoVecCreate_FxP(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* maxValue, int nbComponent);
void fifoVecSimulate_float_FxP(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent);


// scaling VECTOR
StorageFifoVec fifoVecCreate_scaling(uint32_t nb_element, uint32_t storageSize, uint32_t workingSize, float* minValue, float* maxValue, int nbComponent);
float fifoVecGet_float_scaling(StorageFifoVec* fifo, int index, int component);
void fifoVecSet_float_scaling(StorageFifoVec* fifo, float data_in, int index, int component);
void fifoVecSimulate_float_scaling(int NB_ELEMENT, void* buffer, int storageSize, int nbComponent);

/*

<TYPE> fifoGet_<TYPE>(StorageFifo* fifo, int index);
void fifoSet_<TYPE>(StorageFifo* fifo, <TYPE> data_in, int index);

 */

#ifdef __cplusplus
}
#endif

#endif //TEST_FIFOGET_FIFOFUNCTION_H
